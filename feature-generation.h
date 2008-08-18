/* Copyright (C) 2005  Christoph Helma <helma@in-silico.de> 

   
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/


#include "feature-db.h"

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/groupcontrib.h>
#include <openbabel/obconversion.h>
#include <unistd.h>

using namespace std;
using namespace OpenBabel;


//! a class that knows several feature generation techniques
template <class MolType, class FeatureType, class ActivityType>
class FeatGen {

	public:
		typedef FeatMol < MolType, FeatureType, ActivityType > * MolRef ;
		typedef Feature<OBLinFrag> * OBLinFragRef;
		typedef Feature<FeatureType> * FeatRef;

	private:

		MolVect< MolType, FeatureType, ActivityType > * structures;
		vector<OBLinFragRef> alphabet;
		vector<OBLinFragRef> level;
		vector<OBLinFragRef> next_level;
		Out * out;

	public:

		FeatGen< MolType, FeatureType, ActivityType >(MolVect< MolType, FeatureType, ActivityType > * s, Out * out): structures(s), out(out) { };

		FeatGen< MolType, FeatureType, ActivityType >(char * alphabet_file, MolVect< MolType, FeatureType, ActivityType > * s, Out * out): structures(s), out(out) {
			*out << "Reading alphabet from " << alphabet_file << endl;
			out->print_err();
			this->read_smarts(alphabet_file,true,false);
		};

		// AM: from LOO 
		FeatGen< MolType, FeatureType, ActivityType >(char * alphabet_file, MolVect< MolType, FeatureType, ActivityType > * s, MolRef mol, Out * out): structures(s), out(out) {
			this->read_smarts(alphabet_file,false,false);
		};

		//! generate linear fragments from all structures
		void generate_linfrag();

		//! generate linear fragments that occur in test_mol and train_structures
		void generate_linfrag(FeatMolVect< MolType, FeatureType, ActivityType > * train_structures, MolRef test_mol);

		//! Prints testset for random selection
		void generate_testset(int p, Out* out);

		//! Prints testset for a fold, i.e. non-random consecutive.
		void generate_testset(int p, int pos, Out* out);

		//! generate physicochemical properties
		void generate_pcprop(Out* out);
		
		//! generate rex fragments
		void generate_rex(int max_l);

		//! read (and match) smarts from a file
		void read_smarts(char * smarts_file, bool print, bool split_bonds);

		//! generate smallest set of smallest rings
		void sssr();

		void refine_linfrag(unsigned int min_freq);

		void refine_rex(int level);

		void add_cand(string can_sma,LinFrag frag, vector<int> pot_matches);

		int level_size() { return( level.size() ); };

		void level_add(FeatRef feat_ptr) { level.push_back(feat_ptr); };

		vector<OBLinFragRef> get_level();

		void match_level(bool print);

		void match(OBLinFragRef feat_ptr);

		void match(MolRef test_comp);

		void copy_level(FeatMolVect<MolType, FeatureType, ActivityType> * f, MolRef test_comp);

		void alphabet2level();
		
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::generate_linfrag() {

	clock_t t;
	int l = 0;

	level = alphabet;

	while (level.size() > 0) {
		
		l++;
		*out << "LEVEL " << l << endl;
		out->print_err();

		t = clock();
		this->match_level(true);
		t = (clock() - t)/1000;
		*out << "Matching [ " << level.size() << " features, "<< t << " k ticks ]" << endl;
		out->print_err();

		t = clock();
		this->refine_linfrag(1);
		t = (clock() - t)/1000;
		*out << "Refinement [ "<< t << " k ticks ]" << endl;
		out->print_err();

	}
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::generate_linfrag(FeatMolVect< MolType, FeatureType, ActivityType > * train_structures, MolRef test_comp) {

	level = alphabet;

	while (this->level_size() > 0) {
		this->match(test_comp);
		this->match_level(false);
		this->copy_level(train_structures,test_comp);
		this->refine_linfrag(0);
	}
	
};

//! Prints testset for random selection
template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::generate_testset(int p, Out* out) {

	typename vector<MolRef>::iterator vmr_it;
        vector<MolRef> s = structures->get_compounds();
	vector<MolRef> drawn;
	MolRef m = NULL;
	double dpos = 0.0;

	float frac = (100.0 / (float) p);
	int times = (int) (s.size()/frac);
	// write validation file
	srand(static_cast<unsigned int>(clock()));
	for (int i=0; i<times; i++) {
		dpos = double(rand()) / (double(RAND_MAX) + 1.0);
		m = s[(int) (dpos*s.size())];
		// see if already drawn. if not, remember as drawn and output
		for (vmr_it = drawn.begin(); vmr_it != drawn.end(); vmr_it++) {
			if ((*vmr_it) == m) break;
		}
		if (vmr_it == drawn.end()) {
			drawn.push_back(m);
			*out << m->get_id() << "\t";
			*out << m->get_smiles() << endl;
		}
		else i--;
	}
	out->print();

};


//! Prints testset for a fold, i.e. non-random consecutive.
template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::generate_testset(int p, int pos, Out* out) {

        vector<MolRef> s = structures->get_compounds();
	MolRef m = NULL;
	
	float frac = (100.0 / (float) p);
	int calculated_set_size = (int) (s.size()/frac);
	
	// write validation file
	for (int i=0; i<calculated_set_size; i++) {
		try {
			m = s.at(pos);
		}
		catch (...) {
			cerr << "End of structures reached after " << i << " of " << calculated_set_size << " structures." << endl;
			break;
		}
		// see if already drawn. if not, remember as drawn and output
		*out << m->get_id() << "\t";
		*out << m->get_smiles() << endl;
		pos++;
	}
	out->print();
};


//! Prints physicochemical properties in the form id \t molweight \t logP \t psa \t mr
template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::generate_pcprop(Out* out) {

	typename vector<MolRef>::iterator structure_it;
	vector<MolRef> s = structures->get_compounds();
	
	for (structure_it = s.begin(); structure_it!=s.end(); structure_it++) {
		*out << (*structure_it)->get_id() << "\t";
		OBMol* mol;
		mol = (*structure_it)->get_mol_ref();
		*out << mol->GetMolWt() << "\t";
		OBLogP logP;
		OBPSA psa;
		OBMR mr;
		*out << logP.Predict(*mol) << "\t";
		*out << psa.Predict(*mol) << "\t";
		*out << mr.Predict(*mol) << "\t" << endl;
		out->print();
	}

};


template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::generate_rex(int max_l) {

	int l = 0;

	level = alphabet;
	this->match_level(true);
	alphabet = level;
	
	while (l < max_l) {
		
		l++;
		*out << "LEVEL " << l << endl;
		out->print_err();

		this->refine_rex(l);
		this->match_level(true);

	}
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::read_smarts(char * file, bool print, bool split_bonds) {

	ifstream input;
	input.open(file);
	if (!input) {
		*out << "Cannot open " << file << endl;
		out->print();
		exit(1);
	}
	
	string smarts;
	string line;
	while ( getline(input, line) ) {
		istringstream iss(line); 
		int i = 0;
		while(getline(iss, smarts, '\t') && i == 1) {	// split at tab and read only the first field (remove comments after the tab)
			i++;
		}
		Feature<OBLinFrag> * new_feat_ptr;
		new_feat_ptr = new Feature<OBLinFrag>(smarts, split_bonds);
		alphabet.push_back(new_feat_ptr);
	}
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::sssr() {

	vector<string> unique_smarts;
	vector<string> ringsmarts;
	vector<string> diff;

	typename vector<MolRef>::iterator cur_mol;
	typename vector<string>::iterator cur_sma;

	vector<MolRef> compounds = structures->get_compounds();
	for (cur_mol=compounds.begin();cur_mol!=compounds.end();cur_mol++) {

		ringsmarts.clear();
		diff.clear();

		ringsmarts = (*cur_mol)->sssr();

		sort(unique_smarts.begin(),unique_smarts.end());
		sort(ringsmarts.begin(),ringsmarts.end());

		set_difference(ringsmarts.begin(),ringsmarts.end(),
			unique_smarts.begin(),unique_smarts.end(),
			insert_iterator<vector<string> >(diff,diff.begin()));

		for (vector<string>::iterator cur_sma = diff.begin(); cur_sma != diff.end(); cur_sma++) {
			Feature<OBLinFrag> sma(*cur_sma,0);
			this->match(&sma);
			sma.print_matches(out);
			unique_smarts.push_back(*cur_sma);
		}

	}
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::refine_linfrag(unsigned int min_freq) {

	vector<int> * f1_matches;
	vector<int> * f2_matches;
	vector<int> pot_matches;
	typename vector<OBLinFragRef>::iterator frag1;
	typename vector<OBLinFragRef>::iterator frag2;
	typename vector<OBLinFragRef>::iterator cur_frag;

	for (frag1 = level.begin(); frag1 != level.end(); frag1++) {
		
		f1_matches = (*frag1)->get_matches_ptr();
		
		for (frag2 = frag1; frag2 != level.end(); frag2++) {

			f2_matches = (*frag2)->get_matches_ptr();
			pot_matches.clear();
			set_intersection(f1_matches->begin(),f1_matches->end(),
					f2_matches->begin(),f2_matches->end(),
					insert_iterator<vector<int> >(pot_matches,pot_matches.begin()));

			// if we have only elements (first level), we have to add the bonds
			if ((*frag1)->size() == 1) {
				string bonds[] = {"-","=","#",":"};
				for (int n = 0; n < 4; n++) {
					LinFrag newfrag ( (*frag2)->get_name(), 0);
					newfrag.expand(bonds[n]);
					newfrag.expand((*frag1)->get_name());
					string can_sma = newfrag.canonify();
					this->add_cand(can_sma, newfrag, pot_matches);
				}
			}

			else if ((pot_matches.size() > 0)&&(f1_matches->size()>min_freq)&&(f2_matches->size()>min_freq)) { // generate only the most general fragments with frequency==min_freq

				vector<string> f1_top = (*frag1)->get_fragment();
				f1_top.erase(f1_top.begin(),f1_top.begin()+2);

				vector<string> f1_rev_top = (*frag1)->get_fragment();
				reverse(f1_rev_top.begin(),f1_rev_top.end());
				f1_rev_top.erase(f1_rev_top.begin(),f1_rev_top.begin()+2);

				vector<string> f2_bottom = (*frag2)->get_fragment();
				f2_bottom.erase(f2_bottom.end()-2,f2_bottom.end());

				vector<string> f2_rev_bottom = (*frag2)->get_fragment();
				reverse(f2_rev_bottom.begin(),f2_rev_bottom.end());
				f2_rev_bottom.erase(f2_rev_bottom.end()-2,f2_rev_bottom.end());
		
				if (f1_top == f2_bottom) {
					LinFrag newfrag = **frag2;
					newfrag.rev();
					newfrag.expand((*frag1)->first_bond());
					newfrag.expand((*frag1)->first_atom());
					string can_sma = newfrag.canonify();
					this->add_cand(can_sma, newfrag, pot_matches);
				}

				if (f1_top == f2_rev_bottom) {
					LinFrag newfrag = **frag2;
					newfrag.expand((*frag1)->first_bond());
					newfrag.expand((*frag1)->first_atom());
					string can_sma = newfrag.canonify();
					this->add_cand(can_sma, newfrag, pot_matches);
				}

				if (f1_rev_top == f2_bottom) {
					LinFrag newfrag = **frag2;
					newfrag.rev();
					newfrag.expand((*frag1)->last_bond());
					newfrag.expand((*frag1)->last_atom());
					string can_sma = newfrag.canonify();
					this->add_cand(can_sma, newfrag, pot_matches);
				}

				if (f1_rev_top == f2_rev_bottom) {
					LinFrag newfrag = **frag2;
					newfrag.expand((*frag1)->last_bond());
					newfrag.expand((*frag1)->last_atom());
					string can_sma = newfrag.canonify();
					this->add_cand(can_sma, newfrag, pot_matches);
				}
			}
		}
		
		// delete fragments from previous level
		delete *frag1;
		*frag1 = NULL;
	}


	level.clear();
	for (cur_frag=next_level.begin(); cur_frag != next_level.end(); cur_frag++) {
		level.push_back(*cur_frag);
	}
	next_level.clear();
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::refine_rex(int l) {

	vector<int> * f1_matches;
	vector<int> * f2_matches;
	vector<int> pot_matches;
	typename vector<OBLinFragRef>::iterator frag1;
	typename vector<OBLinFragRef>::iterator frag2;
	typename vector<OBLinFragRef>::iterator cur_frag;

	for (frag1 = alphabet.begin(); frag1 != alphabet.end(); frag1++) {
		
		f1_matches = (*frag1)->get_matches_ptr();
		
		for (frag2 = frag1; frag2 != alphabet.end(); frag2++) {

			f2_matches = (*frag2)->get_matches_ptr();
			pot_matches.clear();
			set_intersection(f1_matches->begin(),f1_matches->end(),
					f2_matches->begin(),f2_matches->end(),
					insert_iterator<vector<int> >(pot_matches,pot_matches.begin()));

			int n = 1;
			LinFrag newfrag ( (*frag2)->get_name(), 0);
			while (n < l) {
					newfrag.expand("~");
					newfrag.expand("*");
					n++;
			}
			newfrag.expand("~");
			newfrag.expand((*frag1)->get_name());
			string can_sma = newfrag.canonify();
			this->add_cand(can_sma, newfrag, pot_matches);

		}
	}


	level.clear();
	for (cur_frag=next_level.begin(); cur_frag != next_level.end(); cur_frag++) {
		level.push_back(*cur_frag);
	}
	next_level.clear();
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::match(MolRef test_comp) {

	typename vector<OBLinFragRef>::iterator frag;
	for (frag = level.begin(); frag != level.end(); frag++) {

		if (!test_comp->match(*frag)) {
			delete *frag; // free memory
			*frag = NULL;
			level.erase(frag);
			frag--;	// set the iterator back to avoid skipping the next compound
		}
	}
};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::match(OBLinFragRef feat_ptr) {

	vector<MolRef> compounds = structures->get_compounds();
	typename vector<MolRef>::iterator cur_mol;
	int comp_nr = 0;
	OBSmartsPattern * sp = feat_ptr->get_smarts_pattern();
	
	for (cur_mol = compounds.begin(); cur_mol != compounds.end(); cur_mol++) {

		if ( sp->Match(*(*cur_mol)->get_mol_ref(),true) ) {
			feat_ptr->add_match(comp_nr);
		}

		comp_nr++;

	}

};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::match_level(bool print) {

	typename vector<OBLinFragRef>::iterator frag;

	for (frag = level.begin(); frag != level.end(); frag++) {

		this->match(*frag);

		if ((*frag)->nr_matches() == 0) {
			delete *frag; // free memory
			*frag = NULL;
			level.erase(frag);
			frag--;	// set the iterator back to avoid skipping the next compound
		}
		
		else if (print) {
			(*frag)->print_matches(out);
		}
		
	}

};

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::add_cand(string can_sma, LinFrag newfrag, vector<int> pot_matches) {

	bool isnew = true;
	typename vector<OBLinFragRef>::iterator cur_frag;

	for (cur_frag = next_level.begin(); cur_frag != next_level.end(); cur_frag++) {
		if ((*cur_frag)->get_name() == can_sma) {
			isnew = false;
			break;
		}
	}

	if (isnew) {
		OBLinFragRef frag_ptr;
		frag_ptr = new Feature<OBLinFrag>(can_sma, newfrag, pot_matches);
		next_level.push_back(frag_ptr);
	}

};
			
template <class MolType, class FeatureType, class ActivityType>
vector< Feature< OBLinFrag> * >  FeatGen<MolType, FeatureType, ActivityType>::get_level() {
	return(level);
}

template <class MolType, class FeatureType, class ActivityType>
void FeatGen<MolType, FeatureType, ActivityType>::copy_level(FeatMolVect<MolType, FeatureType, ActivityType> * f, MolRef test_comp) {

	typename vector<OBLinFragRef>::iterator frag;
	typename map<const string, Feature<FeatureType> *>::iterator pos;
	string name;
	
	for (frag = level.begin(); frag != level.end(); frag++) {

		name = (*frag)->get_name();

		if ((*frag)->nr_matches() == 0) {
			test_comp->add_unknown(name);
			delete *frag; // free memory
			*frag = NULL;
			level.erase(frag);
			frag--;	// set the iterator back to avoid skipping the next compound
		}

		else 	// search in the feature_map
			f->add_feature(test_comp,name);

	}

};

template <class MolType, class FeatureType, class ActivityType>
void  FeatGen<MolType, FeatureType, ActivityType>::alphabet2level() {
	level = alphabet;
}
