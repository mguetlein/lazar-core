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


#include "lazmol.h"
#include <string>

#ifndef LAZMOLVECT_H
#define LAZMOLVECT_H

using namespace std;
using namespace OpenBabel;

extern bool quantitative;

//! container for LazMol objects
template <class MolType, class FeatureType, class ActivityType>
class MolVect {

	public:

		typedef FeatMol < MolType, FeatureType, ActivityType > * MolRef ;

	private:

		vector<MolRef> compounds;
		Out * out;

	public:

		MolVect() {};

		//! MolVect constructor: reads SMILES from file (called by FeatMolVect())
		MolVect(char * structure_file, Out * out);

		//! add a new feature to compound comp_nr
		void add_feature(int comp_nr, Feature<FeatureType> * feat_ptr) {
			compounds[comp_nr]->add_feature(feat_ptr);
		};

		//! determine common features with the training compounds. Results are used to find neighbors later.
		void common_features(MolRef test_compound);

		//! Determine similarity as weighted tanimoto index.
		void relevant_features(MolRef test, string act);

		//! determine unknown features
		void determine_unknown(string act, MolRef test);

		//! remove duplicates of the query structure
		vector<MolRef> remove_duplicates(MolRef test_comp);

		//! Get Neighbors by using at least five compounds if any neighbors available
		vector<MolRef> get_neighbors(string act);

		vector<MolRef> get_compounds() { return(compounds); };

		MolRef get_compound(int n) { return(compounds[n]); };

		int get_size() { return(compounds.size()); };

		MolRef get_molfromid(string id);
		vector<string> get_idfromsmi(string smi);
		vector<string> get_idfrominchi(string inchi);


};

template <class MolType, class FeatureType, class ActivityType>
MolVect<MolType, FeatureType, ActivityType>::MolVect(char * structure_file, Out * out): out(out) {

	string line;
	string tmp_field;
	string id;
	string smi;
	string inchi;
	vector<string> ids, dup_ids;
	vector<string>::iterator dup_id;
	vector<string> inchis;
	vector<string>::iterator dup_inchi;

	MolRef mol_ptr;
	int line_nr = 0;

	ifstream input;
	input.open(structure_file);

	if (!input) {
		*out << "Cannot open " << structure_file << endl;
		out->print_err();
		exit(1);
	}

	*out << "Reading structures from " << structure_file << endl;
	out->print_err();

	line_nr = 0;
	int oldnr = 2;

	while (getline(input, line)) {

		istringstream iss(line);

		int field_nr = 0;

		while(getline(iss, tmp_field, '\t')) {	// split at tabs

			if (field_nr == 0) {		// ID

				id = tmp_field;
				dup_id = find(ids.begin(), ids.end(), id);

				if (dup_id == ids.end())
					ids.push_back(id);

				else {
					*out << id << " (line " << line_nr << ") is not a unique ID ... exiting.\n";
					out->print_err();
					exit(1);
				}


			}

			else if (field_nr == 1)	{ // SMILES
				smi = tmp_field;
				remove_dos_cr(&smi);

				mol_ptr = new FeatMol<MolType,FeatureType,ActivityType>(line_nr, id, smi,out);

				inchi = mol_ptr->get_inchi();
				dup_inchi = find(inchis.begin(), inchis.end(), inchi);

				if (dup_inchi == inchis.end())
					inchis.push_back(inchi);

				else {
					dup_ids =  this->get_idfrominchi(inchi);
					*out << "Compounds " << id ;
					for (dup_id=dup_ids.begin();dup_id!=dup_ids.end();dup_id++) {
						*out << " and " << *dup_id;
					}
					*out << " have identical structures.\n";
					out->print_err();
				}

			}


			oldnr = field_nr;
			field_nr++;
		}

		compounds.push_back(mol_ptr);
		line_nr++;

	}

	input.close();

};


template <class MolType, class FeatureType, class ActivityType>
vector<FeatMol < MolType, FeatureType, ActivityType > * > MolVect<MolType, FeatureType, ActivityType>::remove_duplicates(MolRef test_comp) {

	vector<MolRef> duplicates;
	typename vector<MolRef>::iterator cur_mol;

	test_comp->clear_db_activities();

	for (cur_mol=compounds.begin();cur_mol!=compounds.end();cur_mol++) {
		if ((*cur_mol)->equal(test_comp)) {
			(*cur_mol)->copy_activities(test_comp); // copy database activities to test structure (for external predictions)
			duplicates.push_back(*cur_mol);
			(*cur_mol)->remove();
		}
	}
	return(duplicates);
};

template <class MolType, class FeatureType, class ActivityType>
void MolVect<MolType, FeatureType, ActivityType>::common_features(MolRef test) {

	typename vector<MolRef>::iterator cur_mol;

	for (cur_mol=compounds.begin();cur_mol!=compounds.end();cur_mol++) {
		(*cur_mol)->common_features(test);
	}

};

template <class MolType, class FeatureType, class ActivityType>
void MolVect<MolType, FeatureType, ActivityType>::relevant_features(MolRef test, string act) {

	typename vector<MolRef>::iterator cur_mol;

	for (cur_mol=compounds.begin();cur_mol!=compounds.end();cur_mol++) {
		(*cur_mol)->relevant_features(test, act);
	}

};

template <class MolType, class FeatureType, class ActivityType>
vector<FeatMol < MolType, FeatureType, ActivityType > * >  MolVect<MolType, FeatureType, ActivityType>::get_neighbors(string act) {

	vector<MolRef> neighbors;
	typename vector<MolRef>::iterator cur_mol;
	neighbors.clear();


    multimap<float,MolRef> sim_sorted_neighbors;
    sim_sorted_neighbors.clear();
    for (cur_mol=compounds.begin();cur_mol!=compounds.end();cur_mol++) {
        if ((*cur_mol)->is_available(act)) {
            float sim = (*cur_mol)->get_similarity();
            sim_sorted_neighbors.insert(pair<float, MolRef>(sim,(*cur_mol)));
        }
    }
    typename multimap<float,MolRef>::iterator cur_sn;

    // cutoff 0.3 ~ (1/100)^(1/4), i.e. 100 compounds of similarity 0.3 are needed to compensate 1 compound of sim
    cur_sn = sim_sorted_neighbors.end();
    cur_sn--;
    while ((cur_sn != sim_sorted_neighbors.begin()) && (cur_sn->second->get_similarity()>0.3)) {
        if (cur_sn->second->is_available(act)) {
            neighbors.push_back(cur_sn->second);
        }
        cur_sn--;
    }

    // * Comment this out to revert to original Lazar implementation *
    unsigned int min_n = 5;
    if (!quantitative) min_n=5;

    if ((cur_sn != sim_sorted_neighbors.begin()) && (neighbors.size() < min_n) && (neighbors.size() > 0)) {
        do {
            if (cur_sn->second->is_available(act)) {
                neighbors.push_back(cur_sn->second);
            }
            cur_sn--;
        } while ((neighbors.size() < min_n) && (cur_sn != sim_sorted_neighbors.begin()));
    }
    // *                                                             *

	return(neighbors);
};


template <class MolType, class FeatureType, class ActivityType>
void MolVect<MolType, FeatureType, ActivityType>::determine_unknown(string act, MolRef test) {

			vector<Feature<FeatureType> *> feats = test->get_features();
			typename vector<Feature<FeatureType> *>::iterator cur_feat;
			bool act_m;
			vector<int> matches;
			vector<int>::iterator cur_m;

			test->delete_unknown();
			//test->delete_infrequent();

			for (cur_feat=feats.begin();cur_feat!=feats.end();cur_feat++) {

				// find features without activity values for the current activity
				act_m = false;
				matches = (*cur_feat)->get_matches();
				for (cur_m=matches.begin();cur_m!=matches.end();cur_m++) {
					if (compounds[*cur_m]->is_available(act)) {
						act_m = true;
						break;
					}
				}

				if ( !act_m || (*cur_feat)->get_too_infrequent(act) )
					test->add_unknown((*cur_feat)->get_name());

			}

		};

template <class MolType, class FeatureType, class ActivityType>
vector<string>  MolVect<MolType, FeatureType, ActivityType>::get_idfromsmi(string smi) {

	typename vector<MolRef>::iterator cur_mol;
	vector<string> ids;

	for (cur_mol=compounds.begin(); cur_mol!= compounds.end(); cur_mol++) {
		if ((*cur_mol)->get_smiles()==smi)
			ids.push_back((*cur_mol)->get_id());
	}
	return(ids);
};

template <class MolType, class FeatureType, class ActivityType>
vector<string>  MolVect<MolType, FeatureType, ActivityType>::get_idfrominchi(string inchi) {

	typename vector<MolRef>::iterator cur_mol;
	vector<string> ids;

	for (cur_mol=compounds.begin(); cur_mol!= compounds.end(); cur_mol++) {
		if ((*cur_mol)->get_inchi()==inchi)
			ids.push_back((*cur_mol)->get_id());
	}
	return(ids);
};

template <class MolType, class FeatureType, class ActivityType>
FeatMol < MolType, FeatureType, ActivityType > *  MolVect<MolType, FeatureType, ActivityType>::get_molfromid(string id) {

	typename vector<MolRef>::iterator cur_mol;

	for (cur_mol=compounds.begin(); cur_mol!= compounds.end(); cur_mol++) {
		if ((*cur_mol)->get_id()==id)
			break;
	}

	if (cur_mol == compounds.end())
		return(NULL);
	else
		return(*cur_mol);

};

#endif // LAZMOLVECT_H
