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


#include "lazmolvect.h"

#ifndef F_DB_H
#define F_DB_H

using namespace std;

//! compounds with features
template <class MolType, class FeatureType, class ActivityType>
class FeatMolVect: public MolVect<MolType,FeatureType,ActivityType> {

	public:

		typedef FeatMol < MolType, FeatureType, ActivityType > * MolRef ;
		typedef Feature<FeatureType> * FeatRef;

	private:

		map<const string, FeatRef> feature_map;		// lookup features by name
		vector<FeatRef> features;
		Out * out;

	public:

		FeatMolVect< MolType, FeatureType, ActivityType >(char * feat_file, char * structure_file, Out * out); 
        ~FeatMolVect() {
            for (unsigned int i=0; i<features.size(); i++) {
                delete features[i];
            }
        }

		//! read features from a file
		void read_features(char * feat_file);

		void add_feature(MolRef s, string name);

		void copy_level(vector<Feature<FeatureType> *> * level, MolRef test_comp);

		vector<FeatRef> * get_features() { return(&features); };
};

// read a feature file
template <class MolType, class FeatureType, class ActivityType>
FeatMolVect<MolType, FeatureType, ActivityType>::FeatMolVect(char * feat_file, char * structure_file, Out * out): MolVect< MolType, FeatureType, ActivityType >(structure_file,out), out(out) {

	ifstream input;
	input.open(feat_file);

	if (!input) {
		*out << "Cannot open " << feat_file << endl;
		out->print_err();
		exit(1);
	}

	string line;
	string tmp_field;
	int line_nr =1;
	FeatRef feat_ptr;

	*out << "Reading features from " << feat_file << endl;
	out->print_err();
	while (getline(input, line)) {
		
		istringstream iss(line); 
		int i =0;

		while(getline(iss, tmp_field, '\t') && i < 2) {	// split at tab and read exactly 2 fields
			remove_dos_cr(&tmp_field);

			if (i==0) {		// SMARTS

				feat_ptr = new Feature<FeatureType>(tmp_field); // initialize Feature with smarts 
				feature_map[tmp_field] = feat_ptr;
				features.push_back(feat_ptr);

			}

			else {			// MATCHES

				size_t startpos = 0, endpos = 0;
				vector<string> tokens;

				for (;;) {

					startpos = tmp_field.find_first_not_of(" \t\n",startpos);
					endpos   = tmp_field.find_first_of(" \t\n",startpos);

					if (endpos < tmp_field.size() && startpos <= tmp_field.size())
						tokens.push_back(tmp_field.substr(startpos,endpos-startpos));
					else
						break;

					startpos = endpos + 1;

				}

				for (unsigned int i = 0 ; i < tokens.size() ; i++ ) {

					if ( tokens[i] == "[" )
						continue;
					else if ( tokens[i] == "]" )
						break;
					  
					int comp_nr = atoi(tokens[i].c_str());
					feat_ptr->add_match(comp_nr);	// comp_nr is line number
					this->get_compound(comp_nr)->add_feature(feat_ptr);

				}

			}

			i++;

		}

		line_nr++;

	}

	input.close();

};

template <class MolType, class FeatureType, class ActivityType>
void FeatMolVect<MolType, FeatureType, ActivityType>::add_feature(MolRef s, string name) {

	typename map<const string, FeatRef>::iterator pos;
	pos = feature_map.find(name);

	if (pos != feature_map.end())
		s->add_feature(pos->second);
};

template <class MolType, class FeatureType, class ActivityType>
void FeatMolVect<MolType, FeatureType, ActivityType>::copy_level(vector<Feature<FeatureType> *> * level, MolRef test_comp) {

	typename vector<Feature<FeatureType> *>::iterator frag;
	typename map<const string, Feature<FeatureType> *>::iterator pos;
	string name;
	
	for (frag = level->begin(); frag != level->end(); frag++) {

		name = (*frag)->get_name();

		if ((*frag)->nr_matches() == 0) {
			test_comp->add_unknown(name);
			delete *frag; // free memory
			*frag = NULL;
			level->erase(frag);
			frag--;	// set the iterator back to avoid skipping the next compound
		}

		else 	// search in the feature_map
			this->add_feature(test_comp,name);

	}

};
#endif // F_DB_H
