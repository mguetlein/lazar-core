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

using namespace std;

extern bool quantitative;

//! compounds with activities and features
template <class MolType, class FeatureType, class ActivityType>
class ActMolVect: public FeatMolVect< MolType, FeatureType, ActivityType > {

	private:

		vector<string> activity_names;
		Out * out;

	public:

		typedef FeatMol < MolType, FeatureType, ActivityType > * MolRef ;
		typedef Feature<FeatureType> * FeatRef;

		//! ActMolVect constructor, called directly by Predictor(). Reads in activity values after (implicitly) calling super class constructor FeatMolVect()
		ActMolVect< MolType, FeatureType, ActivityType >(char * act_file,char * feat_file, char  * structure_file, Out * out);

		//! read activities
		void read_act(char * act_file);

		//! get activity values for activity act
		vector<ActivityType> get_activity_values(string act);

		//! get activity values for a subset of structures and activity act
		vector<ActivityType> get_activity_values(vector<int> compound_numbers, string act);

		vector<string> get_activity_names() { return(activity_names); }

		// MG : precompute significance
		void precompute_feature_significance(string act, vector<bool> activity_values);
		void precompute_feature_significance(string act, vector<float> activity_values);
		// MG

		//! determine significance of boolean training set features using chi-sq test
		void feature_significance(string act, vector<bool> activity_values);

		//! determine significance of quantitative training set features using KS test
		void feature_significance(string act, vector<float> activity_values);

		void print_sig_features(float limit, char* smarts);
		void print_sorted_features(float limit, char* smarts);

};


// read activity file
template <class MolType, class FeatureType, class ActivityType>
ActMolVect<MolType, FeatureType, ActivityType>::ActMolVect(char* act_file,char* feat_file, char* structure_file, Out* out): FeatMolVect< MolType, FeatureType, ActivityType >(feat_file,structure_file,out), out(out) {

	string line;
	string tmp_field;
	string id;
	string no_id;
	string act_name;
	MolRef mol_ptr = NULL;
	int line_nr = 0;

	ifstream input;
	input.open(act_file);

	*out << "Reading activities from " << act_file << endl;
	out->print_err();

	if (!input) {
		*out << "Cannot open " << act_file << endl;
		out->print_err();
		exit(1);
	}

	line_nr = 0;

	while (getline(input, line)) {

		istringstream iss(line);
		int field_nr = 0;

		while(getline(iss, tmp_field, '\t')) {	// split at tabs

			remove_dos_cr(&tmp_field);

			if (field_nr == 0) {		// ID

				if (tmp_field == no_id)		// ignore compounds without structures
					break;
				else {
					mol_ptr = this->get_molfromid(tmp_field);
					if (mol_ptr == NULL) {	// ignore compounds without structures
						no_id = tmp_field;
						*out << "No structure for ID " << tmp_field << ".\n";
						out->print_err();
						break;
					}
				}
			}

			else if (field_nr == 1) {	// ACTIVITY NAME
				activity_names.push_back(tmp_field);
				act_name = tmp_field;
			}

			else if (field_nr == 2) {	// ACTIVITY VALUES

				if (tmp_field == "NA") {
					mol_ptr->set_na(act_name);
				}

				else {
					stringstream str;
					str  << tmp_field;
					ActivityType act_value;
					str >> act_value;

					if (quantitative) act_value = log10(act_value);
					mol_ptr->set_activity(act_name,act_value);
				}

			}

			else {
				*out << "More than 3 columns at line " << line_nr << " ... exiting.\n";
				out->print_err();
				exit(1);
			}

			field_nr++;

		}

		line_nr++;

	}

	input.close();

	// unique activity names
	sort(activity_names.begin(),activity_names.end());
	vector<string>::iterator uend = unique(activity_names.begin(),activity_names.end());
	activity_names.erase(uend,activity_names.end());

};



template <class MolType, class FeatureType, class ActivityType>
vector<ActivityType> ActMolVect<MolType, FeatureType, ActivityType>::get_activity_values(vector<int> comp_nrs, string act) {

	vector<ActivityType> activities;
	vector<ActivityType> tmp;
	vector<MolRef> compounds = this->get_compounds();
	vector<int>::iterator cur_comp;
	typename vector<ActivityType>::iterator cur_a;

	for (cur_comp=comp_nrs.begin();cur_comp!=comp_nrs.end();cur_comp++) {
		if (compounds[*cur_comp]->is_available(act)) {
			tmp = compounds[*cur_comp]->get_act(act);
			for (cur_a=tmp.begin();cur_a!=tmp.end();cur_a++) {
				activities.push_back(*cur_a);
			}
		}
	}

	return(activities);
};

template <class MolType, class FeatureType, class ActivityType>
vector<ActivityType> ActMolVect<MolType, FeatureType, ActivityType>::get_activity_values(string act) {

	vector<ActivityType> activities;
	vector<ActivityType> tmp;
	vector<MolRef> compounds = this->get_compounds();
	typename vector<ActivityType>::iterator cur_a;
	typename vector<MolRef>::iterator cur_comp;

	for (cur_comp=compounds.begin();cur_comp!=compounds.end();cur_comp++) {

		if ((*cur_comp)->is_available(act)) {
			tmp = (*cur_comp)->get_act(act);
			for (cur_a=tmp.begin();cur_a!=tmp.end();cur_a++) {
				activities.push_back(*cur_a);
			}
		}

	}

	return(activities);
};

template <class MolType, class FeatureType, class ActivityType>
void ActMolVect<MolType, FeatureType, ActivityType>::print_sig_features(float limit, char* smarts) {

	vector<ActivityType> activity_values;
	vector<string>::iterator cur_act;
	typename vector<Feature<FeatureType> * >::iterator cur_feat;
	vector<FeatRef> printed_features;
	vector<FeatRef> * features = this->get_features();;

	for (cur_act = activity_names.begin(); cur_act != activity_names.end(); cur_act++) {

		activity_values = this->get_activity_values(*cur_act);
		this->feature_significance(*cur_act, activity_values);

		for (cur_feat=features->begin(); cur_feat!=features->end(); cur_feat++) {

			if ( (*cur_feat)->get_p(*cur_act) > limit ) {

				if (find(printed_features.begin(),printed_features.end(),*cur_feat) == printed_features.end()) {
					(*cur_feat)->print_matches(out, smarts);
					printed_features.push_back(*cur_feat);
				}
			}
		}
	}
};

template <class MolType, class FeatureType, class ActivityType>
void  ActMolVect<MolType, FeatureType, ActivityType>::print_sorted_features(float limit, char* smarts) {

    vector<ActivityType> activity_values;
    vector<string>::iterator cur_act;
    typename vector<Feature<FeatureType> * >::iterator cur_feat;
    vector<FeatRef> printed_features;
    vector<FeatRef> * features = this->get_features();

    for (cur_act = activity_names.begin(); cur_act != activity_names.end(); cur_act++) {

        activity_values = this->get_activity_values(*cur_act);
        this->feature_significance(*cur_act, activity_values);

        for (cur_feat=features->begin(); cur_feat!=features->end(); cur_feat++) {

            if ( (*cur_feat)->get_p(*cur_act) > limit )
		(*cur_feat)->print_all(*cur_act, out, smarts);

        }
    }
};

template <class MolType, class FeatureType, class ActivityType>
void ActMolVect<MolType, FeatureType, ActivityType>::precompute_feature_significance(string act, vector<bool> activity_values) {

	int n_a =0;
	int n_i =0;
	vector<bool>::iterator cur_act_val;
	typename vector<FeatRef>::iterator cur_feat;
	vector<FeatRef> * features = this->get_features();

	// determine global nr of actives/inactives // AM: column sums
	for (cur_act_val=activity_values.begin();cur_act_val!=activity_values.end();cur_act_val++) {

		if (*cur_act_val)
			n_a++;
		else
			n_i++;

	}

	//int count = 0;

	// determine significance of training set features
	for (cur_feat=features->begin(); cur_feat!=features->end(); cur_feat++) {

		//if (count++ % 250 == 0)
		//	printf("precomputing significance %d/%d\n",count,features->size());

		if ((*cur_feat)->nr_matches() > 1) { // remove features that match on a single compound
			activity_values = this->get_activity_values((*cur_feat)->get_matches(), act);

			int f_a=0;
			int f_i=0;
			// AM: cell sums
			for (cur_act_val = activity_values.begin(); cur_act_val != activity_values.end(); cur_act_val++) {

				if (*cur_act_val)
					f_a++;
				else if (!*cur_act_val)
					f_i++;

			}

			(*cur_feat)->precompute_significance(act, n_a, n_i, f_a, f_i);
		}
	}

};

template <class MolType, class FeatureType, class ActivityType>
void ActMolVect<MolType, FeatureType, ActivityType>::precompute_feature_significance(string act, vector<float> activity_values){

	fprintf(stderr, "Not implemented for regression");
	exit(1);
}


template <class MolType, class FeatureType, class ActivityType>
void ActMolVect<MolType, FeatureType, ActivityType>::feature_significance(string act, vector<bool> activity_values) {

	int n_a =0;
	int n_i =0;
	vector<bool>::iterator cur_act_val;
	typename vector<FeatRef>::iterator cur_feat;
	vector<FeatRef> * features = this->get_features();

	// determine global nr of actives/inactives // AM: column sums
	for (cur_act_val=activity_values.begin();cur_act_val!=activity_values.end();cur_act_val++) {

		if (*cur_act_val)
			n_a++;
		else
			n_i++;

	}

	// determine significance of training set features
	for (cur_feat=features->begin(); cur_feat!=features->end(); cur_feat++) {

		if ((*cur_feat)->nr_matches() > 1) { // remove features that match on a single compound
			activity_values = this->get_activity_values((*cur_feat)->get_matches(), act);
			(*cur_feat)->determine_significance(act, n_a, n_i, &activity_values);	// AM: determine significance
		}
	}

};

template <class MolType, class FeatureType, class ActivityType>
void ActMolVect<MolType, FeatureType, ActivityType>::feature_significance(string act, vector<float> all_activity_values) {

	//float global_median;
	//vector<bool>::iterator cur_act_val;
	vector<float> feat_activity_values;
	typename vector<FeatRef>::iterator cur_feat;
	vector<FeatRef> * features = this->get_features();

	// determine significance of training set features
	for (cur_feat=features->begin(); cur_feat!=features->end(); cur_feat++) {
		if ((*cur_feat)->nr_matches() > 1) { // remove features that match on a single compound
			feat_activity_values = this->get_activity_values((*cur_feat)->get_matches(), act);

			/*
			sort(feat_activity_values.begin(), feat_activity_values.end());
			sort(all_activity_values.begin(), all_activity_values.end());

			// remove activity values that come from the current feature
			typename vector<float>::iterator cur_feat_act;
			typename vector<float>::iterator cur_all_act;
			cur_feat_act = feat_activity_values.begin();
			cur_all_act = all_activity_values.begin();


			vector<float> missing_activity_values;						// place to store missing activities

			while (cur_feat_act != feat_activity_values.end()) {
				if (cur_all_act == all_activity_values.end()) {
					break;
				}

				if ((*cur_feat_act) > (*cur_all_act)) {
					missing_activity_values.push_back(*cur_all_act);	// current value doesn't appear in feature activity: store it
					cur_all_act++;
				}

				else if ((*cur_feat_act) < (*cur_all_act)) {
					exit(1); //this shouldn't happen
				}
				else {
					cur_feat_act++;
					cur_all_act++;
				}
			}
			*/

			// run K-S test, sensitive version, but with ALL activity values
			(*cur_feat)->determine_significance(act, all_activity_values, feat_activity_values);
		}

		else {
			(*cur_feat)->set_too_infrequent(act); // mark infrequent features
		}
	}

};


