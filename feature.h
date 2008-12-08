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

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <string>
#include <vector>
#include <iterator>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <map>
#include <sstream>
#include "io.h"

#ifndef FEATURE_H
#define FEATURE_H

using namespace std;
using namespace OpenBabel;


//! the basic feature class
class Feat {

	private:

		string name;

	public:

		Feat() {};
		Feat(string newname): name(newname) {};

		string get_name();
		void set_name(string newname);
		void print(Out * out);
};



//! features for classificatoutn
class ClassFeat: public Feat {

	private:

		map<string, float> na;
		map<string, float> ni;
		map<string, float> fa;
		map<string, float> fi;
		map<string, float> significance;
        map<string, float> p;
		map<string, bool> too_infrequent;

		float cur_sig;
		float cur_p;

		// MG : precompute significance
		static bool cur_str_active;
		bool cur_feat_occurs;
		void precompute_significance(string act);
		string get_map_key(string act);
		// MG

	public:

		ClassFeat(){
			cur_feat_occurs = false;
		};
		ClassFeat(string name): Feat(name) {
			cur_feat_occurs = false;
		};

		//! Determine feature significance using chi-sq test
		void determine_significance(string act, float n_a, float n_i, vector<bool> * activities); // AM: determine significance

		// MG : precompute significance
		void precompute_significance(string act, float n_a, float n_i, float f_a, float f_i);
		static void set_cur_str_active(bool str_active){
			cur_str_active = str_active;
		};
		void set_cur_feat_occurs(bool feat_occurs);
		// MG

		void print_header(Out * out);
		void print(string act,Out * out);
		void print_specifics(string act, Out* out);

		void set_cur_significance(string act) {
			cur_sig = significance[get_map_key(act)]; };
		void set_cur_p(string act) {
			cur_p = p[get_map_key(act)]; };

		float get_significance(string act) {
			return (significance[get_map_key(act)]); };

		float get_cur_significance() { return (cur_sig); };
		float get_cur_p() { return (cur_p); };
		float get_sig_limit();
		float get_p_limit();
		float calc_p(string act);
		float get_p(string act);	//! returns the p value of the feature
		float get_na(string act);
		float get_ni(string act);
		float get_fa(string act);
		float get_fi(string act);

		bool get_too_infrequent(string act) {
			return (too_infrequent[get_map_key(act)]); };
};



//! features for regression
class RegrFeat: public Feat {

	private:

		map<string, float> median;
		map<string, float> global_median;
		map<string, int> feature_larger;
		map<string, int> feature_smaller;
		map<string, int> nr;
		map<string, float> significance;
        map<string, float> p;
		float cur_sig;
		float cur_p;
		map<string, bool> too_infrequent;

	public:

		RegrFeat();
		RegrFeat(string name): Feat(name)  {};

		//! Determine feature significance using KS test
		void determine_significance(string act, float median_all, vector<float> * activities);
		void determine_significance(string act, vector<float> all_activities, vector<float> feat_activities);
		void determine_significance_ks_e(string act, vector<float> all_activities, vector<float> feat_activities);

		//MG:
		void set_cur_feat_occurs(bool feat_occurs){
			fprintf(stderr, "Not implemented for regression");
			exit(1);
		}
		//MG

		void print_header(Out * out);
		void print(string act,Out * out);
		void print_specifics(string act, Out* out);

		float get_significance(string act) { return(significance[act]); };	//! returns the p value of the feature
		void set_cur_significance(string act) { cur_sig = significance[act]; };
		void set_cur_p(string act) { cur_p = p[act]; };
		float get_cur_significance() { return (cur_sig); };
		float get_cur_p() { return (cur_p); };
//		float get_sig_limit();
		float get_p_limit();
		float get_global_median(string act);
		float get_median(string act);
		float calc_p(string act);
		float get_p(string act);	//! returns the p value of the feature

		bool get_too_infrequent(string act) { return (too_infrequent[act]); };
		void set_too_infrequent(string act) { too_infrequent[act] = true; };

};

//! a feature class that is capable to match on OBMol objects
class OBSmartsFrag: public Feat {

	private:

		OBSmartsPattern smarts_pattern;

	public:

		OBSmartsFrag() {};
		OBSmartsFrag(string smarts);
		OBSmartsPattern * get_smarts_pattern() { return(&smarts_pattern); };

};

//! the linear fragment class with refinement methods
class LinFrag {

	private:

		// a vector based representation of SMARTS, needed for refinement
		vector<string> fragment;

	public:

		LinFrag() {};
		LinFrag(string smarts, bool split_bonds);
		LinFrag(vector<string> fragment);

		vector<string> get_fragment();	//!< return the fragment in the vector representation
		int size();

		//! methods for the refinement operator

		string canonify();	//! find a cononical form

		void insert(string element);	//! add new {bonds|atoms} at the beginning of the fragment

		void expand(string element);	//! add new {bonds|atoms} at the end of the fragment

		void rev();	//! reverse the fragment

		string first_atom();
		string first_bond();
		string last_atom();
		string last_bond();

		bool more_specific(LinFrag * g);

		//! methods for rex fragments

		string init_wildcard();
};

//! the linear fragment class with the ability to match LazMol objects
class OBLinFrag: public OBSmartsFrag, public LinFrag {

	private:

		vector<int> pot_matches;

	public:

		OBLinFrag() {};
		OBLinFrag(string sma, bool split_string): OBSmartsFrag(sma), LinFrag(sma, split_string) {};
		OBLinFrag(string sma, LinFrag newfrag, vector<int> pot_matches): OBSmartsFrag(sma), LinFrag(newfrag), pot_matches(pot_matches) { };

};

//! template class for features of various types
template <class FeatureType>
class Feature: public FeatureType {

	private:

		vector<int> matches;
//		map<const int, int> match_freq;

	public:

		Feature() {};
		Feature(string name): FeatureType(name) { };
		Feature(string name, bool split_string): FeatureType(name, split_string) { };
		Feature(string can_sma, LinFrag fragment, vector<int> pot_matches): FeatureType(can_sma,fragment,pot_matches) {};

		void add_match(int comp_nr) { matches.push_back(comp_nr); };

		vector<int> get_matches() { return(matches); };
		vector<int> * get_matches_ptr() { return(&matches); };

		int nr_matches() { return(matches.size()); };

		void print_matches(Out * out, char* smarts);
		void print_all(string act, Out * out, char* smarts);

		bool more_specific(FeatureType * f2);

		void clear_matches() { matches.clear(); }


};

template <class FeatureType>
void Feature<FeatureType>::print_matches(Out * out, char* smarts=NULL) {
	if ((smarts == NULL) || (smarts==this->get_name())) {
		*out << this->get_name() << "\t[ ";
		copy(matches.begin(), matches.end(), ostream_iterator<int>(*out, " "));
		*out << "]\n";
		out->print();
	}
};



template <class FeatureType>
void Feature<FeatureType>::print_all(string act, Out * out, char* smarts=NULL) {
	if ((smarts == NULL) || (smarts==this->get_name())) {
		this->print_specifics(act, out);
		*out << "\t[ ";
		copy(matches.begin(), matches.end(), ostream_iterator<int>(*out, " "));
		*out << "]\n";
		out->print();
	}
};

/*
template <class FeatureType>
void Feature<FeatureType>::print_all(string act, Out * out) {

		float p = this->get_p(act);
		*out << this->get_name() << "\t" << p << "\t" << this->get_significance(act) << "\t" << (int) this->get_fa(act) << "\t" << (int) this->get_fi(act) << "\t" << (int) this->get_na(act) << "\t" << (int)this->get_ni(act) << "\t";
		out->print();

		if (this->get_fa(act)/(this->get_fa(act)+this->get_fi(act)) > this->get_na(act)/(this->get_na(act)+this->get_ni(act)))
				*out << "a";
		else
				*out << "i";
		*out << "\t[ ";
		copy(matches.begin(), matches.end(), ostream_iterator<int>(
*out, " "));
		*out << "]\n";
		out->print();
};
*/

template <class FeatureType>
bool Feature<FeatureType>::more_specific(FeatureType * f2) {

	LinFrag tmp1(this->get_name(),true);
	LinFrag tmp2(f2->get_name(),true);
	return(tmp1.more_specific(&tmp2));

};

template <class FeatureType>
class greater_sig {

typedef Feature<FeatureType> * FeatRef;

public:

    bool operator() (const FeatRef f1,const FeatRef f2) {
            return (f1->get_cur_significance() > f2->get_cur_significance());
    }
};

template <class FeatureType>
class greater_p {

typedef Feature<FeatureType> * FeatRef;

public:

    bool operator() (const FeatRef f1,const FeatRef f2) {
            return (f1->get_cur_p() > f2->get_cur_p());
    }
};

#endif // FEATURE_H
