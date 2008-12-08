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

#include <math.h>
#include "io.h"
#include "stats.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include "openbabel/obconversion.h"
#include "feature.h"
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <list>
#include <time.h>

#include "utils.h"
#include "rutils.h"

#ifndef LAZMOL_H
#define LAZMOL_H

using namespace std;
using namespace OpenBabel;

extern bool quantitative;

//! the basic molecule class
class LazMol {

	private:

		int line_nr;	// line nr -1 of the input file
		string id;	// e.g. name, CAS, NSC, ...
		string smiles;	// SMILES string
		string inchi; // InChI string for the determination of structural identity (use only the basic InChI layer, because stereochemistry is not considered)
		Out * out;

	public:

		LazMol();
		LazMol(int nr, string new_id, string new_smiles, Out * out);

		void set_id(string id);
		void set_inchi(string new_inchi);

		int get_line_nr();
		string get_id();
		string get_inchi();		//!< return InChI string

		void print(); //!< print id and line number
		string get_smiles();	//!< return SMILES string
		void set_output(Out * newout);

};

//! molecule class with an OBMol object for pattern matching
class OBLazMol: public LazMol {

	private:

		OBMol mol;	// OBMol object
		Out * out;

	public:

		OBLazMol(int nr, string id, string new_smiles, Out * out);

		bool match(OBSmartsPattern * smarts_pattern);	//!< match a OBSmartsPattern
		int match_freq(OBSmartsPattern * smarts_pattern);	//!< match a OBSmartsPattern and return the number of matches
		OBMol * get_mol_ref();	//!< return the reference to the corresponding OBMol object
		vector<string> sssr();	//!< identify the smallest set of smallest rings
		void set_output(Out * newout) { out = newout; };

};

template <typename MolType, typename FeatureType, typename ActivityType>
class FeatMol: public MolType {

	typedef vector<Feature<FeatureType> *> FeatVect;
	typedef FeatMol<MolType,FeatureType,ActivityType> * MolRef;
	typedef vector<FeatMol < MolType, ClassFeat, bool > * > ClassMolVect;
	typedef vector<FeatMol < MolType, RegrFeat, float > * > RegrMolVect;
	typedef vector<FeatMol<MolType,FeatureType,ActivityType>*> MolVect;

	private:

		FeatVect features;

		FeatVect pred_features;
		FeatVect common_feat;
		FeatVect sig_pred;
		FeatVect sig_common;

		vector<string> unknown_features;

		map<string, vector<ActivityType> > activities;
		map<string, vector<ActivityType> > db_activities;
		map<string, bool > available;
		map<string, bool > available_bak;

		//! tanimoto distance
		float similarity;

		ActivityType prediction;

		Out * out;

	public:

		FeatMol(int nr): MolType(nr), similarity(0) {};
		FeatMol(int i, string id, string smi): MolType(i, id, smi), similarity(0) {};
		FeatMol(int i, string id, string smi, Out * out): MolType(i, id, smi, out), similarity(0), out(out) {};

		bool find_f_in_n(RegrFeat* f, FeatMol<MolType,RegrFeat,float>* n);

		void extend_matrix(gsl_matrix** X_p, gsl_vector* v);
		void extend_vector(gsl_vector** x_p, float f);

		void run_mahal(gsl_matrix** X_p, gsl_vector** x_p, float* qdist, list<float>* dists); //!< Calculate mahalanobis distance metrics

		void run_pca(gsl_matrix** X_p, gsl_vector** x_p, unsigned int no_c);

		bool equal(gsl_vector* x, gsl_vector* y); //!< Objective feature selection: equality of vectors

		bool singular(gsl_vector* x); //!< Objective feature selection: singularity of feature

		void build_fv(float* qp, gsl_vector* fv, RegrMolVect* n, RegrFeat* f);

		bool build_descriptors_pca(FeatVect* lrf, RegrMolVect* n, unsigned int no_c, gsl_matrix** X_p, gsl_vector** x_p, string act, float* qdist, float* med_ndist, float* std_ndist, float* max_ndist); //!< Build data matrix X using objective feature selection and principal components analysis

		float gauss(float sim, float sigma);  //!< Compute gaussian smoothed similarity

//		void calculate_prediction(ClassMolVect * neighbors, string act);  //!< Calculate prediction based on the set of neighbors

		float calculate_confidence(MolVect* n, string act); //!< Calculate confidence using neighbor similarities and activities

		void unite_features(vector<Feature<FeatureType> *>* lr_features, RegrMolVect * neighbors);

		void extract_neighbors(RegrMolVect* np, multimap<float,MolRef>* sim_sorted_neighbors);

		void set_y(FeatMol<MolType,RegrFeat,float>* cur_n, gsl_vector* y, string act, int rc); //!< Prepare regression: set neighbor activities and weights

		void set_y_w(FeatMol<MolType,RegrFeat,float>* cur_n, gsl_vector* y, gsl_vector* w, string act, int rc); //!< Prepare regression: set neighbor activities and weights

		void calculate_gram_matrix(MolVect * neighbors, gsl_matrix* gram_matrix, string act);

		void calculate_pred_matrix(MolVect * neighbors, gsl_matrix* pred_matrix, SEXP svR);

//		void calculate_prediction(RegrMolVect * neighbors, string act); //!< Calculate prediction using the set of neighbors

		bool is_available(string act) { return(available[act]); };

		bool db_act_available(string act) {
			if (db_activities.size() > 0)
				return(true);
			else
				return(false);
		}
		void print_db_activity(string act, bool loo);

		void copy_activities(MolRef test_mol);

		void clear_db_activities() { db_activities.clear(); };

		map<string, vector<ActivityType> > get_activities() { return activities; };
		map<string, vector<ActivityType> > get_db_activities() { return db_activities; };

		void set_db_activities(map<string, vector<ActivityType> > new_act);

		void replace_db_activities(map<string, vector<ActivityType> > new_act);

		void replace_activities(map<string, vector<ActivityType> > new_act);

		void set_available(map<string, bool> new_avail);

		void clear_act() { activities.clear(); }

		void print_neighbor(string act);

		vector<string> get_unknown() { return(unknown_features); };

		void delete_unknown() { unknown_features.clear(); };

		void print_unknown(string act);

		FeatVect get_features() { return(features); };

		void print_features(string act);

		void add_feature(Feature<FeatureType> * feat);

		void clear_features() { features.clear(); };

		void common_features(MolRef test);

		FeatVect get_common_features() { return(common_feat); };

		void common_features(MolRef m1, MolRef m2, FeatVect* inter, FeatVect* uni);

		void relevant_features(MolRef test, string act);

		//! Determine similarity of two compounds as weighted Tanimoto index
		float get_similarity(MolRef m1, MolRef m2, string act);

		vector<ActivityType> get_act(string act) { return(activities[act]); }

		//! needs inchi strings
		bool equal(const MolRef mol);

		float get_similarity();

		void remove();

		void set_na(string act) { available[act] = false; }

		void set_activity(string name, ActivityType act);

		void restore() { available = available_bak; }

		bool matches(Feature<FeatureType> * feat_ptr);

//		FeatVect get_sig_features();

		void add_unknown(string name) { unknown_features.push_back(name); };

		bool match(Feature<OBLinFrag> * feat_ptr);

		void set_output(Out * newout) { out = newout; }

};

template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::find_f_in_n(RegrFeat* f, FeatMol<MolType,RegrFeat,float>* n) {
	FeatVect fv;
	FeatVect* fvp = &fv;
	(*fvp) = n->get_features();
	typename FeatVect::iterator fv_it;
	for (fv_it = fvp->begin(); fv_it != fvp->end(); fv_it++) {
		if ((*fv_it) == f) return true;
	}
	return false;
}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::extend_matrix(gsl_matrix** X_p, gsl_vector* v) {
	if (v->size != (*X_p)->size1) { cerr << "extend_matrix: vector and matrix have different row size!" << endl; exit(1); }
	gsl_matrix* X_new = gsl_matrix_calloc((*X_p)->size1, ((*X_p)->size2)+1);
	gsl_matrix_view mv = gsl_matrix_submatrix(X_new,0,0,(*X_p)->size1,(*X_p)->size2);
	gsl_matrix_memcpy(&mv.matrix, (*X_p));
	gsl_matrix_set_col(X_new, (X_new->size2)-1, v);
	gsl_matrix_free(*X_p);
	(*X_p) = X_new;
}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::extend_vector(gsl_vector** x_p, float f) {
	gsl_vector* x_new = gsl_vector_calloc((*x_p)->size+1);
	gsl_vector_view vv = gsl_vector_subvector(x_new, 0, (*x_p)->size);
	gsl_vector_memcpy(&vv.vector, (*x_p));
	gsl_vector_set(x_new, (x_new->size)-1, f);
	gsl_vector_free(*x_p);
	(*x_p) = x_new;
}



template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::run_mahal(gsl_matrix** X_p, gsl_vector** x_p, float* qdist, list<float>* dists) {

	// calculate covariance matrix
	gsl_matrix* covm = gsl_matrix_calloc((*X_p)->size2, (*X_p)->size2);
	if ((*X_p)->size2 >1) cov((*X_p),covm);

	// get distances of neighbors
	for (unsigned int i=0; i<(*X_p)->size1; i++) {
		gsl_vector_view vv = gsl_matrix_row((*X_p),i);
		dists->push_back(mahal((&vv.vector), (*X_p), covm));
	}

	// get distance of query structure
	float q_dist = mahal((*x_p), (*X_p), covm);
	(*qdist) = q_dist;
	gsl_matrix_free(covm);

}


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::run_pca(gsl_matrix** X_p, gsl_vector** x_p, unsigned int no_c) {

	// attach x_p as additional row
	gsl_matrix* X_new = gsl_matrix_calloc((((*X_p)->size1)+1), (*X_p)->size2);
	gsl_matrix_view mv = gsl_matrix_submatrix(X_new,0,0,(*X_p)->size1,(*X_p)->size2);
	gsl_matrix_memcpy(&mv.matrix, (*X_p));
	gsl_matrix_set_row(X_new, (X_new->size1)-1, (*x_p));
	gsl_matrix_free(*X_p);
	(*X_p) = X_new;

        // rotation matrix
        gsl_vector* means = gsl_vector_calloc((*X_p)->size2);
        gsl_matrix* rot = pca_cols((*X_p), means, no_c);
        // transform data
        gsl_matrix* X_t_transform = transformData((*X_p), rot, means);
	gsl_matrix* X_transform = gsl_matrix_calloc(X_t_transform->size2, X_t_transform->size1);
       	transposeMatrix(X_t_transform, X_transform);

	gsl_matrix_free(*X_p);
        gsl_vector_free(means);
        gsl_matrix_free(rot);
	gsl_matrix_free(X_t_transform);

	(*X_p) = X_transform;

	// detach x_p from X_p
	X_new = gsl_matrix_calloc((((*X_p)->size1)-1), (*X_p)->size2);
	mv = gsl_matrix_submatrix((*X_p),0,0,X_new->size1,X_new->size2);
	gsl_matrix_memcpy(X_new, &mv.matrix);

	gsl_vector* x_new = gsl_vector_calloc((*X_p)->size2);
	gsl_vector_view vv = gsl_matrix_row((*X_p),(((*X_p)->size1)-1));
	gsl_vector_memcpy(x_new,&vv.vector);

	gsl_matrix_free(*X_p);
	(*X_p) = X_new;
	gsl_vector_free(*x_p);
	(*x_p) = x_new;

}


template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::equal(gsl_vector* x, gsl_vector* y) {
	if (x->size != y->size)	{
		cerr << "Comparing vectors of different size!" << endl; exit(1);
	}
	for (unsigned int i = 0; i<x->size; i++) {
		if (gsl_vector_get(x,i) != gsl_vector_get(y,i)) {
			return false;
		}
	}
	return true;
}

template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::singular(gsl_vector* x) {
	if (x->size <= 3) return false;
	unsigned int cnt=0;
	for (unsigned int i=0; i<x->size; i++) {
		if (gsl_vector_get(x,i) != 0.0) cnt++;
	}
	if (cnt > 1) return false;
	return true;
}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::build_fv(float* qp, gsl_vector* fv, RegrMolVect* n, RegrFeat* f) {
	typename RegrMolVect::iterator n_it;
	unsigned int i = 0;

	(*qp) = 0.0; if (find_f_in_n(f,this)) (*qp) = 1.0;

	for (n_it = n->begin(); n_it != n->end(); n_it++) {
		if (find_f_in_n(f,(*n_it))) gsl_vector_set(fv, i, 1.0);
		else gsl_vector_set(fv, i, 0.0);
		i++;
	}
}


template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::build_descriptors_pca(FeatVect* lrf, RegrMolVect* n, unsigned int no_c, gsl_matrix** X_p, gsl_vector** x_p, string act, float* qdist, float* med_ndist, float* std_ndist, float* max_ndist) {

	#define	FEATURE_POOL_SIZE 100000

	typename FeatVect::iterator cur_feat;
	bool tset_interpolates = true;


	//sort features
	multimap<float, Feature<RegrFeat>*> ms_features;
	typename multimap<float, Feature<RegrFeat>*>::reverse_iterator sf;

	ms_features.clear();
	for (cur_feat = lrf->begin(); cur_feat != lrf->end(); cur_feat++) {
		if (!(*cur_feat)->get_too_infrequent(act)) {
			ms_features.insert(make_pair(((*cur_feat)->get_significance(act)), (*cur_feat)));
		}
	}

	lrf->clear();
	for (sf = ms_features.rbegin(); sf != ms_features.rend(); sf++) {
		lrf->push_back(sf->second);
	}

	unsigned int f_cnt = 0;




	// insert features
	gsl_vector* fv;
	if (no_c) {

		fv = gsl_vector_calloc(n->size());
		if (fv->size != (*X_p)->size1) { cerr << "ofs: feature vector with wrong size!" << endl; exit(1); }

		bool insert;
		float q = 0.0;




		// insert first column
		typename FeatVect::iterator f_it = lrf->begin();
		do {
			if ((f_it == lrf->end()) || (f_cnt == FEATURE_POOL_SIZE)) break;
			build_fv(&q, fv, n, (*f_it));
			insert = true;
			if (gsl_vector_isnull(fv)) {
				insert=false; }
			if (singular(fv)) {
				insert=false; }
			f_it++;
			f_cnt++;

		} while (!insert);

		gsl_matrix_set_col((*X_p), 0, fv);
		gsl_vector_set((*x_p), 0, q);





		// insert the rest, implementing ofs (null, singular, equal)
		if (((*x_p)->size < no_c) && (f_cnt < no_c)) {

			do {
				if ((f_it == lrf->end()) || (f_cnt == FEATURE_POOL_SIZE)) break;

				gsl_vector_set_all(fv,0.0);
				build_fv(&q, fv, n, (*f_it));

				bool extend = true;
				for (unsigned int i=0; i<(*X_p)->size2; i++) {
					gsl_vector_view vv = gsl_matrix_column((*X_p),i);
					if (equal(fv, &vv.vector)) { extend=false; //cerr << "e" << i << "(" << (*X_p)->size2 << ")";
						break; }
				}
				if (gsl_vector_isnull(fv)) { extend = false; //cerr << "n" << "(" << (*X_p)->size2 << ")";
					}
				if (singular(fv)) { extend = false; //cerr << "s" << "(" << (*X_p)->size2 << ")";
					}
				if (extend) { extend_matrix(X_p, fv); extend_vector(x_p, q); //cerr << (*x_p)->size;
					}

				f_it++; //cerr << "-";
				f_cnt++;

			} while ((*x_p)->size < no_c);

		}




		// orthogonalize and de-noise feature space using pca
		unsigned int final_no_c = (unsigned int) (no_c/7);
		if (final_no_c == 0) final_no_c = 1;
		run_pca(X_p, x_p, final_no_c);




		// check if query structure is an outlier using leverage
		if ((*X_p)->size1 > 1) {
			list<float> ds;
			list<float>* ndists = &ds;

			// calculate mahalanobis distances for neigbors and query structure and mean neighbor distance
			run_mahal(X_p, x_p, qdist, ndists);

			list<float>::iterator d_it;
			float dsum, dmedian, dmean, dvar, ddev, dskew, dkurt;
			ndists->sort();
			if (ndists->size() > 3) { ndists->pop_back(); ndists->pop_front(); }
			computeStats(ndists->begin(), ndists->end(), dsum, dmedian, dmean, dvar, ddev, dskew, dkurt);


			// compute maximum of mahal distances
			float dmax = 0.0;
			d_it = ndists->end(); if (d_it != ndists->begin()) { d_it--; dmax = (*d_it); }

			// store return values
			(*med_ndist) = dmedian;
			(*std_ndist) = ddev;
			(*max_ndist) = dmax;


			// leverage normalizer from: C05, p. 124
			float normalizer = 0.0;
			for (d_it = ndists->begin(); d_it != ndists->end(); d_it++) {
				normalizer = normalizer + ((*d_it)*(*d_it));
			}

			// leverage threshold from: C05, p. 124
			float outlier_threshold = 2.0 * ((*X_p)->size2+1) / ((*X_p)->size1);

			// leverage value for query structure : 1/n + qd²/d²
			(*qdist) = ((*qdist)*(*qdist)) / normalizer;
			(*qdist) = (*qdist) + (1.0 / ndists->size());
			if ((*qdist) > outlier_threshold) tset_interpolates = false;

			// neighbor outlier check
			if ((*X_p)->size2 > 1) {
				for (d_it = ndists->begin(); d_it != ndists->end(); d_it++) {
					float current_n_dist = ((*d_it)*(*d_it)) / normalizer;
					current_n_dist = current_n_dist + (1.0 / ndists->size());
					if (current_n_dist > outlier_threshold) cerr << "OUTLIER!!" << endl;
				}
			}

		}




		// insert last column (y-intercept)
		gsl_vector_set_all (fv, 1.0);
		extend_matrix(X_p, fv); extend_vector(x_p, 1.0);
		gsl_vector_free(fv);


	}

	else {

		fv = gsl_vector_calloc(n->size());
		gsl_vector_set_all (fv, 1.0);
		gsl_matrix_set_col((*X_p), 0, fv);
		gsl_vector_set((*x_p), 0, 1.0);
		gsl_vector_free(fv);

	}


	return tset_interpolates;

}


template <typename MolType, typename FeatureType, typename ActivityType>
float FeatMol<MolType,FeatureType,ActivityType>::gauss(float sim, float sigma = 0.3) {

	const float invc = 1.0;
	float x = (invc - sim);
	if (x > 1.0) x = 1.0;
	if (x < 0.0) x = 0.0;

	// gauss kernel
	float g = 0.0;
	g = exp(-(x*x)/(2*sigma*sigma));
	return g;
};


template <typename MolType, typename FeatureType, typename ActivityType>
float FeatMol<MolType,FeatureType,ActivityType>::calculate_confidence(MolVect* n, string act) {

	vector<float> sims;
	vector<float> acts;
	vector<float> activity;
	float confidence = 0.0;
	float sim = 0.0;
	typename RegrMolVect::iterator cur_n;

	float ssum, smedian, smean, svar, sdev, sskew, skurt;
	float asum, amedian, amean, avar, adev, askew, akurt;

	cur_n = n->end();
	if (cur_n != n->begin()) {
		do {
			cur_n--;

			sim = (*cur_n)->get_similarity();
			sim = gauss(sim);
			sims.push_back(sim);

			activity = (*cur_n)->get_act(act);
			sort(activity.begin(), activity.end());
			computeStats(activity.begin(), activity.end(), asum, amedian, amean, avar, adev, askew, akurt);
			acts.push_back(amedian);

		} while (cur_n != n->begin());
	}

	sort(sims.begin(), sims.end());
	computeStats(sims.begin(), sims.end(), ssum, smedian, smean, svar, sdev, sskew, skurt);
	sort(acts.begin(), acts.end());
	computeStats(acts.begin(), acts.end(), asum, amedian, amean, avar, adev, askew, akurt);

	float act_disc = exp(-adev);
	confidence = smedian * act_disc;
	return(confidence);

}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::unite_features(vector<Feature<FeatureType> *>* lr_features, RegrMolVect * neighbors) {

	typename RegrMolVect::iterator cur_n;
	FeatVect tmp_features;
	FeatVect un_features;

	cur_n = neighbors->begin();
	(*lr_features) = (*cur_n)->get_features();
	cur_n++;
	for (; cur_n != neighbors->end(); cur_n++) {
		tmp_features.clear();
		un_features = (*cur_n)->get_features();
		set_union( (lr_features->begin()),(lr_features->end()), (un_features.begin()),(un_features.end()), insert_iterator<FeatVect>(tmp_features,tmp_features.begin()) );
		(*lr_features) = tmp_features;
	}
	un_features = this->get_features();
	set_union( (lr_features->begin()),(lr_features->end()), (un_features.begin()),(un_features.end()), insert_iterator<FeatVect>(tmp_features,tmp_features.begin()) );
	(*lr_features) = tmp_features;
}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::extract_neighbors(RegrMolVect* np, multimap<float,MolRef>* sim_sorted_neighbors) {

	typename multimap<float,MolRef>::iterator cur_sn;

	np->clear();
	cur_sn = sim_sorted_neighbors->end();
	if (cur_sn != sim_sorted_neighbors->begin()) {
		do {
			cur_sn--;
			np->push_back(cur_sn->second);
		} while (cur_sn != sim_sorted_neighbors->begin());
	}
}


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::set_y(FeatMol<MolType,RegrFeat,float>* cur_n, gsl_vector* y, string act, int rc) {

	vector<float> activity;
	float asum, amedian, amean, avar, adev, askew, akurt;

	activity = cur_n->get_act(act);
	amedian = 0;
	sort(activity.begin(), activity.end());
	computeStats(activity.begin(), activity.end(), asum, amedian, amean, avar, adev, askew, akurt);
	gsl_vector_set(y, (rc-1), amedian);

}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::set_y_w(FeatMol<MolType,RegrFeat,float>* cur_n, gsl_vector* y, gsl_vector* w, string act, int rc) {

	float sim = 0.0;
	vector<float> activity;
	float asum, amedian, amean, avar, adev, askew, akurt;

	sim = cur_n->get_similarity();
	sim = gauss(sim);
	gsl_vector_set(w, (rc-1), sim);

	activity = cur_n->get_act(act);
	amedian = 0;
	sort(activity.begin(), activity.end());
	computeStats(activity.begin(), activity.end(), asum, amedian, amean, avar, adev, askew, akurt);
	gsl_vector_set(y, (rc-1), amedian);

}


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::calculate_gram_matrix(MolVect * neighbors, gsl_matrix* gram_matrix, string act) {
	typename MolVect::iterator cur_n, cur_n2;

	// calculate upper right
	unsigned int n1=0;
	for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
		unsigned int n2=n1;
		for (cur_n2 = cur_n; cur_n2 != neighbors->end(); cur_n2++) {
			float tan = gauss(get_similarity((*cur_n), (*cur_n2), act));
			gsl_matrix_set(gram_matrix,n1,n2,tan);
			n2++;
		}
		n1++;
	}

	// add transpose
	gsl_matrix* ll = gsl_matrix_calloc(neighbors->size(), neighbors->size());
	gsl_matrix_transpose_memcpy (ll,gram_matrix);
	gsl_matrix_add(gram_matrix, ll);

	// diagonal has been set twice
	for (unsigned int i = 0; i<neighbors->size(); i++) {
		gsl_matrix_set(gram_matrix,i,i,(gsl_matrix_get(gram_matrix,i,i)/2.0));
	}

	gsl_matrix_free(ll);

}




template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::calculate_pred_matrix(MolVect * neighbors, gsl_matrix* pred_matrix, SEXP svR) {
	typename MolVect::iterator cur_n;
	unsigned int i=0;
	unsigned int j=0;

	SEXP lR;
	PROTECT(lR = R_exec("length", svR));

	for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
		i++;
		if (((unsigned) INTEGER(svR)[j] == i) && (j < ((unsigned) INTEGER(lR)[0]))) {
			j++;
			gsl_matrix_set(pred_matrix,0,j-1,gauss((*cur_n)->get_similarity()));
		}
	}
}


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::print_db_activity(string act, bool loo ) {


	if (loo) {
		if (available_bak[act]) {
			vector<ActivityType> vals = activities[act];
			typename vector<ActivityType>::iterator cur_v;
            *out << "db_activity: [";
			for (cur_v = vals.begin(); cur_v != vals.end(); cur_v++) {
                if (cur_v != vals.begin()) *out << ", ";
				*out << *cur_v;
    		}
			*out << "]\n";
		}
	}


	else {
		if (db_activities.size() > 0) {
			vector<ActivityType> vals = db_activities[act];
			typename vector<ActivityType>::iterator cur_v;
            *out << "db_activity: [";
			for (cur_v = vals.begin(); cur_v != vals.end(); cur_v++) {
                if (cur_v != vals.begin()) *out << ", ";
			    *out << *cur_v;
			}
			*out << "]\n";
		}
	}

	out->print();

};


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::copy_activities(MolRef test_mol) {
	test_mol->set_db_activities(activities);
};

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::set_db_activities(map<string, vector<ActivityType> > new_act) {

	typename map<string, vector<ActivityType> >::iterator cur_map;
	typename vector<ActivityType>::iterator cur_act;
	for (cur_map=new_act.begin();cur_map!=new_act.end();cur_map++) {
		for (cur_act=((*cur_map).second).begin();cur_act!=((*cur_map).second).end();cur_act++) {
			db_activities[(*cur_map).first].push_back(*cur_act);
		}
	}
};

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::replace_db_activities(map<string, vector<ActivityType> > new_act) {
	typename map<string, vector<ActivityType> >::iterator cur_map;
	typename vector<ActivityType>::iterator cur_act;
	db_activities.clear();
	for (cur_map=new_act.begin();cur_map!=new_act.end();cur_map++) {
		for (cur_act=((*cur_map).second).begin();cur_act!=((*cur_map).second).end();cur_act++) {
			db_activities[(*cur_map).first].push_back(*cur_act);
		}
	}
};


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::replace_activities(map<string, vector<ActivityType> > new_act) {
	typename map<string, vector<ActivityType> >::iterator cur_map;
	typename vector<ActivityType>::iterator cur_act;
	activities.clear();
	for (cur_map=new_act.begin();cur_map!=new_act.end();cur_map++) {
		for (cur_act=((*cur_map).second).begin();cur_act!=((*cur_map).second).end();cur_act++) {
			activities[(*cur_map).first].push_back(*cur_act);
		}
	}
};


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::set_available(map<string, bool> new_avail) {
	available = new_avail;
	available_bak = new_avail;
};

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::print_neighbor(string act) {

	typename vector<ActivityType>::iterator cur_act;

	*out << "    - line_nr: " << this->get_line_nr()+1 << "\n";
	*out << "      id: " << this->get_id() << "\n";
	*out << "      smiles: '" << this->get_smiles() << "'\n";
	*out << "      inchi: '" << this->get_inchi() << "'\n";
    *out << "      activity: [";
	for (cur_act = activities[act].begin(); cur_act != activities[act].end(); cur_act++) {
        if (cur_act != activities[act].begin()) *out << ", ";
		*out << *cur_act;
	}
    *out << "]\n";
	*out << "      similarity: " << this->get_similarity() << "\n";
	out->print();
};

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::print_unknown(string act) {

	// print infrequent features
	/*
	if (infrequent_features.size()>0) {
		vector<string>::iterator cur_str;
		for (cur_str = infrequent_features.begin(); cur_str != infrequent_features.end(); cur_str++) {
			unknown_features.push_back(*cur_str);
		}
	}
	*/

	if (unknown_features.size()>0) {

		// determine the most general unknown linear fragments //

		vector<string>::iterator cur_str1;
		vector<string>::iterator cur_str2;
		for (cur_str1 = unknown_features.begin(); cur_str1 != unknown_features.end(); cur_str1++) {
		bool general = true;
			for (cur_str2 = unknown_features.begin(); cur_str2 != unknown_features.end(); cur_str2++) {
				if (cur_str1->size() > cur_str2->size()) {
					string rev_str2 = *cur_str2;
					reverse(rev_str2.begin(),rev_str2.end());
					if ((cur_str1->find(*cur_str2) != string::npos) || (cur_str1->find(rev_str2) != string::npos)) {
						general = false;
						break;
					}
				}
			}
			if (general) {
				*out << "    - '" << *cur_str1 << "'\n";
				out->print();
			}
		}
	}

};

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::print_features(string act) {

	if (features.size()>0) {

		bool redundant;
		FeatVect nonred;
		vector<int> * f1_matches;
		vector<int> * f2_matches;
		vector<int> common_matches;
		typename FeatVect::iterator cur_feat;
		typename FeatVect::iterator cur_nonr;

		sort(features.begin(),features.end(),greater_p<FeatureType>());

		// determine nonredundant features
		for (cur_feat=features.begin();cur_feat!=features.end();cur_feat++) {

			f1_matches = (*cur_feat)->get_matches_ptr() ;
			redundant = false;
			sort(f1_matches->begin(),f1_matches->end());

			for (cur_nonr = nonred.begin(); cur_nonr != nonred.end(); cur_nonr++) {

				f2_matches = (*cur_nonr)->get_matches_ptr();
				sort(f2_matches->begin(),f2_matches->end());
				common_matches.clear();
				set_union(f1_matches->begin(),f1_matches->end(),
						f2_matches->begin(), f2_matches->end(),
						insert_iterator<vector<int> >(common_matches,common_matches.begin()));

				if ((common_matches == *f1_matches) | (common_matches == *f2_matches)) {
					redundant = true;
					break;
				}
			}

			if (!redundant) {
				nonred.push_back(*cur_feat);
				if ((*cur_feat)->get_p(act) >= (*cur_feat)->get_p_limit()) {
					(*cur_feat)->print(act,out);
				}
			}
		}
	}
}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::add_feature(Feature<FeatureType> * feat) {
	features.push_back(feat);
}


template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::common_features(MolRef m1, MolRef m2, FeatVect* inter, FeatVect* uni) {

	inter->clear();
	uni->clear();

	FeatVect m1_features = m1->get_features();
	FeatVect m2_features = m2->get_features();

	sort(m1_features.begin(), m1_features.end());
	sort(m2_features.begin(), m2_features.end());

	set_union(m1_features.begin(),m1_features.end(),
			m2_features.begin(),m2_features.end(),
			insert_iterator<FeatVect>((*uni),(*uni).begin()));

	set_intersection(m1_features.begin(),m1_features.end(),
			m2_features.begin(),m2_features.end(),
			insert_iterator<FeatVect>((*inter),(*inter).begin()));

};

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::common_features(MolRef test) {

	pred_features.clear();
	common_feat.clear();

	FeatVect test_features = test->get_features();
	FeatVect train_features = this->get_features();

	sort(test_features.begin(),test_features.end());
	sort(train_features.begin(),train_features.end());

	set_union(train_features.begin(),train_features.end(),
			test_features.begin(),test_features.end(),
			insert_iterator<FeatVect >(pred_features,pred_features.begin()));

	set_intersection(train_features.begin(),train_features.end(),
			test_features.begin(),test_features.end(),
			insert_iterator<FeatVect >(common_feat,common_feat.begin()));

};

template <typename MolType, typename FeatureType, typename ActivityType>
float FeatMol<MolType,FeatureType,ActivityType>::get_similarity(MolRef m1, MolRef m2, string act) {
	float tanimoto;
	typename FeatVect::iterator cur_feat;

	float c=0;
	float u=0;
	float p=0;

	typedef vector<Feature<FeatureType> *> FeatVect;
	FeatVect suni, sinter, uv, iv;

	// compose union and intersect sets
	if (m1 != this) {
		// sim between two training compounds: calculate sets
		common_features(m1,m2,&iv,&uv);
	}
	else {
		// sim between one training compound and test compound: re-use sets
		uv = this->pred_features;
		iv = this->common_feat;
	}

	// set significances in union feature set
	for (cur_feat = uv.begin(); cur_feat != uv.end(); cur_feat++) {
		(*cur_feat)->set_cur_p(act);
	}

	// determine u = sum of significances of union feature set
	for (cur_feat = uv.begin(); cur_feat != uv.end(); cur_feat++) {
		p = (*cur_feat)->get_cur_p();
		if (p >= (*cur_feat)->get_p_limit()) { // (*)
			p = gauss(p);
			suni.push_back(*cur_feat);
			u = u + p;
		}
	}

	tanimoto = 0.0;
	if (suni.size() > 1) {
		// intersect significance endowed features with common features
		sort(suni.begin(),suni.end(),greater_p<FeatureType>());
		sort(iv.begin(),iv.end(),greater_p<FeatureType>());
		set_intersection(iv.begin(),iv.end(),
				suni.begin(),suni.end(),
				insert_iterator<FeatVect >(sinter,sinter.begin()),
				greater_p<FeatureType>());

		// determine c = sum of significances of common feature set
		for (cur_feat = sinter.begin(); cur_feat != sinter.end(); cur_feat++) {
		    p = (*cur_feat)->get_cur_p();
			if (p >= (*cur_feat)->get_p_limit()) { // test is needed if (*) is not applied to suni.push_back()
				p = gauss(p);
				c = c + p;
			}
		}

		// determine Tanimoto index and store as similarity
		if (u>0) {
			tanimoto = c/u;
		}
	}

	// for similarity to test structure, additionally punish unknown fragments
	if (m1 == this) {
		float known_fraction = float(features.size()) / float(features.size() + unknown_features.size());
		tanimoto = known_fraction * tanimoto;
	}

	return tanimoto;

};



template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::relevant_features(MolRef test, string act) {

	similarity = get_similarity(this, test, act);

};

template <typename MolType, typename FeatureType, typename ActivityType>
float FeatMol<MolType,FeatureType,ActivityType>::get_similarity() {

	if (similarity > 0)
		return(similarity);
	else
		return(0);

};

template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::equal(const MolRef mol) {

	if (this->get_inchi() != "" && mol->get_inchi() != "")
			return(this->get_inchi() == mol->get_inchi());
	else
		return(this->get_smiles() == mol->get_smiles()); // just an attempt to guess identity in the case, that unique smiles have been provided

}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::remove() {

	map<string, bool>::iterator cur_a;

	available_bak = available;

	for (cur_a = available.begin(); cur_a != available.end(); cur_a++) {
		cur_a->second = false;
	}
}

template <typename MolType, typename FeatureType, typename ActivityType>
void FeatMol<MolType,FeatureType,ActivityType>::set_activity(string name, ActivityType act) {
	activities[name].push_back(act);
	available[name] = true;
};

template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::matches(Feature<FeatureType> * feat_ptr) {

	bool match = false;
	typename FeatVect::iterator cur_feat;

	for (cur_feat = features.begin(); cur_feat != features.end(); cur_feat++) {

		if ((*cur_feat) == feat_ptr) {
			match = true;
			break;
		}

	}

	return(match);

};

/*
template <typename MolType, typename FeatureType, typename ActivityType>
vector<Feature<FeatureType> *> FeatMol<MolType,FeatureType,ActivityType>::get_sig_features() {

	typedef vector<Feature<FeatureType> *> FeatVect;
	FeatVect sig_features;
	typename FeatVect::iterator cur_feat;

	for (cur_feat=features.begin();cur_feat!=features.end();cur_feat++) {

		if ((*cur_feat)->get_significance() >= (*cur_feat)->get_sig_limit()) {
			sig_features.push_back(*cur_feat);
		}

	}

	return(sig_features);

};
*/

template <typename MolType, typename FeatureType, typename ActivityType>
bool FeatMol<MolType,FeatureType,ActivityType>::match(Feature<OBLinFrag> * feat_ptr) {
	OBSmartsPattern * sp = feat_ptr->get_smarts_pattern();
	OBMol * mol = this->get_mol_ref();
	return ( sp->Match(*mol, true) );
};

template <typename MolType, typename FeatureType, typename ActivityType>
class greater_sim {

	typedef FeatMol<MolType,FeatureType,ActivityType> * MolRef;

	public:

        bool operator() (const MolRef m1,const MolRef m2) {
			return (m1->get_similarity() > m2->get_similarity());
        }
};

#endif // LAZMOL_H
