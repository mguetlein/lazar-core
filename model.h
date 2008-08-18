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
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "utils.h"
#include "feature.h"
#include "lazmol.h"

#include <time.h>

#ifndef MODEL_H
#define MODEL_H

using namespace std;
using namespace OpenBabel;

float gauss(float sim, float sigma = 0.3);

template <typename MolType, typename FeatureType, typename ActivityType>
class MetaModel {
	typedef vector<Feature<FeatureType> *> FeatVect;
	typedef FeatMol<MolType,FeatureType,ActivityType> * MolRef;
	typedef vector<FeatMol < MolType, ClassFeat, bool > * > ClassMolVect; 
	typedef vector<FeatMol < MolType, RegrFeat, float > * > RegrMolVect; 
	typedef vector<FeatMol<MolType,FeatureType,ActivityType>*> MolVect; 
        
    public:
        virtual ~MetaModel() {};
        virtual void calculate_prediction(FeatMol < MolType, ClassFeat, bool >* test, ClassMolVect * neighbors, string act){};
        virtual void calculate_prediction(FeatMol < MolType, RegrFeat, float >* test, RegrMolVect * neighbors, string act){};
};

template <typename MolType, typename FeatureType, typename ActivityType>
class Model: public MetaModel<MolType,FeatureType,ActivityType> {
	typedef vector<Feature<FeatureType> *> FeatVect;
	typedef FeatMol<MolType,FeatureType,ActivityType> * MolRef;
	typedef vector<FeatMol < MolType, ClassFeat, bool > * > ClassMolVect; 
	typedef vector<FeatMol < MolType, RegrFeat, float > * > RegrMolVect; 
	typedef vector<FeatMol<MolType,FeatureType,ActivityType>*> MolVect; 

    private:
        vector<string> unknown_features;
        Out* out;

    public:
        Model(Out* out): out(out) {};
        virtual ~Model() {};
        virtual void calculate_prediction(FeatMol<MolType,ClassFeat,bool>* t, ClassMolVect* neighbors, string act);
        virtual void calculate_prediction(FeatMol<MolType,RegrFeat,float>* test, RegrMolVect* neighbors, string act);
};


template <typename MolType, typename FeatureType, typename ActivityType>
class KernelModel: public MetaModel<MolType,FeatureType,ActivityType> {
	typedef vector<Feature<FeatureType> *> FeatVect;
	typedef FeatMol<MolType,FeatureType,ActivityType> * MolRef;
	typedef vector<FeatMol < MolType, ClassFeat, bool > * > ClassMolVect; 
	typedef vector<FeatMol < MolType, RegrFeat, float > * > RegrMolVect; 
	typedef vector<FeatMol<MolType,FeatureType,ActivityType>*> MolVect; 

    private:
        vector<string> unknown_features;
        ActivityType prediction;
        Out* out;

    public:
        KernelModel(Out* out): out(out){};
        virtual ~KernelModel() {};
        virtual void calculate_prediction(FeatMol<MolType,ClassFeat,bool>* t, ClassMolVect * neighbors, string act);
        virtual void calculate_prediction(FeatMol<MolType,RegrFeat,float>* test, RegrMolVect * neighbors, string act);
};



// Implementations
template <typename MolType, typename FeatureType, typename ActivityType>
void Model<MolType,FeatureType,ActivityType>::calculate_prediction(FeatMol<MolType,ClassFeat,bool>* test, ClassMolVect * neighbors, string act) {
    vector<Feature<ClassFeat> *> features = test->get_features();
    this->unknown_features = test->get_unknown();

	float prediction = 0;
	float sim;
	float known_fraction = float(features.size()) / float(features.size() + unknown_features.size()); 

	typename vector<FeatMol<MolType,ClassFeat,bool> *>::iterator cur_n;
	vector<bool> activity;
	vector<bool>::iterator a;

	if (neighbors->size()>1) {

		for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {

			// prediction weighted by fraction of known structure
			sim = (*cur_n)->get_similarity()*known_fraction;
			sim = gauss(sim);
			vector<bool> activity = (*cur_n)->get_act(act);

			for (a = activity.begin(); a != activity.end(); a++) {

				if (*a)
					prediction = prediction + sim;
				else
					prediction = prediction - sim;

			}
		}

		prediction = prediction/neighbors->size();
		*(this->out) << "known_fraction: " << known_fraction << "\n";
		this->out->print();

		if (prediction>0) {
			*out << "prediction: 1\n";
			*out << "confidence: " << prediction << "\n";
			out->print();
		}
		else {
			*out << "prediction: 0\n";
			*out << "confidence: " << prediction << "\n";
			out->print();
		}
	}

	else {
		*(this->out) << "prediction: \n";
		*(this->out) << "confidence: \n";
		this->out->print();
	}
};

template <typename MolType, typename FeatureType, typename ActivityType>
void Model<MolType,FeatureType,ActivityType>::calculate_prediction(FeatMol<MolType,RegrFeat,float>* test, RegrMolVect * neighbors, string act) {
	// data storage
    vector<Feature<RegrFeat>*> lrf;							
	vector<Feature<RegrFeat>*>* lr_features = &lrf;				
	multimap<float, RegrFeat*> msf;					
	multimap<float,FeatMol<MolType,RegrFeat,float>*> ssn;
	multimap<float,FeatMol<MolType,RegrFeat,float>*>* sim_sorted_neighbors = &ssn;

	// iterators
	typename RegrMolVect::iterator cur_n;
	typename FeatVect::iterator cur_feat;
	typename multimap<float, RegrFeat*>::iterator sig_feat;
	
	// linear regression
	gsl_vector* x;
	gsl_matrix* X;
	gsl_matrix* cov;
	gsl_vector* c;
	gsl_multifit_linear_workspace* p_workspace;
	gsl_vector* y;
	gsl_vector* w;

	// primitive data types
	unsigned int no_r, no_c;
	float sim;
	double chisq, y_est, y_err;

	if (neighbors->size()) {

		// sort neighbors by similarity
		sim_sorted_neighbors->clear();
		for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
			sim = (*cur_n)->get_similarity();
			sim_sorted_neighbors->insert(pair<float, FeatMol<MolType,RegrFeat,float>*>(sim,*cur_n));	
		}
		test->extract_neighbors(neighbors, sim_sorted_neighbors);

		// calculate confidence 
		float confidence = test->calculate_confidence(neighbors, act);

		// unite neighbor's with test's features 
		test->unite_features(lr_features, neighbors);

		// determine matrix size
		no_c = (unsigned int) neighbors->size(); 
		
		if (lr_features->size() < no_c) {
			no_c = lr_features->size();
		}

		no_r = neighbors->size();

		if ((no_r >= no_c) && (confidence < 0.995) && (confidence > 0.0)) {

			x = gsl_vector_calloc(1);
			X = gsl_matrix_calloc(no_r,1);

			gsl_vector** x_p = &x;
			gsl_matrix** X_p = &X;
			
			float qdist = 0.0;
			float med_ndist = 0.0;
			float std_ndist = 0.0;
			float max_ndist = 0.0;
			
			//bool tset_interpolates = build_descriptors_pca(lr_features, neighbors, no_c-1, X_p, x_p, act, &qdist, &med_ndist, &std_ndist, &max_ndist);
			test->build_descriptors_pca(lr_features, neighbors, no_c-1, X_p, x_p, act, &qdist, &med_ndist, &std_ndist, &max_ndist);
		
			// apply mahalanobis correction for confidence with weight 0.25
			float norm_med_ndist = 0.0;
			if (max_ndist > 0.0) norm_med_ndist = med_ndist / max_ndist;
			float norm_std_ndist = 0.0;
			if (max_ndist > 0.0) norm_std_ndist = std_ndist / max_ndist;

			// 1. 0.5
			float x = (norm_med_ndist + 0.5 * norm_std_ndist) / 1.5;

			if (x == 0.0) x = 1.0;
			// 1. 0.25
			//confidence = (confidence + 0.25*x) / 1.25;

			if (confidence < 0.0) confidence = 0.0;
			if (confidence > 1.0) confidence = 1.0;
	
			y = gsl_vector_calloc((*X_p)->size1);
			w = gsl_vector_calloc((*X_p)->size1);

			unsigned int rc	= 1;
			for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
				test->set_y_w((*cur_n), y, w, act, rc);
				rc++;
			}

			cov = gsl_matrix_calloc((*X_p)->size2, (*X_p)->size2);
			c = gsl_vector_calloc((*X_p)->size2);
			p_workspace = gsl_multifit_linear_alloc((*X_p)->size1, (*X_p)->size2);

			// do regression
			if ((*X_p)->size1 && (*X_p)->size2) {
				y_est = 0.0; y_err = 0.0; chisq = 0.0;
				
				// learn model and predict activity

				for (unsigned int i=0; i<(*X_p)->size1; i++) {
					for (unsigned int j=0; j<(*X_p)->size2; j++) {
					}
				}
				for (unsigned int j=0; j<(*x_p)->size; j++) {
				}
			
				gsl_multifit_wlinear((*X_p), w, y, c, cov, &chisq, p_workspace);
				gsl_multifit_linear_est((*x_p), c, cov, &y_est, &y_err);
	
				// free learning matrix
				gsl_matrix_free(*X_p);
	
				// output from here
				int df = no_r-1;
				if (df>0) chisq = chisq / df;
				
				*(this->out) << "med_ndist: " << norm_med_ndist << "\n";
                *(this->out) << "std_ndist: " << norm_std_ndist << "\n";
                *(this->out) << "prediction: " << y_est << "\n";
			    *(this->out) << "confidence: " << confidence << "\n";
			}

			else { // end if vector set
                *(this->out) << "prediction: \n";
                *(this->out) << "confidence: \n";
			}

			gsl_vector_free(y);
			gsl_vector_free(*x_p);	
			gsl_multifit_linear_free(p_workspace);
			gsl_vector_free(c);
			gsl_matrix_free(cov);
	
		} // end if at least one significant feature

		else {
            *(this->out) << "prediction: \n";
            *(this->out) << "confidence: \n";
		}
		


	} // end if at least one neighbor
	
	else {
        *(this->out) << "prediction: \n";
        *(this->out) << "confidence: \n";
	}

};

template <typename MolType, typename FeatureType, typename ActivityType>
void KernelModel<MolType,FeatureType,ActivityType>::calculate_prediction(FeatMol<MolType,ClassFeat,bool>* test, ClassMolVect * neighbors, string act) {
    vector<Feature<ClassFeat> *> features = test->get_features();
    this->unknown_features = test->get_unknown();
 
	float confidence = 0.0;
	float sim = 0.0;
   
	float known_fraction = float(features.size()) / float(features.size() + unknown_features.size()); 

	typename vector<FeatMol<MolType,ClassFeat,bool> *>::iterator cur_n;
	vector<bool> activity;
	vector<bool>::iterator a;

	if (neighbors->size()>1) {

		// calculate confidence 
		for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
			sim = (*cur_n)->get_similarity()*known_fraction;
			sim = gauss(sim);
			activity = (*cur_n)->get_act(act);
			for (a = activity.begin(); a != activity.end(); a++) {
				if (*a)	confidence = confidence + sim;
				else confidence = confidence - sim;
			}
		}
		confidence = confidence/neighbors->size();
	
			
		// calculate activities
		gsl_vector* y = gsl_vector_calloc(neighbors->size());
		SEXP yR; PROTECT(yR = allocVector(INTSXP, neighbors->size()));
		unsigned int rc	= 1;

		for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
			activity = (*cur_n)->get_act(act);

			for (a = activity.begin(); a != activity.end(); a++) {
				gsl_vector_set(y, (rc-1), (*a));	// set gsl vector
				INTEGER(yR)[rc-1] = (*a); 
			}

			rc++;
		}
		PROTECT(yR = R_exec("as.factor",yR));
		//R_exec("print", yR);

		gsl_vector* y_bar = gsl_vector_calloc(y->size);
		for (unsigned int i=0; i<y->size; i++) gsl_vector_set(y_bar,i,(gsl_vector_get(y,0)-gsl_vector_get(y,i)));

		if (gsl_vector_isnull(y_bar)) {
            if (gsl_vector_get(y,0) == 0) *(this->out) << "prediction: 0\n";
            else *(this->out) << "prediction: 1\n";
			*(this->out) << "confidence: " << confidence<< "\n";
			*(this->out) << "known_fraction: " << known_fraction <<"\n";
			this->out->print();
            UNPROTECT(2);
		}

		else {
			// calculate gram matrix
			gsl_matrix* gram_matrix = gsl_matrix_calloc(neighbors->size(), neighbors->size());
			test->calculate_gram_matrix(neighbors, gram_matrix, act);

			// convert gram matrix to R kernelMatrix using R util function
			SEXP mr;
			SEXP* gramR = &mr;
			matrix_gsl2R(gramR, gram_matrix);
			PROTECT((*gramR) = R_exec("as.kernelMatrix", (*gramR)));
				
			// learn kernel model
			SEXP e;
			SEXP fun;
			PROTECT(fun = Rf_findFun(Rf_install("ksvm"), R_GlobalEnv));
			if(fun == R_NilValue) {
				fprintf(stderr, "No definition for function.\n");
				UNPROTECT(1);
				exit(1);
			}
			PROTECT(e = allocVector(LANGSXP,6));
			SETCAR(e, fun);
			SETCAR(CDR(e), (*gramR));
			SETCAR(CDR(CDR(e)), yR);
			SETCAR(CDR(CDR(CDR(e))), mkString("matrix"));
			SET_TAG(CDR(CDR(CDR(e))), install("kernel"));
			SETCAR(CDR(CDR(CDR(CDR(e)))), mkString("C-svc"));
			SET_TAG(CDR(CDR(CDR(CDR(e)))), install("type"));
            SETCAR(CDR(CDR(CDR(CDR(CDR(e))))), mkString("1.0"));
            SET_TAG(CDR(CDR(CDR(CDR(CDR(e))))), install("C"));


			SEXP regm;
			PROTECT(regm = R_tryEval(e, R_GlobalEnv, NULL));

			// extract Support Vector indices and create predictive gram matrix
			// get indices of support vectors
			SEXP svR;
			PROTECT(svR = R_exec("SVindex", regm));
				
			// how many sv's
			SEXP lR;
			PROTECT(lR = R_exec("length", svR));
				
			gsl_matrix* pred_matrix = gsl_matrix_calloc(1, INTEGER(lR)[0]);
			test->calculate_pred_matrix(neighbors, pred_matrix, svR);

			// convert predictive gram matrix into correct format
			SEXP  pr;
			SEXP* predR = &pr;
			matrix_gsl2R(predR, pred_matrix);
			SEXP nrR; PROTECT (nrR = allocVector(INTSXP, 1)); INTEGER(nrR)[0] = 1;
			PROTECT((*predR) = R_exec3("matrix",(*predR), nrR));	// enforce row representation
			PROTECT((*predR) = R_exec("as.kernelMatrix", (*predR)));
				
			// predict(regm,pred)
			SEXP pR;
			
			PROTECT(pR = R_exec3("predict", regm, (*predR)));
			PROTECT(pR = R_exec("as.integer", pR));

			UNPROTECT(14);
			
			gsl_matrix_free(pred_matrix);
            if ((INTEGER(pR)[0]-1) == 0) *(this->out) << "prediction: 0\n";
            else *(this->out) << "prediction: 1\n";
			*(this->out) << "confidence: " << confidence << "\n";
			*(this->out) << "known_fraction: " << known_fraction << "\n";
			this->out->print();
		}

		gsl_vector_free(y);
		gsl_vector_free(y_bar);

	}

	else {
		*(this->out) << "prediction: \n";
		*(this->out) << "confidence: \n";
		*(this->out) << "known_fraction: \n";
		this->out->print();
	}

};

template <typename MolType, typename FeatureType, typename ActivityType>
void KernelModel<MolType,FeatureType,ActivityType>::calculate_prediction(FeatMol<MolType,RegrFeat,float>* test, RegrMolVect * neighbors, string act) {

	// data storage
    vector<Feature<RegrFeat>*> lrf;							
	vector<Feature<RegrFeat>*>* lr_features = &lrf;				
	multimap<float, RegrFeat*> msf;
	multimap<float,FeatMol<MolType,RegrFeat,float>*> ssn;
	multimap<float,FeatMol<MolType,RegrFeat,float>*>* sim_sorted_neighbors = &ssn;


	// iterators
	typename RegrMolVect::iterator cur_n, cur_n2;
	typename FeatVect::iterator cur_feat;
	typename multimap<float, RegrFeat*>::iterator sig_feat;
	gsl_vector* y;

	// primitive data types
	float confidence, sim;

	if (neighbors->size()) {

		// sort neighbors by similarity
		sim_sorted_neighbors->clear();
		for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
			sim = (*cur_n)->get_similarity();
            sim_sorted_neighbors->insert(pair<float, FeatMol<MolType,RegrFeat,float>*>(sim,*cur_n));
		}
		test->extract_neighbors(neighbors, sim_sorted_neighbors);

		// calculate confidence 
		confidence = test->calculate_confidence(neighbors, act);

		// unite neighbor's with test's features 
		test->unite_features(lr_features, neighbors);

		if (confidence < 0.995) {

			y = gsl_vector_calloc(neighbors->size());

			unsigned int rc	= 1;
			for (cur_n = neighbors->begin(); cur_n != neighbors->end(); cur_n++) {
				test->set_y((*cur_n), y, act, rc);
				rc++;
			}

			gsl_matrix* gram_matrix = gsl_matrix_calloc(neighbors->size(), neighbors->size());
			test->calculate_gram_matrix(neighbors, gram_matrix, act);
	
			// convert gram matrix to R kernelMatrix using R util function
			SEXP mr;
			SEXP* gramR = &mr;
			matrix_gsl2R(gramR, gram_matrix);
			PROTECT((*gramR) = R_exec("as.kernelMatrix", (*gramR)));
			
			// convert activity values to R vector
			SEXP vr;
			SEXP* aR = &vr;
			vector_gsl2R(aR,y);
			PROTECT((*aR) = R_exec("as.vector", (*aR)));

			// learn kernel model
		        SEXP e;
		        SEXP fun;
			PROTECT(fun = Rf_findFun(Rf_install("ksvm"), R_GlobalEnv));
			if(fun == R_NilValue) {
				fprintf(stderr, "No definition for function.\n");
				UNPROTECT(1);
				exit(1);
			}
                        PROTECT(e = allocVector(LANGSXP,6));
                        SETCAR(e, fun);
                        SETCAR(CDR(e), (*gramR));
                        SETCAR(CDR(CDR(e)), (*aR));
                        SETCAR(CDR(CDR(CDR(e))), mkString("matrix"));
                        SET_TAG(CDR(CDR(CDR(e))), install("kernel"));
                        SETCAR(CDR(CDR(CDR(CDR(e)))), mkString("nu-svr"));
                        SET_TAG(CDR(CDR(CDR(CDR(e)))), install("type"));
                        SETCAR(CDR(CDR(CDR(CDR(CDR(e))))), mkString("0.8"));
                        SET_TAG(CDR(CDR(CDR(CDR(CDR(e))))), install("nu"));



			SEXP regm;
			PROTECT(regm = R_tryEval(e, R_GlobalEnv, NULL));

			// extract Support Vector indices and create predictive gram matrix
			// get indices of support vectors
			SEXP svR;
			PROTECT(svR = R_exec("SVindex", regm));
			
			// how many sv's
			SEXP lR;
			PROTECT(lR = R_exec("length", svR));
			
			// extract predictive gram matrix
			gsl_matrix* pred_matrix = gsl_matrix_calloc(1, INTEGER(lR)[0]);
			test->calculate_pred_matrix(neighbors, pred_matrix, svR);

			// convert predictive gram matrix into correct format
			SEXP  pr;
			SEXP* predR = &pr;
			matrix_gsl2R(predR, pred_matrix);
			SEXP nrR; PROTECT (nrR = allocVector(INTSXP, 1)); INTEGER(nrR)[0] = 1;
			PROTECT((*predR) = R_exec3("matrix",(*predR), nrR));	// enforce row representation
			PROTECT((*predR) = R_exec("as.kernelMatrix", (*predR)));
			
			// predict
			SEXP pR;
			PROTECT(pR = R_exec3("predict", regm, (*predR)));
			prediction = REAL(pR)[0];

			UNPROTECT(11);

			gsl_matrix_free(gram_matrix);
			gsl_matrix_free(pred_matrix);
			gsl_vector_free(y);

			*(this->out) << "prediction: " << prediction << "\n";
			*(this->out) << "confidence: " << confidence << "\n";
			this->out->print();
			
	
		} // end if confidence >...

        else {

            *(this->out) << "prediction: \n";
            *(this->out) << "confidence: \n";
            this->out->print();
    
        }


	} // end if at least one neighbor

    else {
    
        *(this->out) << "prediction: \n";
        *(this->out) << "confidence: \n";
        this->out->print();

    }

};


float gauss(float sim, float sigma) {

	const float invc = 1.0;
	float x = (invc - sim);
	if (x > 1.0) x = 1.0;
	if (x < 0.0) x = 0.0;

	// gauss kernel
	float g = 0.0;
	g = exp(-(x*x)/(2*sigma*sigma));
	return g;
};

#endif
