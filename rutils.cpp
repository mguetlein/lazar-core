#include "rutils.h"

void init_R(int argc, char **argv) {
	int defaultArgc = 1;
	char *defaultArgv[] = {(char*)"Rtest"};

	if(argc == 0 || argv == NULL) {
		argc = defaultArgc;
		argv = defaultArgv;
	}
	Rf_initEmbeddedR(argc, argv);
}


void end_R() {
	Rf_endEmbeddedR(0);
}

SEXP R_exec3 (const char* command, SEXP structure1, SEXP structure2) {
	SEXP e;
	SEXP val = NILSXP;
	int errorOccurred;

	PROTECT(e = lang3(install((char*) command), structure1, structure2));
	val = R_tryEval(e, R_GlobalEnv, &errorOccurred);
	UNPROTECT(1);
	if(!errorOccurred) {
		return(val);
	}
	else {
		return(NILSXP);
	}
}

SEXP R_exec4 (const char* command, SEXP structure1, SEXP structure2, SEXP structure3) {
	SEXP e;
	SEXP val = NILSXP;
	int errorOccurred;

	PROTECT(e = lang4(install((char*) command), structure1, structure2, structure3));
	val = R_tryEval(e, R_GlobalEnv, &errorOccurred);
	UNPROTECT(1);
	if(!errorOccurred) {
		return(val);
	}
	else {
		return(NILSXP);
	}
}

SEXP R_exec (const char* command, SEXP structure) {
	SEXP e;
	SEXP val = NILSXP;
	int errorOccurred;

	PROTECT(e = lang2(install((char*) command), structure));
	val = R_tryEval(e, R_GlobalEnv, &errorOccurred);
	UNPROTECT(1);
	if(!errorOccurred) {
		return(val);
	}
	else {
		return(NILSXP);
	}
}


SEXP get_list_element(SEXP list, char *str) {
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}


void print_matrix(gsl_matrix* m) {
	for(unsigned int i = 0; i < m->size1; i++) {
		for(unsigned int j = 0; j < m->size2; j++) {
			cerr << fixed << showpos << setprecision(3) << gsl_matrix_get(m,i,j) << "\t";
	        }
		cerr << endl;
	}
}

void print_vector(gsl_vector* v) {
	for(unsigned int i = 0; i < v->size; i++) {
		cerr << fixed << showpos << setprecision(3) << gsl_vector_get(v,i) << "\t";
	}
	cerr << endl;
}

void matrix_gsl2R(SEXP* mr, gsl_matrix* m) {

        PROTECT((*mr) = allocMatrix(REALSXP, m->size1, m->size2));
        for(unsigned int i = 0; i < m->size1; i++) {
                for(unsigned int j = 0; j < m->size2; j++) {
                        REAL((*mr))[i+(m->size1)*j] = gsl_matrix_get(m, i, j);
                }
        }
}

void vector_gsl2R(SEXP* vr, gsl_vector* v) {
        PROTECT((*vr) = allocVector(REALSXP, v->size));
        for(unsigned int i = 0; i < v->size; i++) {
              REAL((*vr))[i] = gsl_vector_get(v, i);
        }
}

void vector_gsl2R_str(SEXP* vr, gsl_vector* v) {
        PROTECT((*vr) = allocVector(STRSXP, v->size));
        for(unsigned int i = 0; i < v->size; i++) {
		char s[1];
		sprintf(s,"%d",(int)gsl_vector_get(v,i));
        	SET_STRING_ELT((*vr),i, mkChar(s));
        }
}

void matrix_R2gsl(gsl_matrix* m, SEXP* mr) {
	for(unsigned int i = 0; i < m->size1; i++) {
                for(unsigned int j = 0; j < m->size2; j++) {
                        gsl_matrix_set(m, i, j, REAL((*mr))[i+(m->size1)*j]);
                }
        }
}

void vector_R2gsl(gsl_vector* v, SEXP* vr) {
	for(unsigned int i = 0; i < v->size; i++) {
        	gsl_vector_set(v, i, REAL((*vr))[i]);
        }
}
float getVectorMean(gsl_vector* v) {
	float m = 0.0;
	for (unsigned int i =0; i < v->size; i++) {
		m += gsl_vector_get(v, i);
	}
	m /= v->size;
	return(m);
}

void transposeMatrix(gsl_matrix* m, gsl_matrix* mt) {
	for (unsigned int i = 0; i < m->size1; i++ ) {
		for (unsigned int j = 0; j < m->size2; j++ ) {
			gsl_matrix_set(mt, j, i, gsl_matrix_get(m,i,j));
		}
	}
}



void cov(gsl_matrix* m, gsl_matrix* covm) {

	// initialise R matrix
	SEXP mr;
	PROTECT(mr = allocMatrix(REALSXP, m->size1, m->size2));
	for(unsigned int i = 0; i < m->size1; i++) {
                for(unsigned int j = 0; j < m->size2; j++) {
                        REAL(mr)[i+(m->size1)*j] = gsl_matrix_get(m, i, j);
                }
        }

	// calculate covariance
	SEXP covmr;
	PROTECT(covmr = R_exec("cov", mr));

	// claim covariance matrix
	//printf("Covariance matrix: \n");
	for(unsigned int i = 0; i < covm->size1; i++) {
                for(unsigned int j = 0; j < covm->size2; j++) {
                        gsl_matrix_set(covm, i, j, REAL(covmr)[i+(covm->size1)*j]);
			//printf("%f ", REAL(covmr)[i+(covm->size1)*j]);
                }
		//printf("\n");
        }
	//fflush(stdout);

	UNPROTECT(2);

}

void get_mean(gsl_matrix* m, gsl_vector* mean) {

	// initialise R matrix
	SEXP mr;
	PROTECT(mr = allocMatrix(REALSXP, m->size1, m->size2));
	for(unsigned int i = 0; i < m->size1; i++) {
                for(unsigned int j = 0; j < m->size2; j++) {
                        REAL(mr)[i+(m->size1)*j] = gsl_matrix_get(m, i, j);
                }
        }

	// calculate mean
	SEXP meanr;
	//printf("Centroid: \n");
	PROTECT(meanr = R_exec("colMeans",mr)); for(unsigned int i = 0; i < mean->size; i++) { gsl_vector_set(mean, i, REAL(meanr)[i]); //printf("%f ", REAL(meanr)[i]);
						} ;
	//printf("\n");
	fflush(stdout);
	UNPROTECT(2);
}


/* Computes mahalanobis distance  between centroid of x and m */

float mahal(gsl_vector* x, gsl_matrix* m, gsl_matrix* covm) {

	gsl_vector* mean = gsl_vector_calloc(m->size2);

	// calculate mean vector
	get_mean(m,mean);

	float d = 0.0;

	// cov can be used
	if (m->size2 > 1) {

		SEXP xr, meanr, covmr;

		//gsl_matrix* covm = gsl_matrix_calloc(m->size2, m->size2);
		// calculate covariance
		//cov(m,covm);

		// set x vector
		PROTECT(xr = allocVector(REALSXP, x->size)); for (unsigned int i=0; i<x->size; i++) REAL(xr)[i] = gsl_vector_get(x,i);
		// set mean vector
		PROTECT(meanr = allocVector(REALSXP, mean->size)); for (unsigned int i=0; i<x->size; i++) REAL(meanr)[i] = gsl_vector_get(mean,i);
		// set cov matrix
		PROTECT(covmr = allocMatrix(REALSXP, covm->size1, covm->size2));
		for(unsigned int i = 0; i < covm->size1; i++) {
			for(unsigned int j = 0; j < covm->size2; j++) {
				REAL(covmr)[i+(covm->size1)*j] = gsl_matrix_get(covm, i, j);
			}
		}


		// R CALC
		SEXP dist;
		PROTECT(dist = R_exec4("mahalanobis", xr, meanr, covmr));
		d = sqrt(REAL(dist)[0]);
		UNPROTECT(4);
		//gsl_matrix_free(covm);
	}

	// cov can not be used -- use euclidean distance for 1-D
	else {
		d = fabs(gsl_vector_get(x,0) - gsl_vector_get(mean,0));
	}

	gsl_vector_free(mean);
	return d;
}


/*
int main(int argc, char *argv[]) {
	char *localArgs[] = {"R", "--silent"};
	init_R(sizeof(localArgs)/sizeof(localArgs[0]), localArgs);

	gsl_vector* x = gsl_vector_calloc(2);
	gsl_vector* center = gsl_vector_calloc(2);
	gsl_matrix* cov = gsl_matrix_calloc(2,2);

	gsl_vector_set(x,0,2);
	gsl_vector_set(x,1,1);

	gsl_vector_set(center,0,1);
	gsl_vector_set(center,1,1);

	gsl_matrix_set(cov,0,0,0.2);
	gsl_matrix_set(cov,1,0,0.2);
	gsl_matrix_set(cov,0,1,0.2);
	gsl_matrix_set(cov,1,1,4.0);

	printf("%f\n\n", mahal(x, center, cov));

	gsl_vector_free(x);
	gsl_vector_free(center);
	gsl_matrix_free(cov);

	end_R();
	return(0);
}
*/





gsl_matrix* pca_cols(gsl_matrix* feature_matrix, gsl_vector* means, unsigned int no_c) {

	// subtract means of columns
	for (unsigned int j = 0; j < feature_matrix->size2; j++) {
		gsl_vector_view vv = gsl_matrix_column(feature_matrix,j);
		gsl_vector* v = &vv.vector;
		gsl_vector_set(means, j, getVectorMean(v));
		gsl_vector_add_constant(v, (-1.0) * gsl_vector_get(means,j));
	}

	// initialise matrix
	//fprintf(stderr, "feature_matrix: %ix%i\n", feature_matrix->size1, feature_matrix->size2);
	SEXP m;
	double* matrix;
	PROTECT(m = allocMatrix(REALSXP, feature_matrix->size1, feature_matrix->size2));
	matrix = REAL(m);
	for(unsigned int i = 0; i < feature_matrix->size1; i++) {
		for(unsigned int j = 0; j < feature_matrix->size2; j++) {
			matrix[i+(feature_matrix->size1)*j] = gsl_matrix_get( feature_matrix, i, j);
	        }
	}

	// do principal components analysis, using R
	//fprintf(stderr, "PCA\n"); fflush(stdout);
	SEXP pca;
	PROTECT(pca = R_exec("prcomp", m));							//R_exec("print", pca);
	SEXP summary;
	PROTECT(summary = R_exec("summary", pca));						//R_exec("print", summary);


	// get proportion of variance
	SEXP ev;
	PROTECT(ev = get_list_element(pca,(char*)"sdev"));						//R_exec("print",ev);
	unsigned int dim = length(ev);								//fprintf(stderr, "dim: %i\n", dim);

	float sum_var = 0.0; float c_ev = 0.0;
	for (unsigned int i = 0; i < dim; i++) {
		c_ev = (REAL(ev)[i]) * (REAL(ev)[i]);
		sum_var += c_ev;
	}

	// using all eigenvectors?
	bool using_all = false;
	if (dim <= no_c) { no_c = dim; using_all = true; }

	float cum_var = 0.0; c_ev = 0.0;
	for (unsigned int i = 0; i < no_c; i++) {
		c_ev = (REAL(ev)[i]) * (REAL(ev)[i]);
		cum_var += c_ev;
	}
	//if (using_all) fprintf(stderr, "Cumulative variance of %g reached by using ALL %i eigen vector(s).\n" , (cum_var/sum_var), no_c);
	//else fprintf(stderr, "Cumulative variance of %g reached by using %i eigen vector(s).\n" , (cum_var/sum_var), no_c);

	// get loads (eigenvectors)
	SEXP loads;
	PROTECT(loads = get_list_element(pca, (char*)"rotation"));					//R_exec("print", loads);
	gsl_matrix* rot = gsl_matrix_alloc(dim, no_c);
	for(unsigned int i = 0; i < dim; i++) {
		for(unsigned int j = 0; j < no_c; j++) {
			gsl_matrix_set(rot, i, j, REAL(loads)[i+dim*j]);			//printf("%g \n", REAL(loads)[i+dim*j]);
	        }
	}

	// de-initialise R
	UNPROTECT(4);
	end_R();
	return(rot);
}


gsl_matrix* pca(gsl_matrix* feature_matrix, gsl_vector* means, float sig_limit) {

	// subtract means of columns
	for (unsigned int j = 0; j < feature_matrix->size2; j++) {
		gsl_vector_view vv = gsl_matrix_column(feature_matrix,j);
		gsl_vector* v = &vv.vector;
		gsl_vector_set(means, j, getVectorMean(v));
		gsl_vector_add_constant(v, (-1.0) * gsl_vector_get(means,j));
	}

	// initialise matrix
	SEXP m;
	double* matrix;
	PROTECT(m = allocMatrix(REALSXP, feature_matrix->size1, feature_matrix->size2));
	matrix = REAL(m);
	for(unsigned int i = 0; i < feature_matrix->size1; i++) {
		for(unsigned int j = 0; j < feature_matrix->size2; j++) {
			matrix[i+(feature_matrix->size1)*j] = gsl_matrix_get( feature_matrix, i, j);
	        }
	}

	// do principal components analysis, using R
	//fprintf(stderr, "PCA\n"); fflush(stdout);
	SEXP pca;
	PROTECT(pca = R_exec("prcomp", m));							//R_exec("print", pca);
	SEXP summary;
	PROTECT(summary = R_exec("summary", pca));						//R_exec("print", summary);


	// get proportion of variance
	SEXP ev;
	PROTECT(ev = get_list_element(pca,(char*)"sdev"));						//R_exec("print",ev);
	unsigned int dim = length(ev);									//printf("dim: %i\n", dim);

	float sum_var = 0.0; float c_ev = 0.0;
	for (unsigned int i = 0; i < dim; i++) {
		c_ev = (REAL(ev)[i]) * (REAL(ev)[i]);
		sum_var += c_ev;
	}

	float cum_var = 0.0; unsigned int sig_cnt = 0;
	for (unsigned int i = 0; i < dim; i++) {
		c_ev = (REAL(ev)[i]) * (REAL(ev)[i]);
		cum_var += c_ev;							 	//printf("ev%i: %.7g\n", i, REAL(ev)[i]);
		sig_cnt++;
		if ((cum_var/sum_var) > sig_limit) break;
	}

	//fprintf(stderr, "Cumulative variance of %g reached by using %i eigen vector(s).\n" , (cum_var/sum_var), sig_cnt);

	// get loads (eigenvectors)
	SEXP loads;
	PROTECT(loads = get_list_element(pca, (char*)"rotation"));					//R_exec("print", loads);
	gsl_matrix* rot = gsl_matrix_alloc(dim, sig_cnt);
	for(unsigned int i = 0; i < dim; i++) {
		for(unsigned int j = 0; j < sig_cnt; j++) {
			gsl_matrix_set(rot, i, j, REAL(loads)[i+dim*j]);			//printf("%g \n", REAL(loads)[i+dim*j]);
	        }
	}

	// de-initialise R
	UNPROTECT(4);
	end_R();
	return(rot);
}

gsl_matrix* transformData (gsl_matrix* data_c, gsl_matrix* rot, gsl_vector* means) {
	// get transpose of rotation matrix
	gsl_matrix* rot_t = gsl_matrix_alloc(rot->size2, rot->size1);
	//fprintf(stderr, "rot: %ix%i\n", rot->size1, rot->size2);
	transposeMatrix(rot, rot_t);

	// get transpose of zero centered data
	gsl_matrix* data_ct = gsl_matrix_alloc(data_c->size2, data_c->size1);
	transposeMatrix(data_c, data_ct);

	// t_data = rot_t * data_ct
	gsl_matrix* t_data = gsl_matrix_alloc(rot_t->size1, data_ct->size2);
	//fprintf(stderr, "rot_t: %ix%i, data_ct: %ix%i\n", rot_t->size1, rot_t->size2, data_ct->size1, data_ct->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rot_t, data_ct, 0.0, t_data);

	gsl_matrix_free(rot_t);
	gsl_matrix_free(data_ct);
	return(t_data);
}

gsl_matrix* reconstructData (gsl_matrix* t_data, gsl_matrix* rot, gsl_vector* means) {

	// rec_data = rot * t_data
	gsl_matrix* rec_data_t = gsl_matrix_alloc(rot->size1, t_data->size2);
	gsl_blas_dgemm( CblasNoTrans,
			CblasNoTrans,
			1.0, rot, t_data, 0.0, rec_data_t);


	gsl_matrix* rec_data = gsl_matrix_alloc(rec_data_t->size2, rec_data_t->size1);
	transposeMatrix(rec_data_t, rec_data);

	gsl_matrix_free(rec_data_t);

	// rec_data + means
	for (unsigned int j = 0; j < (rec_data->size2); j++) {
		gsl_vector_view vv = gsl_matrix_column(rec_data,j);
		gsl_vector* v = &vv.vector;
		gsl_vector_add_constant(v, gsl_vector_get(means,j));
	}

	return(rec_data);
}


/*
int
main() {

	// initialise R
	printf("INIT R, "); fflush(stdout);
	char* R_argv[] = {"REmbeddedPostgres", "--gui=none", "--silent"};
	int R_argc = sizeof(R_argv)/sizeof(R_argv[0]);
	putenv("R_HOME=/usr/local/lib/R");
       	init_R(R_argc, R_argv);


	// prepare data
	gsl_matrix* data;
	int nx = 10;
	int ny = 3;

	printf("\nDATA:\n");
	data = gsl_matrix_alloc(nx, ny);

	gsl_matrix_set(data, 0, 0, 2.5);
	gsl_matrix_set(data, 0, 1, 2.4);
	gsl_matrix_set(data, 0, 2, 2.1);

	gsl_matrix_set(data, 1, 0, 0.5);
	gsl_matrix_set(data, 1, 1, 0.7);
	gsl_matrix_set(data, 1, 2, 0.7);

	gsl_matrix_set(data, 2, 0, 2.2);
	gsl_matrix_set(data, 2, 1, 2.9);
	gsl_matrix_set(data, 2, 2, 2.6);

	gsl_matrix_set(data, 3, 0, 1.9);
	gsl_matrix_set(data, 3, 1, 2.2);
	gsl_matrix_set(data, 3, 2, 2.8);

	gsl_matrix_set(data, 4, 0, 3.1);
	gsl_matrix_set(data, 4, 1, 3.0);
	gsl_matrix_set(data, 4, 2, 3.3);

	gsl_matrix_set(data, 5, 0, 2.3);
	gsl_matrix_set(data, 5, 1, 2.7);
	gsl_matrix_set(data, 5, 2, 2.3);

	gsl_matrix_set(data, 6, 0, 2.0);
	gsl_matrix_set(data, 6, 1, 1.6);
	gsl_matrix_set(data, 6, 2, 3.0);

	gsl_matrix_set(data, 7, 0, 1.0);
	gsl_matrix_set(data, 7, 1, 1.1);
	gsl_matrix_set(data, 7, 2, 1.1);

	gsl_matrix_set(data, 8, 0, 1.5);
	gsl_matrix_set(data, 8, 1, 1.6);
	gsl_matrix_set(data, 8, 2, 1.2);

	gsl_matrix_set(data, 9, 0, 1.1);
	gsl_matrix_set(data, 9, 1, 0.9);
	gsl_matrix_set(data, 9, 2, 0.8);

	printMatrix(data);


	gsl_vector* means = gsl_vector_calloc(data->size2);

	// DO PCA //
	printf("\nROTATION MATRIX (");
	gsl_matrix* rot = pca_cols(data, means, 2);
	printf("):\n");
	printMatrix(rot);

	// TRANSFORM DATA //
	printf("\nTRANSFORMED DATA:\n");
	gsl_matrix* transformed =  transformData(data, rot, means);
	printMatrix(transformed);

	gsl_matrix_free(data);

	// GET BACK ORIGINAL DATA //
	printf("\nRECONSTRUCTED DATA:\n");
	gsl_matrix* rec_data = reconstructData(transformed, rot, means);
	printMatrix(rec_data);

	gsl_matrix_free(rot);
	gsl_matrix_free(transformed);
	gsl_vector_free(means);

	return(0);
}
*/
