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


#include "feature.h"
#include "stats.h"
#include <limits>

extern float sig_thr;

// Feature

string Feat::get_name() { return (name); };

void Feat::set_name(string newname) { name = newname; };

void Feat::print(Out * out) {
	*out << name;
	out->print();
};

// ClassFeat

void ClassFeat::precompute_significance(string act, float n_a, float n_i, float f_a, float f_i) { // AM: determine significance

	fa[act]=f_a;
	fi[act]=f_i;
	na[act]=n_a;
	ni[act]=n_i;

	//pre-computing 4 significance values for all features for loo classification
	//the right value for each feature depends on two criteria:
	// - is the test structure active (reduces either n_a, or n_i)
	// - does the feature occur in the test structure (if yes then reduce either f_a or f_i)

	cur_str_active = true;
	cur_feat_occurs = true;
	precompute_significance(act);

	cur_str_active = true;
	cur_feat_occurs = false;
	precompute_significance(act);

	cur_str_active = false;
	cur_feat_occurs = true;
	precompute_significance(act);

	cur_str_active = false;
	cur_feat_occurs = false;
	precompute_significance(act);
}

bool ClassFeat::cur_str_active = false;

void ClassFeat::set_cur_feat_occurs(bool feat_occurs){
	cur_feat_occurs = feat_occurs;
}

string ClassFeat::get_map_key(string act)
{
	if(cur_str_active)
		if(cur_feat_occurs)
			return act.append("_active_occurs");
		else
			return act.append("_active");
	else
		if (cur_feat_occurs)
			return act.append("_occurs");
		else
			return act;
}

void ClassFeat::precompute_significance(string act) { // AM: determine significance

	float n_a = get_na(act);
	float n_i = get_ni(act);
	float f_a = get_fa(act);
	float f_i = get_fi(act);

	float ea;
	float ei;
	float chisq;

	// fragments with total frequency of 1 are not informative, but confounding
	if (f_a>=0 && f_i>=0 && (f_a+f_i)>1) {

		ea = n_a*(f_a+f_i)/(n_a+n_i); // expected actives
		ei = n_i*(f_a+f_i)/(n_a+n_i);  // expected inactives

		// calculate chi square with Yate's correction
		// i.e. reduce observed frequencies by 0.5
		chisq = (f_a-ea-0.5)*(f_a-ea-0.5)/ea + (f_i-ei-0.5)*(f_i-ei-0.5)/ei;

	}
	else
		chisq = 0;

	significance[get_map_key(act)]=chisq;
    p[get_map_key(act)]=calc_p(act);
	float cur_chisq = 0;
	float cur_ea;
	float cur_ei;
	int i;

	if (n_a > n_i) {	// exchange na and ni
		float tmp  = n_a;
		n_a = n_i;
		n_i = tmp;
	}
	// find minimum frequency for statistical significance
	for (i=1; cur_chisq < 3.84; i++) {
		cur_ea = n_a*i/(n_a+n_i);
		cur_ei = n_i*i/(n_a+n_i);
		cur_chisq = (i-cur_ea-0.5)*(i-cur_ea-0.5)/cur_ea + (cur_ei+0.5)*(cur_ei+0.5)/cur_ei ;
	}
	if (f_a<0 || f_i<0 || (f_a+f_i) < i)
		too_infrequent[get_map_key(act)] = true;
	else
		too_infrequent[get_map_key(act)] = false;
};


void ClassFeat::determine_significance(string act, float n_a, float n_i, vector<bool> * activities) { // AM: determine significance

	vector<bool>::iterator a;

	float ea;
	float ei;
	float chisq;
	float f_a=0;
	float f_i=0;

	// AM: cell sums
	for (a = activities->begin(); a != activities->end(); a++) {

		if (*a)
			f_a++;
		else if (!*a)
			f_i++;

	}

	// fragments with total frequency of 1 are not informative, but confounding
	if ((f_a+f_i)>1) {

		ea = n_a*(f_a+f_i)/(n_a+n_i); // expected actives
		ei = n_i*(f_a+f_i)/(n_a+n_i);  // expected inactives

		// calculate chi square with Yate's correction
		// i.e. reduce observed frequencies by 0.5
		chisq = (f_a-ea-0.5)*(f_a-ea-0.5)/ea + (f_i-ei-0.5)*(f_i-ei-0.5)/ei;

	}

	else
		chisq = 0;

	fa[act]=f_a;
	fi[act]=f_i;
	na[act]=n_a;
	ni[act]=n_i;
	significance[act]=chisq;
    p[act]=calc_p(act);
	float cur_chisq = 0;
	float cur_ea;
	float cur_ei;
	int i;

	if (n_a > n_i) {	// exchange na and ni
		float tmp  = n_a;
		n_a = n_i;
		n_i = tmp;
	}
	// find minimum frequency for statistical significance
	for (i=1; cur_chisq < 3.84; i++) {
		cur_ea = n_a*i/(n_a+n_i);
		cur_ei = n_i*i/(n_a+n_i);
		cur_chisq = (i-cur_ea-0.5)*(i-cur_ea-0.5)/cur_ea + (cur_ei+0.5)*(cur_ei+0.5)/cur_ei ;
	}
	if ((f_a+f_i) < i)
		too_infrequent[act] = true;
	else
		too_infrequent[act] = false;

};

void ClassFeat::print(string act,Out * out) {

	float na = get_na(act);
	float ni = get_ni(act);
	float fa = get_fa(act);
	float fi = get_fi(act);

	*out << "    - smarts: '" << this->get_name() << "'\n";
	*out << "      p_chisq: '" << get_p(act) << "'\n";
	*out << "      chisq: '" << get_significance(act) << "'\n";
	*out << "      fa: '" << fa/na << "'\n";
	*out << "      fi: '" << fi/ni << "'\n";
	*out << "      property: '";

	if (fa/(fa+fi) > na/(na+ni)) {
		*out << "activating";
	}
	else {
		*out << "deactivating";
	}
	*out << "'\n";
	out->print();
};

float ClassFeat::get_sig_limit() { return(3.84); };

float ClassFeat::get_p_limit() { return(0.95); };

float ClassFeat::calc_p(string act) {

	return gsl_cdf_chisq_P(significance[get_map_key(act)], 1);
};

float ClassFeat::get_p(string act) {

	return p[get_map_key(act)];
};

float ClassFeat::get_na(string act) {

	if (cur_str_active)
		return na[act]-1;
	else
		return na[act];
};

float ClassFeat::get_ni(string act) {

	if (cur_str_active)
		return ni[act];
	else
		return ni[act]-1;
};

float ClassFeat::get_fa(string act) {

	if (cur_str_active && cur_feat_occurs)
		return fa[act]-1;
	else
		return fa[act];
};

float ClassFeat::get_fi(string act) {

	if (!cur_str_active && cur_feat_occurs)
		return fi[act]-1;
	else
		return fi[act];
};

void ClassFeat::print_specifics(string act, Out* out) {

	float na = get_na(act);
	float ni = get_ni(act);
	float fa = get_fa(act);
	float fi = get_fi(act);

	*out << this->get_name() << "\t" << get_p(act) << "\t" << get_significance(act)
		<< "\t" << (int)fa << "\t" << (int)fi << "\t" << (int)na << "\t" << (int)ni << "\t";
	out->print();

	if (fa/(fa+fi) >na/(na+ni))
			*out << "a";
	else
			*out << "i";

	out->print();
}





// RegrFeat

void RegrFeat::determine_significance(string act, vector<float> all_activities, vector<float> feat_activities) {

	// Kolmogorov-Smirnov Test
	// numerical recipies in C pp 626, extended version with better sensitivity at the ends
	unsigned int j1=0, j2=0;
	float d,d1,d2,d_1,d_2,dt1,dt2,en1,en2,en,fn1=0,fn2=0,alam;
	d1 = d2 = d_1 = d_2 = 0.0;

	sort(feat_activities.begin(),feat_activities.end());
	sort(all_activities.begin(), all_activities.end());

	en1 = all_activities.size();
	en2 = feat_activities.size();

	while (j1 < en1 && j2 < en2) {
		if ((!isnan(d1=all_activities[j1])) && (!isnan(d2=feat_activities[j2]))) {
			if (d1 <= d2) { // next step is in all_activities
				j1++;
				fn1=j1/en1;
			}

			if (d2 <= d1) { // next step is in feat_activities
				j2++;
				fn2=j2/en2;
			}

		}
		else {
			if (isnan(d1)) {
				j1++;
			}
			if (isnan(d2)){
				j2++;
			}
			continue;
		}
		dt1=fn2-fn1;
		dt2=fn1-fn2;
		if (dt1 > d_1) d_1=dt1;
		if (dt2 > d_2) d_2=dt2;
	}
	d = d_1 + d_2;
	en=sqrt(en1*en2/(en1+en2));
	alam=(en+0.155+0.24/en)*d;

	/*
    // KS probability function
	a2 = -2.0 *alam*alam;
	p = 0;	// if probability function does not converge
	fac2 = 4.0*alam*alam;
	for (j=1;j<=100;j++) {
		term = fac*((fac2*j*j)-1.0)*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= 0.001*termbf || fabs(term) <= 1.0e-8*sum) {
			p=1-sum;
			break;
		}
		termbf=fabs(term);
	}
    */

	// compute and compare medians to determine activation property
	float aasum, aamedian, aamean, aavar, aadev, aaskew, aakurt;
	computeStats(all_activities.begin(), all_activities.end(), aasum, aamedian, aamean, aavar, aadev, aaskew, aakurt);
	float fasum, famedian, famean, favar, fadev, faskew, fakurt;
	computeStats(feat_activities.begin(), feat_activities.end(), fasum, famedian, famean, favar, fadev, faskew, fakurt);

	median[act]=famedian;
	global_median[act]=aamedian;
	significance[act]=alam;
    p[act]=calc_p(act);
};


void RegrFeat::determine_significance(string act, float global_med, vector<float> * activities) {

	// perform a simple sign test

	vector<float> tmp;
	vector<float>::iterator a;

	int f_l = 0;
	int f_s = 0;
	int n = 0;
	float sig =0;
	float med;

	// AM: count: f_l structures below
	//            f_s structures above     median
	for (a = activities->begin(); a != activities->end(); a++) {
		n++;
		if (*a > global_med) {
			f_l++;
			tmp.push_back(*a);
		}
		else if (*a < global_med) {	// ignore zero differences
			f_s++;
			tmp.push_back(*a);
		}
	}

	int k = min(f_l,f_s);

	// calculate cumulative distribution function

	for (int i = 0; i <= k; i++) {
		sig= sig + gsl_ran_binomial_pdf(i,0.5,n);
	}

	if (sig < 0.5)
		sig = 1-2*sig;
	else
		sig =0;

	if (tmp.size()>0) {
		sort(tmp.begin(),tmp.end());
		med = ( *(tmp.begin() + (tmp.size()/2)) + *(tmp.begin() + (tmp.size()/2-1+tmp.size()%2)) )/2;

		if (med == global_med) // set p to 0 for the pathologic case that global_median = min|max(all_values)
			sig=0;
	}
	else {
		med = numeric_limits<float>::quiet_NaN();
		sig = 0;
	}
	significance[act]=sig;
    p[act]=calc_p(act);
	median[act]=med;
	global_median[act]=global_med;
	feature_larger[act]=f_l;
	feature_smaller[act]=f_s;
	nr[act]=n;
    };


void RegrFeat::print_header(Out * out) {
	*out << "Feature\tp_sign_test\tmedian\tmedian_trainset\tlarger\tsmaller\tn\tactivating/inactivating\n";
	out->print();
} ;

void RegrFeat::print(string act,Out * out) {

		*out << "    - smarts: '" << this->get_name() << "'\n";
		*out << "      p_ks: '" << this->p[act] << "'\n";
		*out << "      property: '";
		if (median[act]>global_median[act]) {
			*out << "deactivating";
		}
		else {
			*out << "activating";
		}
		*out << "'\n";
		out->print();
};

//float RegrFeat::get_sig_limit() { return(sig_thr); };

float RegrFeat::get_p_limit() { return(sig_thr); };

float RegrFeat::get_global_median(string act) { return(global_median[act]); };

float RegrFeat::get_median(string act) { return(median[act]); };

float RegrFeat::get_p(string act) { return(p[act]); };

float RegrFeat::calc_p(string act) {

	float p,a2,fac=2,sum=0,term,termbf=0,fac2;
	// KS probability function
	a2 = -2.0 *significance[act]*significance[act];
	p = 0;	// if probability function does not converge
	fac2 = 4.0*significance[act]*significance[act];
	for (int j=1;j<=100;j++) {
		term = fac*((fac2*j*j)-1.0)*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= 0.001*termbf || fabs(term) <= 1.0e-8*sum) {
			p=1-sum;
			break;
		}
		termbf=fabs(term);
	}
    return (p);

};

void RegrFeat::print_specifics(string act, Out* out) {
		float p = this->p[act];
		*out << this->get_name() << "\t" << p << "\t";
		out->print();
}

// OBSmartsFrag

OBSmartsFrag::OBSmartsFrag(string smarts): Feat(smarts) {
	smarts_pattern.Init(smarts);
};

//LinFrag

LinFrag::LinFrag(string smarts, bool split_bonds) {

	fragment.clear();

	if (split_bonds) {
		size_t startpos = 0, endpos = 0;

		// split at bonds
		startpos = smarts.find_first_not_of("-=#:|",startpos);
		endpos   = smarts.find_first_of("-=#:|",startpos);

		while (endpos <= smarts.size() && startpos <= smarts.size()) {
			fragment.push_back(smarts.substr(startpos,endpos-startpos)); // add atom
			fragment.push_back(smarts.substr(endpos,1));	// add bond
			startpos = smarts.find_first_not_of("-=#:",endpos);
			endpos   = smarts.find_first_of("-=#:",startpos);
		}


		fragment.push_back(smarts.substr(startpos,smarts.size()-startpos));	 //add last atom
	}

	else { fragment.push_back(smarts); }

};

LinFrag::LinFrag(vector<string> frag): fragment(frag) {};

vector<string> LinFrag::get_fragment() { return (fragment); };

int LinFrag::size() { return ( fragment.size() ); };

// refinement

string LinFrag::canonify() {
	string sma = "";
	string rev_sma = "";
	string name;
	// create smarts string
	for (vector<string>::iterator iter = fragment.begin(); iter != fragment.end(); ++iter) {
		sma = sma + *iter;
	}
	// create reverse smarts string
	for (vector<string>::reverse_iterator iter = fragment.rbegin(); iter != fragment.rend(); ++iter) {
		rev_sma = rev_sma + *iter;
	}
	if (sma < rev_sma) {
		reverse(fragment.begin(),fragment.end());
		sma = rev_sma;
	}

	return(sma);
};

void LinFrag::insert(string element) {
	fragment.insert(fragment.begin(),element);
};

void LinFrag::expand(string element) {
	fragment.push_back(element);
};

void LinFrag::rev() {
	reverse(fragment.begin(),fragment.end());
}

string LinFrag::first_atom() {
	return ( fragment[0] );
};

string LinFrag::first_bond() {
	return ( fragment[1] );
};

string LinFrag::last_atom() {
	vector<string>::iterator last_atom = fragment.end();
	last_atom--;
	return (*last_atom);
};

string LinFrag::last_bond() {
	vector<string>::iterator last_bond = fragment.end();
	last_bond--;
	last_bond--;
	return (*last_bond);
};

bool LinFrag::more_specific(LinFrag * g) {

	vector<string>::iterator s_result;

	vector<string> g_frag = g->get_fragment();
	s_result = search(fragment.begin(),fragment.end(),g_frag.begin(),g_frag.end());
	if (s_result == fragment.end()) {	// try reverse fragment
		s_result = search(fragment.begin(),fragment.end(),g_frag.rbegin(),g_frag.rend());
	}

	return(s_result != fragment.end());

};

// rex generation

string LinFrag::init_wildcard() {
	vector<string>::iterator it;
	fragment[1] = '*';
	string newname;
	for (vector<string>::iterator iter = fragment.begin(); iter != fragment.end(); ++iter) {
		newname = newname + *iter;
	}
	return(newname);
}
