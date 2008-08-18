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


#include "activity-db.h"
#include <getopt.h>

using namespace std;

float sig_thr=0.9;
bool quantitative = 0;

//! filter features according to their chisquare value
//! option -n can be used to print chisq, fa, fi
int main(int argc, char *argv[]) {

	int status = 0;
	int c;
	bool s_file=false;
	bool t_file=false;
	bool f_file=false;
	bool print_all = false;
	char * train_file = NULL;
	char * smi_file = NULL;
	char * feature_file = NULL;
	char * smarts = NULL;
	float limit = 0;

	while ((c = getopt(argc, argv, "s:t:f:l:m:n")) != -1) {
		switch(c) {
		case 's':
		 	smi_file = optarg; 
			s_file = true;
		    break;
		case 't':
		 	train_file = optarg; 
			t_file = true;
		    break;
		case 'f':
		 	feature_file = optarg; 
			f_file = true;
		    break;
		case 'l':
		 	limit = atof(optarg); 
		    break;
		case 'm':
			smarts = optarg;
		    break;
		case 'n':
			print_all = true;
		    break;
		case ':':
			status = 1;
		    break;
		 case '?':
			status = 1;
		    break;
		 default:
			status = 1;
		    break;
		}
	}

	ConsoleOut * out;
	out = new ConsoleOut();

	if (status | !s_file | !t_file | !f_file) {
		fprintf(stderr, "usage: %s -s structures -t training_set -f feature_set [-l min_chisq] [-m smarts] [-n]\n",argv[0]);
		return(status);
	}
	
	ActMolVect<OBLazMol,ClassFeat,bool> * train_set;
	train_set = new ActMolVect<OBLazMol,ClassFeat,bool>(train_file,feature_file,smi_file,out);

	if (print_all)
		train_set->print_sorted_features(limit,smarts);
	else
		train_set->print_sig_features(limit,smarts);

	return (0);
}
