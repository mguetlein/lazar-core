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


#include "feature-generation.h"
#include <getopt.h>

using namespace std;

float sig_thr=0;

//! match SMARTS features
int main(int argc, char *argv[]) {

	int status = 0;
	bool t_file, s_file = false;
	char * structure_file = NULL;
	char * smarts_file = NULL;
	typedef MolVect<OBLazMol,OBLinFrag,bool> OBLazMolVect ;

	int c;
	while ((c = getopt(argc, argv, "s:f:")) != -1) {
		switch(c) {
		case 's':
		 	structure_file = optarg; 
			t_file = true;
		    break;
		case 'f':
		    smarts_file = optarg;
			s_file = true;
		    break;
		case ':':
			status = 1;
		    break;
		 case '?':
			status = 1;
		    break;
		}
	}

	ConsoleOut * out;
	out = new ConsoleOut();

	if (status | !t_file | !s_file) {
		fprintf(stderr, "usage: %s -s training_structures -f smarts_file\n",argv[0]);
		return(status);
	}

	OBLazMolVect * structures = new OBLazMolVect (structure_file,out);

	FeatGen<OBLazMol,OBLinFrag,bool> * fragments = new FeatGen<OBLazMol,OBLinFrag,bool>(smarts_file,structures,out);
	fragments->alphabet2level();
	fragments->match_level(true);

	return (0);
}
