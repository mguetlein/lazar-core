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

float sig_thr=0.9;

//! generate physicochemical properties in the form
//  id \t molwgt \t logP \t psa \t mr
int main(int argc, char *argv[]) {

	int status = 0;
	int c;
	bool t_file = false;
	char * structure_file = NULL;

	typedef MolVect<OBLazMol,OBLinFrag,bool> OBLazMolVect ;

	while ((c = getopt(argc, argv, "s:")) != -1) {
		switch(c) {
		case 's':
		 	structure_file = optarg; 
			t_file = true;
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

	if (status | !t_file) {
		*out << "usage: " << argv[0] << " -s id_and_smiles\n";
		out->print_err();
		return(status);
	}

	OBLazMolVect * structures = new OBLazMolVect(structure_file, out);

	FeatGen<OBLazMol,OBLinFrag,bool> * fragments = new FeatGen<OBLazMol,OBLinFrag,bool>(structures, out);

	fragments->generate_pcprop(out);

	delete fragments;
	delete structures;
	return (0);

}
