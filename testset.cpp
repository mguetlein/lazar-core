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

//! generate all linear fragments
int main(int argc, char *argv[]) {

	int status = 0;
	int c;
	bool s_file = false;
	bool p_set = false;
	bool i_set = false;
	int i = 0;

	char * structure_file = NULL;
	int percentage = 0;

	typedef MolVect<OBLazMol,OBLinFrag,bool> OBLazMolVect ;

	while ((c = getopt(argc, argv, "s:p:i:")) != -1) {
		switch(c) {
		case 's':
		 	structure_file = optarg; 
			s_file = true;
		    break;
		case 'p':
			percentage = atoi(optarg);
			if ((percentage > 0) && (percentage < 100)) p_set = true;
		    break;
		case 'i':
			i = atoi(optarg);
			if (i>=0) i_set = true;
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

	if (status | !s_file | !p_set ) {
		*out << "Puts p\% of database files into a validation file." << endl;
		*out << "usage: " << argv[0] << " -s input_database_file -p percentage [-i start index]" << endl;
		*out << "When -i is given, the next p percent of consecutive compounds are used." << endl;
		*out << "When -i is not given, p percent are drawn at random." << endl;
		out->print_err();
		return(status);
	}

	OBLazMolVect * structures = new OBLazMolVect(structure_file, out);
	FeatGen<OBLazMol,OBLinFrag,bool> * testset = new FeatGen<OBLazMol,OBLinFrag,bool>(structures, out);

	if (!i_set)
		testset->generate_testset(percentage,  out);
	else {
		testset->generate_testset(percentage, i, out);    
	}

	delete structures;
	delete testset;
	return (0);

}
