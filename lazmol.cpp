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


#include "lazmol.h"

// LazMol

LazMol::LazMol() {};

LazMol::LazMol(int nr, string new_descr, string new_smiles, Out * out): line_nr(nr), id(new_descr), smiles(new_smiles), out(out) {

};

int LazMol::get_line_nr() { return(line_nr); };

void LazMol::set_id(string id) { id = id; };

void LazMol::set_inchi(string new_inchi) { inchi = new_inchi; };

string LazMol::get_id() { return(id); };

void LazMol::print() {
	int l = line_nr + 1;
	*out << "line_nr: " << l << "\n";
	*out << "id: " << id << "\n";
	*out << "smiles: '" << smiles << "'\n";
	*out << "inchi: '" << inchi << "'\n";
	out->print();
};

string LazMol::get_smiles() { return(smiles); };

string LazMol::get_inchi() { return(inchi); };

void LazMol::set_output(Out * newout) { out = newout; };

// OBLazMol

OBLazMol::OBLazMol(int nr, string new_descr, string new_smiles, Out * out):

	LazMol(nr,new_descr,new_smiles,out), out(out) {

		static OBMol mol;
		string tmp_inchi;
		OBConversion conv(&cin,&cout);
		conv.SetInAndOutFormats("SMI","INCHI");
		if (!conv.ReadString(&mol,new_smiles)) {
			*out << "\nError reading molecule nr. " << nr <<  endl;
			out->print_err();
		}
		else {
			// don't warn about undefined stereo
			conv.SetOptions("w",OBConversion::OUTOPTIONS);
			string inchi = conv.WriteString(&mol);
			// remove newline
			string::size_type pos = inchi.find_last_not_of("\n");
			if(pos != string::npos) {
				inchi = inchi.substr(0, pos+1);
			}
			this->set_inchi(inchi);
//			cerr << new_smiles << "\t" << inchi << endl;
		}
};

//OBMol * OBLazMol::get_mol_ref() { return(&mol); };
//
//bool OBLazMol::match(OBSmartsPattern * smarts_pattern) {
//	return (smarts_pattern->Match(mol,true));
//};
//
//int OBLazMol::match_freq(OBSmartsPattern * smarts_pattern) {
//	smarts_pattern->Match(mol,false);
//	vector<vector<int> > maplist;
//	maplist = smarts_pattern->GetUMapList();
//	return (maplist.size());
//};
//
//vector<string> OBLazMol::sssr() {
//
//	vector<string> ringset;
//	OBElementTable element_table;
//
//	// identify the smallest set of smallest rings //
//	vector<OBRing*> ringsystems = mol.GetSSSR();
//
//	vector<OBRing*>::iterator cur_ring;
//	for (cur_ring = ringsystems.begin(); cur_ring != ringsystems.end(); cur_ring++) {
//
//		vector<int>::iterator cur_atom;
//		vector<int>::iterator next_atom;
//
//		stringstream smarts;	// stores smarts string
//
//		// convert ring systems to smarts //
//		for(cur_atom = (*cur_ring)->_path.begin();cur_atom != (*cur_ring)->_path.end();cur_atom++) {
//				int j = *cur_atom;
//				next_atom = cur_atom + 1;
//
//				// identify element symbols for atoms //
//				OBAtom atom = *mol.GetAtom(j);
//				int atom_nr = atom.GetAtomicNum();
//				string symbol = element_table.GetSymbol(atom_nr);
//
//				// convert aromatic elements to lowercase //
//				if (atom.IsAromatic()) transform (symbol.begin(), symbol.end(), symbol.begin(), towlower);
//
//				//Add square brackets for nonstandard elements //
//				if (symbol.size() > 1) smarts << "[";
//
//				smarts << symbol;	// append to smarts
//
//				//Add square brackets for nonstandard elements //
//				if (symbol.size() > 1) smarts << "]";
//
//				// add ring closure for the first atom //
//				if ( cur_atom == (*cur_ring)->_path.begin()) smarts << "1";
//
//				else {
//
//					// identify double bonds in nonaromatic rings //
//					if (!(*cur_ring)->IsAromatic()) {
//
//						// only, if we are not at the last atom
//						if (next_atom != (*cur_ring)->_path.end() ) {
//							OBBond bond = *mol.GetBond(*cur_atom,*next_atom);
//							if (bond.IsDouble()) smarts << "=";
//						}
//					}
//				}
//		}
//		smarts << "1";
//		string tmp;
//		smarts >> tmp;
//		ringset.push_back(tmp);
//	}
//
//	vector<string> unique_rings(ringset.size());
//	sort(ringset.begin(),ringset.end());
//	vector<string>::iterator end_it;
//	end_it = unique_copy(ringset.begin(),ringset.end(),unique_rings.begin());
//	unique_rings.erase(end_it,unique_rings.end());
//	return (unique_rings);
//};
