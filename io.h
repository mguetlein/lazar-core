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

#include <string>	
#include <iostream>
#include <sstream>

#include "ServerSocket.h"		// socket library


#ifndef IO_class
#define IO_class

using namespace std;

class Out: public stringstream {
	public:
		Out() {};
		virtual void print() {};
		virtual void print_err() {};
};

class ConsoleOut: public Out {
	public:
		ConsoleOut() {};
		void print();
		void print_err();
};

class SocketOut: public Out {

	private:
		ServerSocket * socket;
		
	public:
		SocketOut(ServerSocket * socket): socket(socket) {} ;
		void print();
		void print_err();
};
#endif
