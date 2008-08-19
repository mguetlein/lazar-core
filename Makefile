# Edit according to your needs
SHELL=/bin/bash

PROGRAM = lazar 
FEAT_GEN = rex linfrag smarts-features testset
TOOLS = chisq-filter pcprop
INSTALLDIR = /usr/local/bin

OBJ = feature.o lazmol.o io.o utils.o rutils.o
HEADERS = lazmolvect.h feature.h lazmol.h io.h ServerSocket.h Socket.h utils.h feature-generation.h rutils.h
SERVER_OBJ = ServerSocket.o Socket.o
OBJ += $(SERVER_OBJ)

CC            = g++
CXXFLAGS      = -g -O2 -I../openbabel/include/openbabel-2.0/ -I../R/include/ -Wall
LIBS	      = -lm -ldl -lopenbabel -lgsl -lgslcblas -lR
LDFLAGS       = -L../openbabel/lib -L../R/lib
RPATH         = -Wl,-rpath=../openbabel/lib:../R/lib

.PHONY:
all: $(PROGRAM) $(FEAT_GEN) $(TOOLS)

.PHONY:
doc: $(OBJ) $(HEADERS) $(PROGRAM) $(FEAT_GEN) $(TOOLS) Doxyfile
	doxygen Doxyfile

lazar: $(OBJ)  lazar.o 
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o lazar $(OBJ)  lazar.o 

linfrag: $(OBJ) linfrag.o
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o linfrag $(OBJ) linfrag.o

testset: $(OBJ) testset.o
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o testset $(OBJ) testset.o

pcprop: $(OBJ) pcprop.o
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o pcprop $(OBJ) pcprop.o

rex: $(OBJ) rex.o
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o rex $(OBJ) rex.o

chisq-filter: $(OBJ)  chisq-filter.o 
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o chisq-filter $(OBJ)  chisq-filter.o 

smarts-features: $(OBJ)  smarts-features.o 
	$(CC) $(CXXFLAGS) $(INCLUDE) $(LIBS) $(LDFLAGS) $(RPATH) -o smarts-features $(OBJ)  smarts-features.o 

chisq-filter.o: $(HEADERS) activity-db.h feature-db.h 

rex.o: $(HEADERS) feature-generation.h

linfrag.o: $(HEADERS) feature-generation.h

testset.o: $(HEADERS)

pcprop.o: $(HEADERS)

smarts-features.o: $(HEADERS) feature-generation.h

lazar.o: $(HEADERS) predictor.h model.h activity-db.h feature-db.h feature-generation.h

lazmol.o: lazmol.h

feature.o: feature.h 

io.o: io.cpp io.h $(SERVER_OBJ)

rutils.o: rutils.h

utils.o: utils.h

.PHONY:
clean:
	-rm -rf *.o $(PROGRAM) $(TOOLS) $(FEAT_GEN)
