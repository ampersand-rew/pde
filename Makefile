# COMPHY: path of COMPHY directory
COMPHY	=	$(HOME)/comphy
FC	=	gfortran
#LFLAGS 	= 	-O -L$(COMPHY)/lib -lnumer
LFLAGS 	= 	-g -Llib -lnumer
FFLAGS 	= 	-c -g

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) 
CXXFLAGS  += $(ROOTCFLAGS)
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3
GXXd	   = /usr/bin/g++ -Wall -g

ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
P5640FLAGS  = -L${P5640LIB}/lib -lP5640  -I${P5640LIB}
GSLFLAGS    = -I${EBROOTGSL}/include/gsl  -I/usr/include/gsl -lgsl -lgslcblas




all: p1


p1: p1.C
	g++ -g -Wall -o p1 p1.C $(ROOTFLAGS) $(GSLFLAGS)

clean:
	rm -f p1 *.o *.so *.pcm *.d *~
