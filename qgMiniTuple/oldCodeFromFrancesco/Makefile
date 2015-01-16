CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lTMVA




all: fitQGVars drawRoCs


fitQGVars: fitQGVars.cpp localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o fitQGVars fitQGVars.cpp localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

drawRoCs: drawRoCs.cpp localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o drawRoCs drawRoCs.cpp localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

scanQGL: scanQGL.cpp localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o scanQGL scanQGL.cpp localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

checkQGL: checkQGL.cpp localQGLikelihoodCalculator2.o localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o checkQGL checkQGL.cpp localQGLikelihoodCalculator2.o localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

makePDFs: makePDFs.cpp localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o makePDFs makePDFs.cpp localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

makePDFs_onlyPt: makePDFs_onlyPt.cpp localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o makePDFs_onlyPt makePDFs_onlyPt.cpp localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

comparePDFs: comparePDFs.cpp localQGLikelihoodCalculator.o
	$(CC) -Wall $(INCLUDES) -o comparePDFs comparePDFs.cpp localQGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)



localQGLikelihoodCalculator.o: ../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc
	$(CC) $(CFLAGS) -c ../localQGLikelihoodCalculator/localQGLikelihoodCalculator.cc $(ROOTFLAG) $(EXTRALIBS)

localQGLikelihoodCalculator2.o: ../localQGLikelihoodCalculator/localQGLikelihoodCalculator2.cc
	$(CC) $(CFLAGS) -c ../localQGLikelihoodCalculator/localQGLikelihoodCalculator2.cc $(ROOTFLAG) $(EXTRALIBS)
