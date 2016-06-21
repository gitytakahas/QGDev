CPP=g++
CFLAGS=-O3 -I `root-config --incdir` -std=c++11 `root-config --libs`
CFLAGS += -ggdb

.PHONY: all
all: createPDF plotLikelihood plotDistributions plotProfile plotProfile2D plotROC

## createPDF:  createPDF.C treeLooper.h
## 	g++ -O3 -I `root-config --incdir` -o $@ $< `root-config --libs` -std=c++11

%:  %.C treeLooper.h
	$(CPP) $(CFLAGS) -o $@ $< 
