QGDev
=====

Repository for development of 13TeV quark-gluon tagger

Making a new qgMiniTuple
------------------------
* Use CMSSW_7_2_X or higher
* Extra variables can be added in QGDev/qgMiniTuple/plugins/qgMiniTuple.cc
* cmsRun executable is QGDev/qgMiniTuple/test/qgMiniTuple_cfg.py or QGDev/qgMiniTuple/test/qgMiniTupleForMiniAOD_cfg.py

Macros
------
Compile them with **g++ -O3 -I \`root-config --incdir\` -o $1 $1.C \`root-config --libs\` s-td=c++11** with $1 the executable
* **createPDF**: creates the ROOT files with the pdf's
* **plotDistributions**: plots distributions of axis2, ptD, mult and qg-likelihood
* **plotROC**: plots ROC curves


Create new QGL_v*.db file
-------------------------
* Use CMSSW_7_4_X or higher
* cms git-addpkg JetMETCorrections/Modules
* go to JetMETCorrections/Modules/test and copy the ROOT files created by createPDF.C
* update QGLikelihoodDBWriter_cfg.py with new version number
