QGDev
=====

Repository for development of 13TeV quark-gluon tagger

Making a new qgMiniTuple
------------------------
* Use CMSSW_7 or higher
* If CMSSW-version of tagger is used, CMSSW_7_1_0 or higher is needed 
* Extra variables can be added in QGDev/qgMiniTuple/plugins/qgMiniTuple.cc
* cmsRun executable is QGDev/qgMiniTuple/test/qgMiniTuple_cfg.py or QGDev/qgMiniTuple/test/qgMiniTupleForMiniAOD_cfg.py
* For every QGL_v*.db file in QGDev/qgMiniTupled/data there is a cfi fragment in QGDev/qgMiniTuple/python (for CMSSW-version of tagger)
* When running on the grid with crab: if the needed QGL_v*.db file is not yet uploaded to the conditions database, use QGDev/qgMiniTuple/crab/crabCreateWithDB.sh instead of crab -create

Macros
------
Compile them with **g++ -O3 -I \`root-config --incdir\` -o $1 $1.C \`root-config --libs\` s-td=c++0x** with $1 the executable
* **createPDF**: creates the ROOT files with the pdf's
* **plotDistributions**: plots distributions of axis2, ptD, mult and qg-likelihood
* **plotROC**: plots ROC curves


Create new QGL_v*.db file
-------------------------
* Use CMSSW_7_1_X or higher
* cms git-addpkg JetMETCorrections/Modules
* replace JetMETCorrections/Modules/plugins/QGLikelihoodDBWriter.cc with https://raw.githubusercontent.com/UAEDF-tomc/cmssw/updateQGLikelihoodDBWriter/JetMETCorrections/Modules/plugins/QGLikelihoodDBWriter.cc (pull request is made into CMSSW_7_2_X)
* replace JetMETCorrections/Modules/test/QGLikelihoodDBWriter_cfg.py with https://raw.githubusercontent.com/UAEDF-tomc/cmssw/updateQGLikelihoodDBWriter/JetMETCorrections/Modules/test/QGLikelihoodDBWriter_cfg.py
* update QGLikelihoodDBWriter_cfg.py with new version number and add new jet collections if needed
