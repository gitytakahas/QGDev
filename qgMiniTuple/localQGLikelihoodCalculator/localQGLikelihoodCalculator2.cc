#include <vector>
#include <map>
#include <stdexcept>
#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TIterator.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TVector.h>
#include <TString.h>

#include "localQGLikelihoodCalculator2.h"


/// Constructor
QGLikelihoodCalculator2::QGLikelihoodCalculator2(const TString& filename){
  if(filename == "" || !this->init(filename)) throw std::runtime_error("Initialization failed: please check filename");
}


/// Initialisation of the QGLikelihoodCalculator2
bool QGLikelihoodCalculator2::init(const TString& fileName){
  f = new TFile(fileName);
  if(f->IsZombie()) 				return false;
  if(!(this->getBinsFromFile(etaBins, "etaBins")))	return false;
  if(!(this->getBinsFromFile(ptBinsC, "ptBinsC")))	return false;

  TList *keys = f->GetListOfKeys();
  if(!keys) return false;

  TIter nextdir(keys);
  TKey *keydir;
  while((keydir = (TKey*) nextdir())){
    if(!keydir->IsFolder()) continue;
    TDirectory *dir = (TDirectory*) keydir->ReadObj() ;
    TIter nexthist(dir->GetListOfKeys());
    TKey *keyhist;
    while((keyhist = (TKey*) nexthist())){
      pdfs[keyhist->GetName()] = (TH1F*) keyhist->ReadObj();
    }
  }
  return true;
}


/// Compute the QGLikelihood, given the pT, eta and likelihood variables vector
float QGLikelihoodCalculator2::computeQGLikelihood(float pt, float eta, std::vector<float> vars){
  eta = fabs(eta);
  if(!isValidRange(pt, eta)) return -1;

  float Q=1., G=1.;
  //std::cout << "pt: " << pt << std::endl;
  for(unsigned int varIndex = 0; varIndex < vars.size(); ++varIndex){
  
    if(vars[varIndex] < -0.5) continue; //use to inspect variables separately (i.e. skip if feeding -1)

    auto quarkEntry = findEntry(eta, pt, 0, varIndex);
    auto gluonEntry = findEntry(eta, pt, 1, varIndex);
    //std::cout << quarkEntry->GetName() << std::endl;
    //std::cout << gluonEntry->GetName() << std::endl;
    if(!quarkEntry && !gluonEntry) return -2;

    int binQ = quarkEntry->FindBin(vars[varIndex]);
    float Qi = quarkEntry->GetBinContent(binQ)/quarkEntry->Integral("width");
    //float Qw = quarkEntry->GetBinWidth(binQ);

    int binG = gluonEntry->FindBin(vars[varIndex]);
    float Gi = gluonEntry->GetBinContent(binG)/gluonEntry->Integral("width");
    //float Gw = gluonEntry->GetBinWidth(binG);

    //if(Qi <= 0 || Gi <= 0){
    //  bool gluonAboveQuark = (quarkEntry->GetMean() < gluonEntry->GetMean());

    //  if(quarkEntry->Integral(0, binQ) + gluonEntry->Integral(0, binG) <=0){
    //    if(gluonAboveQuark){ Qi = 0.999; Gi = 0.001;} else { Qi = 0.001; Gi = 0.999;}
    //  } else if(quarkEntry->Integral(binQ, quarkEntry->GetNbinsX() + 1) + gluonEntry->Integral(binG, gluonEntry->GetNbinsX() + 1) <=0){
    //    if(gluonAboveQuark){ Qi = 0.001; Gi = 0.999;} else { Qi = 0.999; Gi = 0.001;}
    //  } else {
    //    int q = 1, g = 1;
    //    while(Qi <= 0){
    //      Qi += quarkEntry->GetBinContent(binQ+q) + quarkEntry->GetBinContent(binQ-q);
    //      Qw += quarkEntry->GetBinWidth(binQ+q) + quarkEntry->GetBinWidth(binQ-q);
    //      ++q;
    //    }
    //    while(Gi <= 0){
    //      Gi += gluonEntry->GetBinContent(binG+g) + gluonEntry->GetBinContent(binG-g);
    //      Gw += gluonEntry->GetBinWidth(binG+g) + gluonEntry->GetBinWidth(binG-g);
    //      ++g;
    //    }
    //  }
    //} 
    Q*=Qi;
    G*=Gi;
    //Q*=Qi/Qw;
    //G*=Gi/Gw;
  }

  if(Q==0) return 0;
  return Q/(Q+G);
}




/// Find matching entry for a given eta, pt, qgIndex and varIndex
TH1F* QGLikelihoodCalculator2::findEntry(float eta, float pt, int qgIndex, int varIndex){
  int etaBin, ptBin;
  if(!getBinNumber(etaBins, fabs(eta), etaBin)) 			return nullptr;
  if(!getBinNumber(ptBinsC, pt, ptBin)) 		return nullptr;
  TString histName = (varIndex == 2 ? "axis2" : (varIndex? "ptD" : "mult")) + TString("_") + (qgIndex ? "gluon":"quark")  + TString::Format("_eta%d_pt%d", etaBin, ptBin);
  return pdfs[histName];
}


/// Check the valid range of this qg tagger, using the bin vectors
bool QGLikelihoodCalculator2::isValidRange(float pt, float eta){
  if(pt < ptBinsC.front()) 		return false;
  if(pt > ptBinsC.back()) 		return false;
  if(fabs(eta) < etaBins.front()) 	return false;
  if(fabs(eta) > etaBins.back()) 	return false;
  return true;
}


/// Translates the TVector with the bins to std::vector
bool QGLikelihoodCalculator2::getBinsFromFile(std::vector<float>& bins, const TString& name ) {
  TVectorT<float> *tbins = nullptr;
  f->GetObject(name, tbins);
  if(!tbins) return false;
  for(int i = 0; i < tbins->GetNoElements(); ++i) bins.push_back((*tbins)[i]);
  return true;
}


/// Find the bin number for a value, given the bin vector
bool QGLikelihoodCalculator2::getBinNumber(std::vector<float>& bins, float value, int& bin){
  if(value < bins.front() || value > bins.back()) return false;
  auto binUp = bins.begin() + 1;
  while(value > *binUp) ++binUp;
  bin = binUp - bins.begin() - 1;
  return true;
}


/// Destroy the QGLikelihoodCalculator2
QGLikelihoodCalculator2::~QGLikelihoodCalculator2(){
  for(auto& pdf : pdfs) delete pdf.second;
  delete f;
}
