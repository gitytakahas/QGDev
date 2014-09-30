#include <vector>
#include <map>
#include <stdexcept>

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TIterator.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TVector.h>
#include <TString.h>

#include "localQGLikelihoodCalculator.h"


/// Constructor
QGLikelihoodCalculator::QGLikelihoodCalculator(const TString& filename){
  if(filename == "" || !this->init(filename)) throw std::runtime_error("Initialization failed: please check filename");
}


/// Initialisation of the QGLikelihoodCalculator
bool QGLikelihoodCalculator::init(const TString& fileName){
  f = new TFile(fileName);
  if(f->IsZombie()) 				return false;
  if(!getBinsFromFile(etaBins, "etaBins", f))	return false;
  if(!getBinsFromFile(ptBinsC, "ptBinsC", f))	return false;
  if(!getBinsFromFile(ptBinsF, "ptBinsF", f))	return false;
  if(!getBinsFromFile(rhoBins, "rhoBins", f))	return false;

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


/// Compute the QGLikelihood, given the pT, eta, rho and likelihood variables vector
float QGLikelihoodCalculator::computeQGLikelihood(float pt, float eta, float rho, std::vector<float> vars){
  if(!isValidRange(pt, rho, eta)) return -1;

  float Q=1., G=1.;
  for(unsigned int varIndex = 0; varIndex < vars.size(); ++varIndex){

    auto qgEntry = findEntry(eta, pt, rho, 0, varIndex);
    if(!qgEntry) return -1;
    int binQ = qgEntry->FindBin(vars[varIndex]);
    float Qi = qgEntry->GetBinContent(binQ)*qgEntry->GetBinWidth(binQ);
    float mQ = qgEntry->GetMean();

    qgEntry = findEntry(eta, pt, rho, 1, varIndex); 
    if(!qgEntry) return -1;
    int binG = qgEntry->FindBin(vars[varIndex]);
    float Gi = qgEntry->GetBinContent(binG)*qgEntry->GetBinWidth(binG);
    float mG = qgEntry->GetMean();

    float epsilon=0;
    float delta=0.000001;
    if(Qi <= epsilon && Gi <= epsilon){
      if(mQ>mG){
	if(vars[varIndex] > mQ){ Qi = 1-delta; Gi = delta;}
	else if(vars[varIndex] < mG){ Qi = delta; Gi = 1-delta;}
      }
      else if(mQ<mG){
	if(vars[varIndex]<mQ) { Qi = 1-delta; Gi = delta;}
	else if(vars[varIndex]>mG){Qi = delta;Gi = 1-delta;}
      }
    } 
    Q*=Qi;
    G*=Gi;	
  }

  if(Q==0) return 0;
  return Q/(Q+G);
}


/// Compute the QGLikelihood using CDF, given the pT, eta, rho and likelihood variables vector
float QGLikelihoodCalculator::computeQGLikelihoodCDF(float pt, float eta, float rho, std::vector<float> vars){
  if(!isValidRange(pt, rho, eta)) return -1;

  float Q=1., G=1.;
  for(unsigned int varIndex = 0; varIndex < vars.size(); ++varIndex){

    auto quarkEntry = findEntry(eta, pt, rho, 0, varIndex);
    auto gluonEntry = findEntry(eta, pt, rho, 1, varIndex);
    if(!quarkEntry || !gluonEntry) return -1;

    float cdfQuark = quarkEntry->Integral(0, quarkEntry->FindBin(vars[varIndex]), "width");
    float cdfGluon = gluonEntry->Integral(0, gluonEntry->FindBin(vars[varIndex]), "width");

    float ccdfQuark = 1.-cdfQuark;
    float ccdfGluon = 1.-cdfGluon;

    float Qi, Gi;
    if(quarkEntry->GetMean() < gluonEntry->GetMean()){
      if(cdfQuark+cdfGluon <= 0){
        Qi = 0.99;
        Gi = 0.01;
      } else if(ccdfQuark+ccdfGluon <= 0){
        Qi = 0.01;
        Gi = 0.99;
      } else {
        Qi = cdfQuark/(cdfQuark+cdfGluon) - 0.5;
        Gi = ccdfGluon/(ccdfQuark+ccdfGluon) - 0.5;
      }
    } else {
      if(cdfQuark+cdfGluon <= 0){
        Qi = 0.01;
        Gi = 0.99;
      } else if(ccdfQuark+ccdfGluon <= 0){
        Qi = 0.99;
        Qi = 0.01;
      } else {
        Qi = ccdfQuark/(ccdfQuark+ccdfGluon) - 0.5;
        Gi = cdfGluon/(cdfQuark+cdfGluon) - 0.5;
      }
    }

    Q*=Qi;
    G*=Gi;
  }

  if(Q+G <= 0) return 0.5;
  return Q/(Q+G);
}


/// Find matching entry for a given eta, pt, rho, qgIndex and varIndex
TH1F* QGLikelihoodCalculator::findEntry(float eta, float pt, float rho, int qgIndex, int varIndex){
  int etaBin, ptBin, rhoBin;
  if(!getBinNumber(etaBins, fabs(eta), etaBin)) 			return nullptr;
  if(!getBinNumber(etaBin == 0? ptBinsC : ptBinsF, pt, ptBin)) 		return nullptr;
  if(!getBinNumber(rhoBins, rho, rhoBin)) 				return nullptr;
  TString histName = (varIndex == 2 ? "axis2" : (varIndex? "ptD" : "mult")) + TString("_") + (qgIndex ? "gluon":"quark")  + TString::Format("_eta-%d_pt-%d_rho-%d", etaBin, ptBin, rhoBin);
  return pdfs[histName];
}


/// Check the valid range of this qg tagger, using the bin vectors
bool QGLikelihoodCalculator::isValidRange(float pt, float rho, float eta){
  if(pt < ptBinsC.front()) 		return false;
  if(pt > ptBinsC.back()) 		return false;
  if(rho < rhoBins.front()) 		return false;
  if(rho > rhoBins.back()) 		return false;
  if(fabs(eta) < etaBins.front()) 	return false;
  if(fabs(eta) > etaBins.back()) 	return false;
  return true;
}


/// Translates the TVector with the bins to std::vector
bool QGLikelihoodCalculator::getBinsFromFile(std::vector<float>& bins, const TString& name, TFile* f){
  TVectorT<float> *tbins = nullptr;
  f->GetObject(name, tbins);
  if(!tbins) return false;
  for(int i = 0; i < tbins->GetNoElements(); ++i) bins.push_back((*tbins)[i]);
  return true;
}


/// Find the bin number for a value, given the bin vector
bool QGLikelihoodCalculator::getBinNumber(std::vector<float>& bins, float value, int& bin){
  if(value < bins.front() || value > bins.back()) return false;
  auto binUp = bins.begin() + 1;
  while(value > *binUp) ++binUp;
  bin = binUp - bins.begin() - 1;
  return true;
}


/// Destroy the QGLikelihoodCalculator
QGLikelihoodCalculator::~QGLikelihoodCalculator(){
  for(auto& pdf : pdfs) delete pdf.second;
  delete f;
}
