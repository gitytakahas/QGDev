#ifndef Local_QGLikelihoodCalculator_h
#define Local_QGLikelihoodCalculator_h
#include <math.h>
#include <TFile.h>
#include <TKey.h>
#include <TVector.h>
#include <vector>
#include <map>
#include "binFunctions.h"

/**
 * * A modified version of /RecoJets/JetAlgorithms/src/QGLikelihoodCalculator.cc, working with ROOT files instead of database obect
 * */

class QGLikelihoodCalculator{
  public:
    QGLikelihoodCalculator(){};
    ~QGLikelihoodCalculator();
    bool init(TString fileName);
    float computeQGLikelihood(float pt, float eta, float rho, std::vector<float> vars);

  private:
    TH1F* findEntry(float eta, float pt, float rho, int qgIndex, int varIndex);
    bool isValidRange(float pt, float rho, float eta);

    std::vector<float> etaBins, ptBinsC, ptBinsF, rhoBins;
    std::map<TString, TH1F*> pdfs; 
    TFile* f;
};



float QGLikelihoodCalculator::computeQGLikelihood(float pt, float eta, float rho, std::vector<float> var){
  if(!isValidRange(pt, rho, eta)) return -1;

  float Q=1., G=1.;
  for(unsigned int varIndex = 0; varIndex < vars.size(); ++varIndex){

    auto qgEntry = findEntry(eta, pt, rho, 0, varIndex);
    if(!qgEntry) return -1; 
    float Qi = qgEntry->GetBinContent(qgEntry->FindBin(vars[varIndex]));
    float mQ = qgEntry->GetMean();

    qgEntry = findEntry(eta, pt, rho, 1, varIndex); 
    if(!qgEntry) return -1;
    float Gi = qgEntry->GetBinContent(qgEntry->FindBin(vars[varIndex]));
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


bool QGLikelihoodCalculator::init(TString fileName){
  f = new TFile("../data/" + fileName);
  if(f->IsZombie()) 				return false;
  if(!getBinsFromFile(etaBins, "etaBins", f))	return false;
  if(!getBinsFromFile(ptBinsC, "ptBinsC", f))	return false;
  if(!getBinsFromFile(ptBinsF, "ptBinsF", f))	return false;
  if(!getBinsFromFile(rhoBins, "rhoBins", f))	return false;
  std::cout << "localQGLikelihoodCalculator: Initialized binning of pdfs..." << std::endl;

  TList *keys = f->GetListOfKeys();
  if(!keys) return false;

  TIter nextdir(keys);
  TKey *keydir;
  while((keydir = (TKey*)nextdir())){
    if(!keydir->IsFolder()) continue;
    TDirectory *dir = (TDirectory*)keydir->ReadObj() ;
    TIter nexthist(dir->GetListOfKeys());
    TKey *keyhist;
    while((keyhist = (TKey*)nexthist())){
      pdfs[keyhist->GetName()] = (TH1F*) keyhist->ReadObj(); 
    }
  }
  std::cout << "localQGLikelihoodCalculater: pdfs initialized..." << std::endl;

  return true;
}


QGLikelihoodCalculator::~QGLikelihoodCalculator(){
  for(auto& pdf : pdfs) delete pdf.second;
  delete f;
}
#endif
