#include <vector>
#include <map>
#include <stdexcept>
#include <cmath>

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
    if(vars[varIndex] < -0.5) continue; //use to inspect variables separately (i.e. skip if feeding -1)

    auto quarkEntry = findEntry(eta, pt, rho, 0, varIndex);
    auto gluonEntry = findEntry(eta, pt, rho, 1, varIndex);
    if(!quarkEntry && !gluonEntry) return -2;

    int binQ = quarkEntry->FindBin(vars[varIndex]);
    float Qi = quarkEntry->GetBinContent(binQ);
    float Qw = quarkEntry->GetBinWidth(binQ);

    int binG = gluonEntry->FindBin(vars[varIndex]);
    float Gi = gluonEntry->GetBinContent(binG);
    float Gw = gluonEntry->GetBinWidth(binG);

    if(Qi <= 0 || Gi <= 0){
      bool gluonAboveQuark = (quarkEntry->GetMean() < gluonEntry->GetMean());

      if(quarkEntry->Integral(0, binQ) + gluonEntry->Integral(0, binG) <=0){
        if(gluonAboveQuark){ Qi = 0.999; Gi = 0.001;} else { Qi = 0.001; Gi = 0.999;}
      } else if(quarkEntry->Integral(binQ, quarkEntry->GetNbinsX() + 1) + gluonEntry->Integral(binG, gluonEntry->GetNbinsX() + 1) <=0){
        if(gluonAboveQuark){ Qi = 0.001; Gi = 0.999;} else { Qi = 0.999; Gi = 0.001;}
      } else {
        int q = 1, g = 1;
        while(Qi <= 0){
          Qi += quarkEntry->GetBinContent(binQ+q) + quarkEntry->GetBinContent(binQ-q);
          Qw += quarkEntry->GetBinWidth(binQ+q) + quarkEntry->GetBinWidth(binQ-q);
          ++q;
        }
        while(Gi <= 0){
          Gi += gluonEntry->GetBinContent(binG+g) + gluonEntry->GetBinContent(binG-g);
          Gw += gluonEntry->GetBinWidth(binG+g) + gluonEntry->GetBinWidth(binG-g);
          ++g;
        }
      }
    } 
    Q*=Qi/Qw;
    G*=Gi/Gw;
  }

  if(Q==0) return 0;
  return Q/(Q+G);
}


/// Compute the QGLikelihood using CDF, given the pT, eta, rho and likelihood variables vector
float QGLikelihoodCalculator::computeQGLikelihoodCDF(float pt, float eta, float rho, std::vector<float> vars){
  if(!isValidRange(pt, rho, eta)) return -1;

  float Q=1., G=1.;
  for(unsigned int varIndex = 0; varIndex < vars.size(); ++varIndex){
    if(vars[varIndex] < -0.5) continue; //use to inspect variables separately (i.e. skip if feeding -1)

    auto quarkEntry = findEntry(eta, pt, rho, 0, varIndex);
    auto gluonEntry = findEntry(eta, pt, rho, 1, varIndex);
    if(!quarkEntry || !gluonEntry) return -2;

    float cdfQuark = quarkEntry->Integral(0, quarkEntry->FindBin(vars[varIndex]));
    float cdfGluon = gluonEntry->Integral(0, gluonEntry->FindBin(vars[varIndex]));

    float ccdfQuark = 1.-cdfQuark;
    float ccdfGluon = 1.-cdfGluon;

    float Qi, Gi;
    if(quarkEntry->GetMean() < gluonEntry->GetMean()){
      if(cdfQuark+cdfGluon <= 0){
        Qi = 0.999;
        Gi = 0.001;
      } else if(ccdfQuark+ccdfGluon <= 0){
        Qi = 0.001;
        Gi = 0.999;
      } else {
        Qi = cdfQuark/(cdfQuark+cdfGluon) - 0.5;
        Gi = ccdfGluon/(ccdfQuark+ccdfGluon) - 0.5;
      }
    } else {
      if(cdfQuark+cdfGluon <= 0){
        Qi = 0.001;
        Gi = 0.999;
      } else if(ccdfQuark+ccdfGluon <= 0){
        Qi = 0.999;
        Gi = 0.001;
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
