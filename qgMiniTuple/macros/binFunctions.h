#ifndef binFunctions_h
#define binFunctions_h
#include <vector>
#include <iostream>
#include "TVector.h"

void getBins(std::vector<float>& bins, int nBins, double min, double max, bool log){
  const float dx = (log ? std::pow((max/min), (1./(float)nBins)) : ((max - min)/(float)nBins));
  bins.push_back(min);
  while(bins.back() < max) bins.push_back(log ? bins.back()*dx : bins.back()+dx);
}


void printBins(TString name, std::vector<float>& bins){
  std::cout << "Binning for " << name << ": {";
  for(float i : bins) std::cout << i << ", "; std::cout << "\b\b}" << std::endl;
}


bool getBinNumber(std::vector<float>& bins, float value, int& bin){
  if(value < bins.front() || value > bins.back()) return false;
  auto binUp = bins.begin() + 1;
  while(value > *binUp) ++binUp;
  bin = binUp - bins.begin() - 1;
  return true;
}


int getBinFromString(TString s, TString var){
  TString sub = s(s.Index(var + '-') + var.Length() + 1, s.Length());
  TString bin = sub.Contains('_')? sub(0, sub.Index('_')) : sub;
  return bin.Atoi();
}


int getNBins(std::vector<float>& bins){ return bins.size() - 1;}


void writeBinsToFile(std::vector<float>& bins, TString name){
  TVectorT<float> tbins(bins.size(), &bins[0]);
  tbins.Write(name);
}


bool getBinsFromFile(std::vector<float>& bins, TString name, TFile* f){
  TVectorT<float> *tbins = nullptr;
  f->GetObject(name, tbins);
  if(!tbins) return false;
  for(int i = 0; i < tbins->GetNoElements(); ++i) bins.push_back((*tbins)[i]);
  return true;
}
#endif
