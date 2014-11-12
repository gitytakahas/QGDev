#ifndef binClass_h
#define binClass_h
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TVector.h"

/*
 * Class to take care of binning the pdfs/plots and related functions
 */
class binClass{
  private:
    std::vector<TString> valueNames;
    std::vector<float*> valuePointers;
    std::vector<std::vector<float>> binRanges;
    std::map<TString, std::vector<TString>> linkedBins;
    std::map<TString, TString> masterBins;

  public:
    binClass(){};
    ~binClass(){};
    bool update();
    bool getBinNumber(std::vector<float>&, float*, int&);
    void setBinRange(TString, std::vector<float>);
    void printBinRanges(); 
    void setReference(TString, float*);
    std::vector<float> getBins(int, float, float, bool, std::vector<float>);
    std::vector<TString> getAllBinNames(bool);
    std::vector<TString> getLinkedBins(TString bin){ return (linkedBins.find(bin) == linkedBins.end() ? std::vector<TString>() : linkedBins[bin]);}
    TString masterOfBin(TString bin){ return (masterBins.find(bin) == masterBins.end() ? "" : masterBins[bin]);}
    void setLinks(TString, std::vector<TString>);
    int getBinFromString(TString, TString);
    float getLowerEdge(TString, TString, bool);
    float getUpperEdge(TString, TString);
    void writeBinsToFile();
    void makeTexLoops();

    TString name;
};


// Redirects one bin to another
void binClass::setLinks(TString master, std::vector<TString> toBeLinked){
  auto allBinNames = getAllBinNames(true);
  if(std::find(allBinNames.begin(), allBinNames.end(), master) == allBinNames.end()){ std::cout << master << " does not exist!" << std::endl; exit(1);}
  if(masterOfBin(master) != "") master = masterOfBin(master);
  for(auto i : toBeLinked){
    if(std::find(allBinNames.begin(), allBinNames.end(), i) == allBinNames.end()){ std::cout << i << " does not exist!" << std::endl; continue;}
    if(master == i) continue;
    masterBins[i] = master;
    if(linkedBins.find(master) == linkedBins.end()) linkedBins[master] = std::vector<TString>();
    if(std::find(linkedBins[master].begin(), linkedBins[master].end(), i) == linkedBins[master].end()){
      std::cout << "Merging " << i << " into " << master << std::endl;
      linkedBins[master].push_back(i);
    }
  }
}


// Set bin ranges for variable
void binClass::setBinRange(TString name, std::vector<float> varBinning){
  valueNames.push_back(name);
  valuePointers.push_back(nullptr);
  binRanges.push_back(varBinning);
}


// Connects variable from tree
void binClass::setReference(TString name, float* pointer){
  auto it = std::find(valueNames.begin(), valueNames.end(), name);
  if(it != valueNames.end()) valuePointers[it - valueNames.begin()] = pointer;
}


// Print out all the bin ranges
void binClass::printBinRanges(){
  for(int i = 0; i < valueNames.size(); ++i){
    std::cout << "Binning for " << valueNames[i] << ": {";
    for(float j : binRanges[i]) std::cout << j << ", "; std::cout << "\b\b}" << std::endl;
  }
  std::cout << std::endl;
}


// To be called for every new event: true if bin is found
bool binClass::update(){
  name = "";
  if(!valuePointers.size()) name = "nobins";
  for(int i = 0; i < valuePointers.size(); ++i){
    if(!valuePointers[i]){ std::cout << "No reference value for " << valueNames[i] << " binning!" << std::endl; exit(1);}
    int binNumber;
    if(!getBinNumber(binRanges[i], valuePointers[i], binNumber)) return false;
    name += (i == 0 ? "" : "_") + valueNames[i] + TString::Format("%d",binNumber);
  }
  if(masterOfBin(name) != "") name = masterOfBin(name);
  return true;
}


// Returns the names of all bins
std::vector<TString> binClass::getAllBinNames(bool alsoLinked = false){
  std::vector<TString> allBinNames;
  std::vector<TString> baseNames;
  if(binRanges.size() > 0){
    for(int j = 0; j < binRanges[0].size() - 1; ++j) allBinNames.push_back(valueNames[0] + TString::Format("%d", j));
    for(int i = 1; i < valueNames.size(); ++i){
      baseNames = allBinNames;
      allBinNames.clear();
      for(TString baseName : baseNames){
        for(int j = 0; j < binRanges[i].size() - 1; ++j){
          TString name = baseName + "_" + valueNames[i] + TString::Format("%d",j);
          if(alsoLinked || masterOfBin(name) == "") allBinNames.push_back(name);
        }
      }
    }
  }
  if(!allBinNames.size()) allBinNames.push_back("nobins");
  return allBinNames;
}


// Constructs vector with bin ranges
std::vector<float> binClass::getBins(int nBins, float min, float max, bool log, std::vector<float> additionalBins = {}){
  std::vector<float> bins;
  const float dx = (log ? std::pow((max/min), (1./(float)nBins)) : ((max - min)/(float)nBins));
  bins.push_back(min);
  while(bins.back() < max) bins.push_back(log ? bins.back()*dx : bins.back()+dx);
  bins.insert(bins.end(), additionalBins.begin(), additionalBins.end());
  return bins;
}


// Gets the bin number for a variable
bool binClass::getBinNumber(std::vector<float>& bins, float* value, int& bin){
  if(*value < bins.front() || *value > bins.back()) return false;
  auto binUp = bins.begin() + 1;
  while(*value > *binUp) ++binUp;
  bin = binUp - bins.begin() - 1;
  return true;
}


// Find lower and upper edge from variable in bin
float binClass::getUpperEdge(TString binName, TString varName){ return getLowerEdge(binName, varName, true);}
float binClass::getLowerEdge(TString binName, TString varName, bool upper = false){
  TString sub = binName(binName.Index(varName) + varName.Length(), binName.Length());
  TString bin = sub.Contains('_')? sub(0, sub.Index('_')) : sub;
  int index = bin.Atoi() + (upper? 1 : 0);
  auto it = std::find(valueNames.begin(), valueNames.end(), varName);
  if(it == valueNames.end()) return -1;
  else return binRanges[it - valueNames.begin()][index];
}


// Write bins to open ROOT file
void binClass::writeBinsToFile(){
  for(int i = 0; i < valueNames.size(); ++i){
    TVectorT<float> tbins(binRanges[i].size(), &binRanges[i][0]);
    tbins.Write(valueNames[i] + "Bins");
  }
}


// Construct latex loop containing \bin, \min and \max variables
void binClass::makeTexLoops(){
  for(int i = 0; i < valueNames.size(); ++i){
    ofstream output;
    output.open((valueNames[i] + "Bins.tex").Data());
    output << "\\foreach \\bin in {0,...," << binRanges[i].size() - 2 << "}{" << std::endl;
    for(int j = 0; j < binRanges[i].size() - 1; ++j){
      output << "  \\ifthenelse{\\equal{\\" + valueNames[i] + "bin}{" << i << "}}";
      output << "{\\def\\" << valueNames[i] << "Min{" << binRanges[i][j] << "}";
      output << "\\def\\" << valueNames[i] << "Max{" << binRanges[i][j+1] << "}}{}" << std::endl;
    }
    output << "}" << std::endl;
    output.close();
  }
}
#endif
