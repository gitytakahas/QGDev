#ifndef binClass_h
#define binClass_h
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TVector.h"
#include "TLatex.h"
#include "TRegexp.h"

// Easy way to add integers to TString (useful to construct bin names) 
TString operator+(TString left, int right){ TString result = left; result += right; return result;}

// Fast way to check if vector contains element
template <typename T> const bool contains(std::vector<T>& vec, const T& element){
  return (std::find(vec.begin(), vec.end(), element) != vec.end());
}




/*
 * Class to take care of binning the pdfs/plots and related functions
 */
class binClass{
  private:
    std::vector<TString> 			valueNames;
    std::map<TString, TString> 			valueDisplayNames;
    std::map<TString, float*> 			valuePointers;
    std::map<TString, std::vector<float>> 	binRanges;
    std::map<TString, std::vector<TString>> 	linkedBins;
    std::map<TString, TString> 			masterBins;
    std::map<TString, std::vector<float>> 	varWeights;
    
  public:
    binClass(){};
    ~binClass(){};

    // Functions to define bins, store pointer to the variable related to the bins, merge bins, and set likelihood weights
    std::vector<float> 		getBins(int, float, float, bool, std::vector<float>);
    void 			setBinRange(TString, TString, std::vector<float>);
    void 			setReference(TString, float*);
    void 			setLink(TString, TString);
    void 			setWeights(TString binName, std::vector<float> weights){ varWeights[binName] = weights;}

    // Functions to retrieve all bin names, and the links between them 
    std::vector<TString> 	getAllBinNames(bool);
    std::vector<TString> 	getLinkedBins(TString bin){ return (linkedBins.count(bin) ? linkedBins[bin] : std::vector<TString>());}
    TString 			masterOfBin(TString bin){   return (masterBins.count(bin) ? masterBins[bin] : "");}

    // Functions to retrieve bin name, bin numbers,...
    bool 			getBinName(TString&);
    bool 			getBinNumber(std::vector<float>&, float*, int&);
    float 			getLowerEdge(TString, TString, bool);
    float 			getUpperEdge(TString, TString);
    int				getNBins(TString varName){ return binRanges[varName].size() - 1;};

    // Functions to print or write binning information
    void 			printBinRanges(); 
    void 			printInfoOnPlot(TString, TString);
    void 			writeBinsToFile();
    void 			writeWeightsToFile(TFile *file);
    void 			makeTexLoops();
};




/*
 * Returns the names of all bins
 */
std::vector<TString> binClass::getAllBinNames(bool alsoLinked = false){
  std::vector<TString> allBinNames;
  std::vector<TString> baseNames;
  if(valueNames.size() > 0){
    for(int j = 0; j < binRanges[valueNames[0]].size() - 1; ++j) allBinNames.push_back(valueNames[0] + j);
    for(int i = 1; i < valueNames.size(); ++i){
      baseNames = allBinNames;
      allBinNames.clear();
      for(TString baseName : baseNames){
        for(int j = 0; j < binRanges[valueNames[i]].size() - 1; ++j){
          TString name = baseName + "_" + valueNames[i] + j;
          if(alsoLinked || masterOfBin(name) == "") allBinNames.push_back(name);
        }
      }
    }
  }
  if(!allBinNames.size()) allBinNames.push_back("nobins");
  return allBinNames;
}




/* 
 * Constructs vector with nBins between min and max, with additional bins added if needed
 */
std::vector<float> binClass::getBins(int nBins, float min, float max, bool log, std::vector<float> additionalBins = {}){
  std::vector<float> bins;
  const float dx = (log ? std::pow((max/min), (1./(float)nBins)) : ((max - min)/(float)nBins));
  bins.push_back(min);
  while(bins.back() < max) bins.push_back(log ? bins.back()*dx : bins.back()+dx);
  bins.insert(bins.end(), additionalBins.begin(), additionalBins.end());
  return bins;
}




/* 
 * Gets the bin number for a variable
 */
bool binClass::getBinNumber(std::vector<float>& bins, float* value, int& bin){
  if(*value < bins.front() || *value > bins.back()) return false;
  auto binUp = bins.begin() + 1;
  while(*value > *binUp) ++binUp;
  bin = binUp - bins.begin() - 1;
  return true;
}




/*
 * Find lower and upper edge from variable in bin
 */
float binClass::getUpperEdge(TString binName, TString varName){ return getLowerEdge(binName, varName, true);}
float binClass::getLowerEdge(TString binName, TString varName, bool upper = false){
  TString sub = binName(binName.Index(TRegexp(varName+"[0-9]")) + varName.Length(), binName.Length());
  TString bin = sub.Contains('_')? sub(0, sub.Index('_')) : sub;
  int index = bin.Atoi() + (upper? 1 : 0);
  auto it = std::find(valueNames.begin(), valueNames.end(), varName);
  if(it == valueNames.end()) return -1;
  else return binRanges[*it][index];
}




/*
 * Bin initialization and merging functions
 */
// Set bin ranges for variable
void binClass::setBinRange(TString name, TString displayName, std::vector<float> varBinning){
  valueNames.push_back(name);
  valueDisplayNames[name] = displayName;
  binRanges[name] = varBinning;
}


// Connects variable from tree
void binClass::setReference(TString name, float* pointer){
  valuePointers[name] = pointer;
}


// Redirects one bin to another, taking into account already existing links
void binClass::setLink(TString master, TString toBeLinked){
  auto allBinNames = getAllBinNames(true);
  if(!contains(allBinNames, master)){     std::cout << master     << " does not exist!" << std::endl; exit(1);}		// Check if master and toBeLinked exist and are not the sane
  if(!contains(allBinNames, toBeLinked)){ std::cout << toBeLinked << " does not exist!" << std::endl; return;}
  if(master == toBeLinked) return;

  if(masterOfBin(master) != "") master = masterOfBin(master);								// If master was linked earlier to another bin, use that bin as master
  masterBins[toBeLinked] = master;

  if(linkedBins.find(master) == linkedBins.end()) linkedBins[master] = std::vector<TString>();				// Create new list of linked bins if not yet existing
  std::cout << "Merging " << toBeLinked << " into " << master << std::endl;
  linkedBins[master].push_back(toBeLinked);

  if(linkedBins.count(toBeLinked)){											// Redirect bins which are linked to toBeLinked
    for(auto j : linkedBins[toBeLinked]){
      linkedBins[master].push_back(j);
      masterBins[j] = master;
    }
    linkedBins.erase(toBeLinked);
  }
}


/*
 * Find the bin based on the variable pointers and bins, returns true if found
 * Redirects to master bin if we are dealing with a linked bin
 */
bool binClass::getBinName(TString& binName){
  binName = (valueNames.size()? "" : "nobins");
  for(int i = 0; i < valueNames.size(); ++i){
    if(!valuePointers.count(valueNames[i])){ std::cout << "No reference value for " << valueNames[i] << " binning!" << std::endl; exit(1);}
    int binNumber;
    if(!getBinNumber(binRanges[valueNames[i]], valuePointers[valueNames[i]], binNumber)) return false;
    binName += (i == 0 ? "" : "_") + valueNames[i] + binNumber;
  }
  if(masterOfBin(binName) != "") binName = masterOfBin(binName);
  return true;
}



/*
 * Print and write functions
 */
// Print out all the bin ranges
void binClass::printBinRanges(){
  for(int i = 0; i < valueNames.size(); ++i){
    std::cout << "Binning for " << valueNames[i] << ": {";
    for(float j : binRanges[valueNames[i]]) std::cout << j << ", "; std::cout << "\b\b}" << std::endl;
  }
  std::cout << std::endl;
}


// Write bins to open ROOT file
void binClass::writeBinsToFile(){
  for(int i = 0; i < valueNames.size(); ++i){
    TVectorT<float> tbins(binRanges[valueNames[i]].size(), binRanges[valueNames[i]].data());
    tbins.Write(valueNames[i] + "Bins");
  }
}


// Write weights to open ROOT file
void binClass::writeWeightsToFile(TFile *file){
  file->mkdir("weights");
  file->cd("weights");
  for(TString bin : getAllBinNames(true)){
    if(!varWeights.count(bin)) varWeights[bin] = {1.,1.,1.};
    TVectorT<float> tweights(varWeights[bin].size(), varWeights[bin].data());
    tweights.Write(bin);
  }
  file->cd();
}


// Print bin ranges on plot
void binClass::printInfoOnPlot(TString binName, TString jetType){
  TLatex t;
  t.SetNDC(kTRUE);
  t.SetTextAlign(31);
  t.SetTextSize(0.02);
  float height = 0.025;
  if(valueNames.size() > 3){
    height = 0.08/(float)valueNames.size();
    t.SetTextSize(0.07/(float) valueNames.size());
  }
  for(int i = 0; i < valueNames.size(); ++i){
    t.DrawLatex(0.9,0.91+i*height, TString::Format("%.2f < " + valueDisplayNames[valueNames[i]] + " < %.2f", getLowerEdge(binName, valueNames[i]), getUpperEdge(binName, valueNames[i])));
  }
  t.SetTextAlign(11);
  t.DrawLatex(0.1,0.91,  jetType);
}


// Construct latex loop containing \bin, \min and \max variables
void binClass::makeTexLoops(){
  for(int i = 0; i < valueNames.size(); ++i){
    ofstream output;
    output.open((valueNames[i] + "Bins.tex").Data());
    output << "\\foreach \\bin in {0,...," << binRanges[valueNames[i]].size() - 2 << "}{" << std::endl;
    for(int j = 0; j < binRanges[valueNames[i]].size() - 1; ++j){
      output << "  \\ifthenelse{\\equal{\\" + valueNames[i] + "bin}{" << i << "}}";
      output << "{\\def\\" << valueNames[i] << "Min{" << binRanges[valueNames[i]][j] << "}";
      output << "\\def\\" << valueNames[i] << "Max{" << binRanges[valueNames[i]][j+1] << "}}{}" << std::endl;
    }
    output << "}" << std::endl;
    output.close();
  }
}
#endif
