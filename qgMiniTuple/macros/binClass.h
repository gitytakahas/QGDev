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
    std::vector<TString> 			varNames;
    std::map<TString, TString> 			varDisplayNames;
    std::map<TString, float*> 			varPointers;
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
    void 			setLink(TString, TString);
    void 			setReference(TString name, float* pointer){ 		 varPointers[name] = pointer;}
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
  if(!varNames.size()) return std::vector<TString>(1, "nobins");

  std::vector<TString> allBinNames = std::vector<TString>(1, "");
  std::vector<TString> baseNames;
  for(TString& varName : varNames){
    baseNames = allBinNames;
    allBinNames.clear();
    for(TString baseName : baseNames){
      for(int j = 0; j < getNBins(varName); ++j){
        TString name = baseName + (baseName == "" ? "" : "_") + varName + j;
        if(alsoLinked || masterOfBin(name) == "") allBinNames.push_back(name);
      }
    }
  }
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
bool binClass::getBinNumber(std::vector<float>& bins, float* var, int& bin){
  if(*var < bins.front() || *var > bins.back()) return false;
  auto binUp = bins.begin() + 1;
  while(*var > *binUp) ++binUp;
  bin = binUp - bins.begin() - 1;
  return true;
}




/*
 * Find lower and upper edge from variable in bin
 */
float binClass::getUpperEdge(TString binName, TString varName){ return getLowerEdge(binName, varName, true);}
float binClass::getLowerEdge(TString binName, TString varName, bool upper = false){
  if(!contains(varNames, varName)) return -1;

  TString sub = binName(binName.Index(TRegexp(varName+"[0-9]")) + varName.Length(), binName.Length());
  TString bin = sub.Contains('_')? sub(0, sub.Index('_')) : sub;
  int index = bin.Atoi() + (upper? 1 : 0);
  return binRanges[varName][index];
}




/*
 * Bin initialization and merging functions
 */
// Set bin ranges for variable
void binClass::setBinRange(TString varName, TString displayName, std::vector<float> varBinning){
  varNames.push_back(varName);												// The varNames vector keeps the order
  varDisplayNames[varName] = displayName;
  binRanges[varName]       = varBinning;
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
  binName = (varNames.size()? "" : "nobins");
  for(TString& varName : varNames){
    if(!varPointers.count(varName)){ std::cout << "No reference var for " << varName << " binning!" << std::endl; exit(1);}
    int binNumber;
    if(!getBinNumber(binRanges[varName], varPointers[varName], binNumber)) return false;
    binName += (binName == "" ? "" : "_") + varName + binNumber;
  }
  if(masterOfBin(binName) != "") binName = masterOfBin(binName);
  return true;
}



/*
 * Print and write functions
 */
// Print out all the bin ranges
void binClass::printBinRanges(){
  for(TString& varName : varNames){
    std::cout << "Binning for " << varName << ": {";
    for(float j : binRanges[varName]) std::cout << j << ", "; std::cout << "\b\b}" << std::endl;
  }
  std::cout << std::endl;
}


// Write bins to open ROOT file
void binClass::writeBinsToFile(){
  for(TString& varName : varNames){
    TVectorT<float> tbins(binRanges[varName].size(), binRanges[varName].data());
    tbins.Write(varName + "Bins");
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
  if(varNames.size() > 3){
    height = 0.08/(float)varNames.size();
    t.SetTextSize(0.07/(float) varNames.size());
  }
  for(int i = 0; i < varNames.size(); ++i){
    t.DrawLatex(0.9,0.91+i*height, TString::Format("%.2f < " + varDisplayNames[varNames[i]] + " < %.2f", getLowerEdge(binName, varNames[i]), getUpperEdge(binName, varNames[i])));
  }
  t.SetTextAlign(11);
  t.DrawLatex(0.1,0.91,  jetType);
}


// Construct latex loop containing \bin, \min and \max variables
void binClass::makeTexLoops(){
  for(int i = 0; i < varNames.size(); ++i){
    ofstream output;
    output.open((varNames[i] + "Bins.tex").Data());
    output << "\\foreach \\bin in {0,...," << getNBins(varNames[i]) - 1 << "}{" << std::endl;
    for(int j = 0; j < getNBins(varNames[i]); ++j){
      output << "  \\ifthenelse{\\equal{\\" + varNames[i] + "bin}{" << i << "}}";
      output << "{\\def\\" << varNames[i] << "Min{" << binRanges[varNames[i]][j] << "}";
      output << "\\def\\" << varNames[i] << "Max{" << binRanges[varNames[i]][j+1] << "}}{}" << std::endl;
    }
    output << "}" << std::endl;
    output.close();
  }
}
#endif
