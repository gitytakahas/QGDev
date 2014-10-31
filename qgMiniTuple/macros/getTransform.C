#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "binFunctions.h"
#include "treeLooper.h"

// Function to merge bins
void addTo(std::multimap<TString, std::vector<std::vector<int>>>& mymultimap, TString variableAndType, std::vector<std::vector<int>> binList){
  for(auto& i : mymultimap){
    if(i.first == variableAndType && i.second[0] == binList[0]){
      i.second.insert(i.second.end(), binList.begin() + 1, binList.end());
      return;
    }
  }
  mymultimap.insert(std::pair<TString, std::vector<std::vector<int>>>(variableAndType, binList));
}

// Print matrix
template<typename T> void printMatrix(TString name, std::vector<std::vector<T>> matrix){
    std::cout << name + " matrix:" << std::endl << std::endl;
    for(auto& row : matrix){
      for(double& element : row) std::cout << std::setw(15) << element;
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}

// Multiply nxn matrix
template<typename T> std::vector<std::vector<T>> multiply(std::vector<std::vector<T>> matrix1, std::vector<std::vector<T>> matrix2){
  std::vector<std::vector<T>> resultMatrix(matrix1);
  for(int i = 0; i < matrix1.size(); ++i){
    for(int j = 0; j < matrix1.size(); ++j){
      resultMatrix[i][j] = 0;
      for(int k = 0; k < matrix1.size(); ++k) resultMatrix[i][j] += matrix1[i][k]*matrix2[k][j];
    }
  }
  return resultMatrix;
}

// Main program
int main(int argc, char**argv){
  bool fineBinning = false;

  // Define binning for pdfs
  std::vector<float> etaBins = {0,1.3,1.5,2,2.5,3,4.7};
  std::vector<float> ptBins; getBins(ptBins, 20, 20, 2000, true); ptBins.push_back(6500);
  std::vector<float> rhoBins = {0,9999};

  // Check binning
  printBins("rho", rhoBins);
  printBins("eta", etaBins);
  printBins("pt", ptBins); std::cout << std::endl;

  // Link some bins to be merged (order eta - pt - rho) because of low statistics (for example higher pT bins at large eta)
  // Should be redefined when changes are made to eta-pt-rho grid
  std::multimap<TString, std::vector<std::vector<int>>> associatedBins;
  for(TString type : {"gluon","quark"}){
    for(TString var : {"axis2","mult","ptD"}){
      for(int i=10; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{5,9,0}, {5,i,0}});
      for(int i=14; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{4,13,0}, {4,i,0}});
      for(int i=16; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{3,15,0}, {3,i,0}});
      for(int i=18; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{2,17,0}, {2,i,0}});
      for(int i=19; i <= 21; ++i) addTo(associatedBins, var + "_" + type, {{1,18,0}, {1,i,0}});
    }
  }

  // For different jet types
  for(TString jetType : {"AK4chs"}){//,"AK5","AK5chs","AK7","AK7chs"}){
    treeLooper t("QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8_S14", jetType);
    t.setMaxEntries(1000);

    int nEvents = 0;
    std::vector<double> meanVector(3, 0);
    std::vector<std::vector<double>> covarianceMatrix(3, std::vector<double>(3, 0));
    for(bool calcMean : {true, false}){
      while(t.next()){
        int rhoBin, etaBin, ptBin;
        if(!getBinNumber(rhoBins, t.rho, rhoBin)) 	continue;
        if(!getBinNumber(etaBins, t.eta, etaBin)) 	continue;
        if(!getBinNumber(ptBins,  t.pt,  ptBin)) 	continue;

        if(t.jetIdLevel < 3) continue;											// Select tight jets
        if(!t.matchedJet) continue; 
        if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
        if((fabs(t.partonId) > 3 && t.partonId != 21)) continue;							// Keep only udsg
        if(t.bTag) continue;												// Anti-b tagging
        if(t.axis2 > 8 || t.mult > 100 || t.ptD > 1.001) continue;							// Do not keep the large values

//        std::vector<double> varVector{(1.06188*axis2), (0.0326763*axis2+0.0833892*mult), (-1.40839*axis2+0.0688766*mult+13.5899*ptD)};
        std::vector<double> varVector{(double) t.axis2, (double) t.mult, (double) t.ptD};

        if(calcMean){													// First iteration: calculate the mean
          ++nEvents;
          for(int i = 0; i < varVector.size(); ++i) meanVector[i] += varVector[i]; 
        } else {													// Second itereation: calculate the covariance matrix
          for(int i = 0; i < varVector.size(); ++i){
            for(int j = 0; j < varVector.size(); ++j){ 
              covarianceMatrix[i][j] += (varVector[i] - meanVector[i])*(varVector[j] - meanVector[j]);
            }
          }
        }
      }
      if(calcMean){
        for(double& mean : meanVector) mean /= (double) nEvents;
      } else {
        for(auto& row : covarianceMatrix){
          for(double& element : row) element /= (double) nEvents;
        }
        printMatrix("Covariance", covarianceMatrix);
      }
    }

    // Calculate the correlation matrix
    std::vector<std::vector<double>> covToCorrMatrix(covarianceMatrix);
    for(int i = 0; i < covToCorrMatrix.size(); ++i){
      for(int j = 0; j < covToCorrMatrix.size(); ++j){
        if(i != j) covToCorrMatrix[i][j] = 0;
        else covToCorrMatrix[i][i] = 1./std::sqrt(covarianceMatrix[i][i]);
      }
    }
    printMatrix("Correlation", multiply(covToCorrMatrix, multiply(covarianceMatrix, covToCorrMatrix)));

    // Use Cholesky algorithm to get square-root matrix
    std::vector<std::vector<double>> choleskyMatrix(covarianceMatrix);								// See wikipedia
    for(int j = 0; j < choleskyMatrix.size(); ++j){
      for(int i = 0; i < j; ++i) choleskyMatrix[i][j] = 0;
      for(int k = 0; k < j; ++k) choleskyMatrix[j][j] -= choleskyMatrix[j][k]*choleskyMatrix[j][k];
      choleskyMatrix[j][j] = std::sqrt(choleskyMatrix[j][j]);
      for(int i = j + 1; i < choleskyMatrix.size(); ++i){
        for(int k = 0; k < j; ++k) choleskyMatrix[i][j] -= choleskyMatrix[i][k]*choleskyMatrix[j][k];
        choleskyMatrix[i][j] /= choleskyMatrix[j][j];
      }
    }
    printMatrix("Cholesky", choleskyMatrix);

    // Take inverse of Cholesky matrix to get transformation matrix
    double detCholesky = 1;
    for(int i = 0; i < choleskyMatrix.size(); ++i) detCholesky *= choleskyMatrix[i][i];
    std::vector<std::vector<double>> choleskyAdjugateMatrix(choleskyMatrix);
    choleskyAdjugateMatrix[0][0] = choleskyMatrix[1][1]*choleskyMatrix[2][2];							// Only holds for 3x3
    choleskyAdjugateMatrix[1][1] = choleskyMatrix[2][2]*choleskyMatrix[0][0];
    choleskyAdjugateMatrix[2][2] = choleskyMatrix[0][0]*choleskyMatrix[1][1];
    choleskyAdjugateMatrix[1][0] = -choleskyMatrix[1][0]*choleskyMatrix[2][2];
    choleskyAdjugateMatrix[2][1] = -choleskyMatrix[0][0]*choleskyMatrix[2][1];
    choleskyAdjugateMatrix[2][0] = choleskyMatrix[1][0]*choleskyMatrix[2][1]-choleskyMatrix[1][1]*choleskyMatrix[2][0];
    std::vector<std::vector<double>> transformMatrix(choleskyAdjugateMatrix);
    for(auto& row : transformMatrix){
      for(double& element : row) element /= detCholesky;
    }
    printMatrix("Transform", transformMatrix);
  }

  return 0;
}
