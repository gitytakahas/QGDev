#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <TMatrixD.h>
#include <TFile.h>
#include <TKey.h>
#include "binClass.h"
#include "treeLooper.h"

/*
 * Old file with some helper functions to work with decorrelation matrices
 * As the likelihood with weights seemed to be a better solution than a likelihood with decorrelated variables, the code was not used anymore recently
 */

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


// Code to calculate decorrelation matrices
std::map<TString, std::vector<std::vector<double>>> getTransform(treeLooper& t, binClass& bins){

  std::map<TString, int> nEvents;
  std::map<TString, std::vector<double>> meanVectors;
  std::map<TString, std::vector<std::vector<double>>> covarianceMatrices;

  for(TString i : bins.getAllBinNames()){
    nEvents[i]			= 0;
    meanVectors[i] 		= std::vector<double>(3, 0);
    covarianceMatrices[i]	= std::vector<std::vector<double>>(3, std::vector<double>(3, 0));
  }

  // Calculate mean (first loop) and covariance matrix (second loop)
  for(bool calcMean : {true, false}){
    while(t.next()){
      if(!bins.update()) 	continue;										// Find bin and return false if outside range
      if(t.jetIdLevel < 3) 	continue;										// Select tight jets
      if(!t.matchedJet) 	continue; 										// Only matched jets
      if(t.nGenJetsInCone != 1 || t.nJetsForGenParticle != 1 || t.nGenJetsForGenParticle != 1) continue;		// Use only jets matched to exactly one gen jet and gen particle, and no other jet candidates
      if((fabs(t.partonId) > 3 && t.partonId != 21)) continue;								// Keep only udsg
      if(t.bTag) continue;												// Anti-b tagging
      if(t.axis2 > 8 || t.mult > 140 || t.ptD > 1.0001) continue;							// Do not keep the large values
      if(t.mult < 3)		continue;

      std::vector<double> varVector{(double) t.axis2, (double) t.ptD, (double) t.mult};

      if(calcMean){													// First iteration: calculate the mean
        ++nEvents[bins.name];
        for(int i = 0; i < varVector.size(); ++i) meanVectors[bins.name][i] += varVector[i]; 
      } else {														// Second itereation: calculate the covariance matrix
        for(int i = 0; i < varVector.size(); ++i){
          for(int j = 0; j < varVector.size(); ++j){ 
            covarianceMatrices[bins.name][i][j] += (varVector[i] - meanVectors[bins.name][i])*(varVector[j] - meanVectors[bins.name][j]);
          }
        }
      }
    }
    if(calcMean){
      for(auto& meanVector : meanVectors){
        for(double& mean : meanVector.second) mean /= (double) nEvents[meanVector.first];
      }
    } else {
      for(auto& covarianceMatrix : covarianceMatrices){
        for(auto& row : covarianceMatrix.second){
          for(double& element : row) element /= (double) nEvents[covarianceMatrix.first];
        }
      }
    }
  }

  std::map<TString, std::vector<std::vector<double>>> transformMatrices;
  for(auto& covarianceMatrix : covarianceMatrices){
    std::cout << "Calculating correlation and transformation matrix for " << covarianceMatrix.first << std::endl;

    // Calculate the correlation matrix
    std::vector<std::vector<double>> covToCorrMatrix(covarianceMatrix.second);
    for(int i = 0; i < covToCorrMatrix.size(); ++i){
      for(int j = 0; j < covToCorrMatrix.size(); ++j){
        if(i != j) covToCorrMatrix[i][j] = 0;
        else covToCorrMatrix[i][i] = 1./std::sqrt(covarianceMatrix.second[i][i]);
      }
    }
    printMatrix("Correlation", multiply(covToCorrMatrix, multiply(covarianceMatrix.second, covToCorrMatrix)));

    // Use Cholesky algorithm to get square-root matrix (see wikipedia for algorithm)
    std::vector<std::vector<double>> choleskyMatrix(covarianceMatrix.second);
    for(int j = 0; j < choleskyMatrix.size(); ++j){
      for(int i = 0; i < j; ++i) choleskyMatrix[i][j] = 0;
      for(int k = 0; k < j; ++k) choleskyMatrix[j][j] -= choleskyMatrix[j][k]*choleskyMatrix[j][k];
      choleskyMatrix[j][j] = std::sqrt(choleskyMatrix[j][j]);
      for(int i = j + 1; i < choleskyMatrix.size(); ++i){
        for(int k = 0; k < j; ++k) choleskyMatrix[i][j] -= choleskyMatrix[i][k]*choleskyMatrix[j][k];
        choleskyMatrix[i][j] /= choleskyMatrix[j][j];
      }
    }
    // printMatrix("Cholesky", choleskyMatrix);

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
    transformMatrices[covarianceMatrix.first] = transformMatrix;
  }

  return transformMatrices;
}


// Function to transform variables to uncorrelated variables
std::vector<double> decorrelate(std::vector<std::vector<double>>& transformMatrix, std::vector<double>& varVector){
  std::vector<double> uncorrelatedVarVector(varVector);
  for(int i = 0; i < varVector.size(); ++i){
    uncorrelatedVarVector[i] = 0; 
    for(int j = 0; j <= i; ++j) uncorrelatedVarVector[i] += transformMatrix[i][j]*varVector[j];
  }
  return uncorrelatedVarVector;
}


// Calculates the range of the uncorrelated variables
std::vector<std::vector<double>> calcRangeTransformation(std::vector<std::vector<double>>& transformMatrix, std::vector<double> varDown, std::vector<double> varUp){
  std::vector<std::vector<double>> ranges(2, std::vector<double>(varDown));
  for(int i = 0; i < varDown.size(); ++i){
    for(int k : {0,1}) ranges[k][i] = 0;
    for(int j = 0; j <= i; ++j){
      if(transformMatrix[i][j] > 0){
        ranges[0][i] += transformMatrix[i][j]*varDown[j];
        ranges[1][i] += transformMatrix[i][j]*varUp[j];
      } else {
        ranges[0][i] += transformMatrix[i][j]*varUp[j];
        ranges[1][i] += transformMatrix[i][j]*varDown[j];
      }
    }
  }
  return ranges;
}


// Writes matrices to file
void writeMatricesToFile(std::map<TString, std::vector<std::vector<double>>>& matrices, binClass& bins){
  for(auto& matrix : matrices){
    TMatrixD tmatrix(matrix.second.size(), matrix.second.size());
    for(int i = 0; i < matrix.second.size(); ++i){
      for(int j = 0; j < matrix.second.size(); ++j){
        tmatrix[i][j] = matrix.second[i][j];
      }
    }
    tmatrix.Write(matrix.first);
    for(auto i : bins.getLinkedBins(matrix.first)) tmatrix.Write(i);
  }
}


// Reads the matrices from file
bool getMatricesFromFile(std::map<TString, std::vector<std::vector<double>>>& matrices, TString fileName){
  TFile *f = new TFile(fileName);
  if(f->IsZombie()) return false;

  TList *keys = f->GetListOfKeys();
  if(!keys) return false;

  TIter nextdir(keys);
  TKey *keydir;
  while((keydir = (TKey*) nextdir())){
    if(!keydir->IsFolder()) continue;
    TDirectory *dir = (TDirectory*) keydir->ReadObj();
    if(TString(dir->GetName()) == "decorrelationMatrices"){
      TIter nextmatrix(dir->GetListOfKeys());
      TKey *keymatrix;
      while((keymatrix = (TKey*) nextmatrix())){
        TMatrixD* tmatrix = (TMatrixD*) keymatrix->ReadObj();
        matrices[keymatrix->GetName()] = std::vector<std::vector<double>>(3, std::vector<double>(3, 0));
        for(int i = 0; i < 3; ++i){
          for(int j = 0; j < 3; ++j){
            matrices[keymatrix->GetName()][i][j] = (*tmatrix)[i][j];
          }
        }
      }
    }
  }
  delete f;
  return true;
}
