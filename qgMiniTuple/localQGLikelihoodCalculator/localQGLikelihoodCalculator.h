#ifndef Local_QGLikelihoodCalculator_h
#define Local_QGLikelihoodCalculator_h
#include <vector>
#include <map>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

/*
 * A modified version of /RecoJets/JetAlgorithms/interface/QGLikelihoodCalculator.h, working with ROOT files instead of database obect
 */

class QGLikelihoodCalculator{
  public:
    QGLikelihoodCalculator(const TString& fileName);
    ~QGLikelihoodCalculator();
    float computeQGLikelihood(float pt, float eta, float rho, std::vector<float> vars);
    float computeQGLikelihoodCDF(float pt, float eta, float rho, std::vector<float> vars);

  private:
    bool init(const TString& fileName);
    TH1F* findEntry(float eta, float pt, float rho, int qgIndex, int varIndex);
    bool isValidRange(float pt, float rho, float eta);
    bool getBinsFromFile(std::vector<float>& bins, const TString& name, TFile* f);
    bool getBinNumber(std::vector<float>& bins, float value, int& bin);

    std::vector<float> etaBins, ptBinsC, ptBinsF, rhoBins;
    std::map<TString, TH1F*> pdfs; 
    TFile* f;
};
#endif
