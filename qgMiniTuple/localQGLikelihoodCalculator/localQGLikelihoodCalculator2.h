#ifndef Local_QGLikelihoodCalculator2_h
#define Local_QGLikelihoodCalculator2_h
#include <vector>
#include <map>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

/*
 * A modified version of /RecoJets/JetAlgorithms/interface/QGLikelihoodCalculator.h, working with ROOT files instead of database obect
 */

class QGLikelihoodCalculator2{
  public:
    QGLikelihoodCalculator2(const TString& fileName);
    ~QGLikelihoodCalculator2();
    float computeQGLikelihood(float pt, float eta, std::vector<float> vars);
    bool getBinsFromFile(std::vector<float>& bins, const TString& name );

  private:
    bool init(const TString& fileName);
    TH1F* findEntry(float eta, float pt, int qgIndex, int varIndex);
    bool isValidRange(float pt, float eta);
    bool getBinNumber(std::vector<float>& bins, float value, int& bin);

    std::vector<float> etaBins, ptBinsC;
    std::map<TString, TH1F*> pdfs; 
    TFile* f;
};
#endif
