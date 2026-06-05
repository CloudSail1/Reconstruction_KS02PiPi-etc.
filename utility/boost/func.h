#include <Math/GenVector/Boost.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/Vector4Dfwd.h>
#include <vector>
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector> LorentzVector;
using ROOT::Math::Boost;

// in func/boost.cxx
struct GammaData {
    double ang_gamma_pip;
    double ang_gamma_pin;
    double e_gamma;
    double d3_gamma;
    LorentzVector gamma;
    LorentzVector pip;
    LorentzVector pin;
};
GammaData BoostGammaInto(double* gamma, double* pip, double* pim, std::vector<double*> frame_fractions, bool is_inverse=false);
GammaData BoostGammaIntoPiPiFrame(double* gamma, double* pip, double* pim);
GammaData BoostGammaIntoPiPiGammaFrame(double* gamma, double* pip, double* pim);
GammaData BoostGammaIntoLabFrame(double* gamma, double* pip, double* pim, double* total);

// in func/dang.cxx
double MinDAng(int size, double* dang);

// in func/recoil.cxx
struct RecoilData {
    double m;
    double e;
    double d3;
    LorentzVector p4;
};
RecoilData CalRecoil(std::vector<double*> terms, double* total);
RecoilData CalRecoilPip(double *pip, double *pim, double *gamma, double *pi, double *k, double *total);
RecoilData CalRecoilPipPim(double *pip, double *pim, double *gamma, double *pi, double *k, double *total);
RecoilData CalRecoilPi(double *pip, double *pim, double *gamma, double *pi, double *k, double *total);
RecoilData CalRecoilK(double *pip, double *pim, double *gamma, double *pi, double *k, double *total);