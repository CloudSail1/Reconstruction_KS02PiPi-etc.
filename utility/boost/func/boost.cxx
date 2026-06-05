#include "../func.h"
GammaData BoostGammaIntoPiPiFrame(double* gamma, double* pip, double* pim) {
    std::vector<double*> frame_fractions = {pip, pim};
    return BoostGammaInto(gamma, pip, pim, frame_fractions);
}
GammaData BoostGammaIntoPiPiGammaFrame(double* gamma, double* pip, double* pim) {
    std::vector<double*> frame_fractions = {pip, pim, gamma};
    return BoostGammaInto(gamma, pip, pim, frame_fractions);
}
GammaData BoostGammaIntoLabFrame(double* gamma, double* pip, double* pim, double* lab) {
    std::vector<double*> frame_fractions = {lab};
    return BoostGammaInto(gamma, pip, pim, frame_fractions, true);
}
GammaData BoostGammaInto(double* gamma, double* pip, double* pim, std::vector<double*> frame_fractions, bool is_inverse) {
    LorentzVector p4_frame;
    for (const auto& raw_fraction : frame_fractions) {
        LorentzVector p4_fraction;
        p4_fraction.SetCoordinates(raw_fraction);
        p4_frame += p4_fraction;
    }
    Boost boost_frame(p4_frame.BoostToCM());
    if (is_inverse) {boost_frame.Inverse();}
    LorentzVector p4_gamma;
    LorentzVector p4_pip;
    LorentzVector p4_pim;
    p4_gamma.SetCoordinates(gamma);
    p4_pip.  SetCoordinates(pip);
    p4_pim.  SetCoordinates(pim);
    p4_gamma = boost_frame (p4_gamma);
    p4_pip   = boost_frame (p4_pip);
    p4_pim   = boost_frame (p4_pim);
    return GammaData {
        p4_gamma.Vect().Unit().Dot(p4_pip.Vect().Unit()),
        p4_gamma.Vect().Unit().Dot(p4_pim.Vect().Unit()),
        p4_gamma.E(),
        p4_gamma.P(),
        p4_gamma,
        p4_pip,
        p4_pim,
    };
}