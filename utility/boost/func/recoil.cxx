#include "../func.h"
RecoilData CalRecoil(std::vector<double*> terms, double* total) {
    LorentzVector sum_terms;
    for (const auto& term : terms) {
        LorentzVector p4_term;
        p4_term.SetCoordinates(term);
        sum_terms += p4_term;
    }
    LorentzVector p4_total;
    p4_total.SetCoordinates(total);
    LorentzVector p4_recoil = p4_total-sum_terms;
    return RecoilData {
        p4_recoil.M(),
        p4_recoil.E(),
        p4_recoil.P(),
        p4_recoil,
    };
}

RecoilData CalRecoilPip(double *pip, double *pim, double *gamma, double *pi, double *k, double *total) {
    return CalRecoil(std::vector{pim,gamma,pi,k}, total);
}
RecoilData CalRecoilPipPim(double *pip, double *pim, double *gamma, double *pi, double *k, double *total) {
    return CalRecoil(std::vector{gamma,pi,k}, total);
}
RecoilData CalRecoilPi(double *pip, double *pim, double *gamma, double *pi, double *k, double *total) {
    return CalRecoil(std::vector{pip,pim,gamma,k}, total);
}
RecoilData CalRecoilK(double *pip, double *pim, double *gamma, double *pi, double *k, double *total) {
    return CalRecoil(std::vector{pip,pim,gamma,pi}, total);
}