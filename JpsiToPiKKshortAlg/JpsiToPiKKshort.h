#ifndef Physics_Analysis_JpsiToPiKKshort_H
#define Physics_Analysis_JpsiToPiKKshort_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h" 
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

class JpsiToPiKKshort : public Algorithm
{
 public:
   JpsiToPiKKshort(const std::string& name, ISvcLocator* pSvcLocator);
   StatusCode initialize();
   StatusCode execute();
   StatusCode finalize();
   void corgen(HepMatrix &, HepVector &, int );
   void corset(HepSymMatrix &, HepMatrix &, int );
   void calibration(RecMdcKalTrack * , HepVector &, int );   
 private:
   int m_evnt;
   double m_ecms;
   double m_vr0cut;
   double m_vz0cut;
   double m_pass[10];
   double m_chisqCut;
   double m_energyThreshold;
   double m_massRangeLower;
   double m_massRangeHigher;
   double m_gammaPhiCut;
   double m_gammaThetaCut;
   double m_gammaAngCut;
   double m_Endcap_th_1;
   double m_Endcap_th_2;
   double m_energyThreshold_b;
   double m_energyThreshold_e;
   double m_Barrel_th;
   bool m_checkMC;
   int m_idNo;
   std::string m_log;
   NTuple::Tuple*  m_tuple_angular;
   NTuple::Tuple*  m_tuple;
   
   NTuple::Array<double>	m_mctruth_P4_Pip_Kshort;
   NTuple::Array<double>   	m_mctruth_P4_Pin_Kshort;
   NTuple::Array<double>   	m_mctruth_P4_Pin_Kstar;
   NTuple::Array<double>   	m_mctruth_P4_Pip_Kstar;
   NTuple::Array<double>   	m_mctruth_P4_Kp_Kstar;
   NTuple::Array<double>   	m_mctruth_P4_Kn_Kstar;
   NTuple::Array<double>   	m_mctruth_P4_Kshort;
   NTuple::Array<double>   	m_mctruth_P4_Kstar;

   NTuple::Item<double>		m_mctruth_P_Pip_Kshort;
   NTuple::Item<double>   	m_mctruth_P_Pin_Kshort;
   NTuple::Item<double>   	m_mctruth_P_Pin_Kstar;
   NTuple::Item<double>   	m_mctruth_P_Pip_Kstar;
   NTuple::Item<double>   	m_mctruth_P_Kp_Kstar;
   NTuple::Item<double>   	m_mctruth_P_Kn_Kstar;
   NTuple::Item<double>   	m_mctruth_P_Kshort;
   NTuple::Item<double>   	m_mctruth_P_Kstar;

   NTuple::Item<double>   	m_mctruth_M_Pip_Kshort;
   NTuple::Item<double>   	m_mctruth_M_Pin_Kshort;
   NTuple::Item<double>   	m_mctruth_M_Pin_Kstar;
   NTuple::Item<double>   	m_mctruth_M_Pip_Kstar;
   NTuple::Item<double>   	m_mctruth_M_Kp_Kstar;
   NTuple::Item<double>   	m_mctruth_M_Kn_Kstar;
   NTuple::Item<double>   	m_mctruth_M_Kshort;
   NTuple::Item<double>   	m_mctruth_M_Kstar;

   NTuple::Item<double>   	m_mctruth_cos_Pip_Kshort;
   NTuple::Item<double>   	m_mctruth_cos_Pin_Kshort;
   NTuple::Item<double>   	m_mctruth_cos_Pin_Kstar;
   NTuple::Item<double>   	m_mctruth_cos_Pip_Kstar;
   NTuple::Item<double>   	m_mctruth_cos_Kp_Kstar;
   NTuple::Item<double>   	m_mctruth_cos_Kn_Kstar;
   NTuple::Item<double>   	m_mctruth_cos_Kshort;
   NTuple::Item<double>   	m_mctruth_cos_Kstar;



   NTuple::Array<double>   mctrue_eta0_P4;
   NTuple::Array<double>   mctrue_eta1_P4;
   NTuple::Array<double>   mctrue_gam0_eta0_P4;
   NTuple::Array<double>   mctrue_gam1_eta0_P4;
   NTuple::Array<double>   mctrue_gam0_eta1_P4;
   NTuple::Array<double>   mctrue_gam1_eta1_P4;
   NTuple::Array<double>   mctrue_Ppp_P4;
   NTuple::Array<double>   mctrue_Ppm_P4;

   NTuple::Item<double>    mctrue_Peta0;
   NTuple::Item<double>    mctrue_Peta1;
   NTuple::Item<double>    mctrue_Pgam0_eta0;
   NTuple::Item<double>    mctrue_Pgam1_eta0;
   NTuple::Item<double>    mctrue_Pgam0_eta1;
   NTuple::Item<double>    mctrue_Pgam1_eta1;
   NTuple::Item<double>    mctrue_Ppp;
   NTuple::Item<double>    mctrue_Ppm;
   NTuple::Item<double>    mctrue_meta0;
   NTuple::Item<double>    mctrue_meta1;
   NTuple::Item<double>    mctrue_mgam0_eta0;
   NTuple::Item<double>    mctrue_mgam1_eta0;
   NTuple::Item<double>    mctrue_mgam0_eta1;
   NTuple::Item<double>    mctrue_mgam1_eta1;
   NTuple::Item<double>    mctrue_mpp;
   NTuple::Item<double>    mctrue_mpm;

   NTuple::Item<double>    m_mctruth_cos_gam1_Eta0;
   NTuple::Item<double>    m_mctruth_cos_gam2_Eta0;
   NTuple::Item<double>    m_mctruth_cos_gam1_Eta1;
   NTuple::Item<double>    m_mctruth_cos_gam2_Eta1;
   NTuple::Item<double>    m_mctruth_cos_Eta0;
   NTuple::Item<double>    m_mctruth_cos_Eta1;
   NTuple::Item<double>    m_mctruth_cos_pp;
   NTuple::Item<double>    m_mctruth_cos_pm;


   NTuple::Item<long>  m_cc;
   NTuple::Array<double>  m_vx0;
   NTuple::Array<double>  m_vy0;
   NTuple::Array<double>  m_vz0;
   NTuple::Array<double>  m_vr0;
   NTuple::Array<double>  m_rvz0;
   NTuple::Array<double>  m_rvxy0;
   NTuple::Array<double>  m_costheta;
   NTuple::Array<double>  m_rvphi0;

   NTuple::Item<long>  m_nn;
   NTuple::Array<double>  m_time;
   NTuple::Array<double>  m_eraw;
   NTuple::Array<double>  m_dang;
   NTuple::Item<double>  m_dang2;
   NTuple::Item<double>  m_dang_proton;
   NTuple::Array<double>  m_costhetaN;

   NTuple::Item<int>   m_ntotCh;
   NTuple::Item<int>   m_ntotNe;
   NTuple::Item<int>   m_run;
   NTuple::Item<int>   m_rec;
   NTuple::Item<int>   m_idxmc;
   NTuple::Array<int>  m_pdgid;
   NTuple::Array<int>  m_trkIdx;
   NTuple::Array<int>  m_motheridx;  


   NTuple::Item<long>  m_nGam;
   
   NTuple::Item<double>  m_chi2;

   NTuple::Array<double>  m_p_total;
   NTuple::Array<double>  m_p_Pip_Kshort;
   NTuple::Array<double>  m_p_Pin_Kshort;
   NTuple::Array<double>  m_p_Pin_Kstar;
   NTuple::Array<double>  m_p_Pip_Kstar;
   NTuple::Array<double>  m_p_Kp_Kstar;
   NTuple::Array<double>  m_p_Kn_Kstar;
   NTuple::Array<double>  m_p_Kshort;
   NTuple::Array<double>  m_p_Kstar;

   NTuple::Item<double>   m_d3_Pip_Kshort;
   NTuple::Item<double>   m_d3_Pin_Kshort;
   NTuple::Item<double>   m_d3_Pin_Kstar;
   NTuple::Item<double>   m_d3_Pip_Kstar;
   NTuple::Item<double>   m_d3_Kp_Kstar;
   NTuple::Item<double>   m_d3_Kn_Kstar;
   NTuple::Item<double>   m_d3_Kshort;
   NTuple::Item<double>   m_d3_Kstar;

   NTuple::Item<double>   m_m_Pip_Kshort;
   NTuple::Item<double>   m_m_Pin_Kshort;
   NTuple::Item<double>   m_m_Pin_Kstar;
   NTuple::Item<double>   m_m_Pip_Kstar;
   NTuple::Item<double>   m_m_Kp_Kstar;
   NTuple::Item<double>   m_m_Kn_Kstar;
   NTuple::Item<double>   m_m_Kshort;
   NTuple::Item<double>   m_m_Kstar;

   NTuple::Item<double>   m_cos_Pip_Kshort	;
   NTuple::Item<double>   m_cos_Pin_Kshort	;
   NTuple::Item<double>   m_cos_Pin_Kstar	;
   NTuple::Item<double>   m_cos_Pip_Kstar	;
   NTuple::Item<double>   m_cos_Kp_Kstar	;
   NTuple::Item<double>   m_cos_Kn_Kstar	;
   NTuple::Item<double>   m_cos_Kshort		;
   NTuple::Item<double>   m_cos_Kstar		;

   NTuple::Array<double>  m_p_other_Kshort;
   NTuple::Item<double>   m_d3_other_Kshort;
   NTuple::Item<double>   m_m_other_Kshort;
   NTuple::Item<double>   m_cos_other_Kshort;
   NTuple::Item<double>   m_bkg_m_other_Kshort; 

   NTuple::Array<double>  m_p_Pip_other_Kshort;
   NTuple::Item<double>   m_d3_Pip_other_Kshort;
   NTuple::Item<double>   m_m_Pip_other_Kshort;
   NTuple::Item<double>   m_cos_Pip_other_Kshort;
   NTuple::Item<double>   m_bkg_m_Pip_other_Kshort; 
 
   NTuple::Array<double>  m_p_Pin_other_Kshort;
   NTuple::Item<double>   m_d3_Pin_other_Kshort;
   NTuple::Item<double>   m_m_Pin_other_Kshort;
   NTuple::Item<double>   m_cos_Pin_other_Kshort;
   NTuple::Item<double>   m_bkg_m_Pin_other_Kshort; 

   NTuple::Array<double>  m_p_other_Kstar;
   NTuple::Item<double>   m_d3_other_Kstar;
   NTuple::Item<double>   m_m_other_Kstar;
   NTuple::Item<double>   m_cos_other_Kstar;
   NTuple::Item<double>   m_bkg_m_other_Kstar; 

   NTuple::Array<double>  m_p_Kp_other_Kstar;
   NTuple::Item<double>   m_d3_Kp_other_Kstar;
   NTuple::Item<double>   m_m_Kp_other_Kstar;
   NTuple::Item<double>   m_cos_Kp_other_Kstar;
   NTuple::Item<double>   m_bkg_m_Kp_other_Kstar; 
 
   NTuple::Array<double>  m_p_Pin_other_Kstar;
   NTuple::Item<double>   m_d3_Pin_other_Kstar;
   NTuple::Item<double>   m_m_Pin_other_Kstar;
   NTuple::Item<double>   m_cos_Pin_other_Kstar;
   NTuple::Item<double>   m_bkg_m_Pin_other_Kstar; 

   NTuple::Item<double>   m_chisq_Kshort;
   NTuple::Item<double>   m_chisq_kmfit;
   NTuple::Item<double>   m_chisq_other_Kshort;
   NTuple::Item<double>   m_chisq_other_Kstar;
   NTuple::Item<double>   m_ratio_Kshort;
   NTuple::Item<double>   m_ratio_other;
   NTuple::Item<double>   m_chisq_km_Kshort;


   NTuple::Item<double>  m_bkg_m_Pip_Kshort     ; 
   NTuple::Item<double>  m_bkg_m_Pin_Kshort     ; 
   NTuple::Item<double>  m_bkg_m_Pin_Kstar      ; 
   NTuple::Item<double>  m_bkg_m_Pip_Kstar      ; 
   NTuple::Item<double>  m_bkg_m_Kp_Kstar       ; 
   NTuple::Item<double>  m_bkg_m_Kn_Kstar       ; 
   NTuple::Item<double>  m_bkg_m_Kshort         ; 
   NTuple::Item<double>  m_bkg_m_Kstar          ; 
   NTuple::Item<double>  m_bkg_m_Pin_Pin        ; 
   NTuple::Item<double>  m_bkg_m_Pip_Kp         ; 
   NTuple::Item<double>  m_bkg_m_Pin_Kp         ; 
   NTuple::Item<double>  m_bkg_m_Pip_Pin        ; 



   NTuple::Item<double>  m_mPi0_bkg0_6C;
   NTuple::Item<double>  m_mPi0_bkg1_6C;

   NTuple::Item<double>  m_mPI0_bkg_6C;
   NTuple::Item<double>  m_mETA_bkg_6C;



   NTuple::Item<double>  m_chi2_6C;
   NTuple::Item<double>  m_chi2_2pi0_6C;
   NTuple::Item<double>  m_chi2_pi0Eta_6C;




   NTuple::Item<double>  m_m12g_4C;
   NTuple::Item<double>  m_m13g_4C;
   NTuple::Item<double>  m_m14g_4C;
   NTuple::Item<double>  m_m23g_4C;
   NTuple::Item<double>  m_m24g_4C;
   NTuple::Item<double>  m_m34g_4C;

   NTuple::Item<double>  m_cosEta0;
   NTuple::Item<double>  m_cosEta1;


   NTuple::Item<double>  m_chi23g;
   NTuple::Item<double>  m_chi25g;
   NTuple::Item<double>  m_chi26g;

   NTuple::Item<double>  m_mEta0_4C;
   NTuple::Item<double>  m_mEta1_4C;
   NTuple::Item<double>  m_m2Etarec_bkg;

   NTuple::Item<double>  m_pEta0_4C;
   NTuple::Item<double>  m_pEta1_4C;
   NTuple::Item<double>  m_eEta0_4C;
   NTuple::Item<double>  m_eEta1_4C;
  
   NTuple::Item<double>  m_mgamrec1_bkg;
   NTuple::Item<double>  m_mgamrec2_bkg;
   NTuple::Item<double>  m_mgamrec3_bkg;
   NTuple::Item<double>  m_mgamrec4_bkg;
   
   NTuple::Item<double>  m_eGam;
   
   NTuple::Item<double>  m_eraw1;
   NTuple::Item<double>  m_eraw2;
   NTuple::Item<double>  m_eoverp1;
   NTuple::Item<double>  m_eoverp2;
   NTuple::Item<double>  m_ptot1;
   NTuple::Item<double>  m_ptot2;
   NTuple::Item<double>  m_depth1;
   NTuple::Item<double>  m_depth2;
  
   NTuple::Item<double>  m_npip;
   NTuple::Item<double>  m_npin;
   NTuple::Item<double>  m_nkp;
   NTuple::Item<double>  m_nkn;
 
   NTuple::Item<double>  m_npp;// add for PID
   NTuple::Item<double>  m_npm;//PID
   //NTuple::Item<double>  m_npim;//PID
   //NTuple::Item<double>  m_npip;//PID
   
  ///for bkg
   
   NTuple::Item<double>  m_mEtarec0_bkg;
   NTuple::Item<double>  m_mEtarec1_bkg;

   NTuple::Item<double>  m_m2gamrec1_bkg;
   NTuple::Item<double>  m_m2gamrec2_bkg;
   NTuple::Item<double>  m_m2gamrec3_bkg;
   NTuple::Item<double>  m_m2gamrec4_bkg;
   NTuple::Item<double>  m_m2gamrec5_bkg;
   NTuple::Item<double>  m_m2gamrec6_bkg;



   NTuple::Item<double>  m_mpEta0;
   NTuple::Item<double>  m_mpEta1;
   NTuple::Item<double>  m_mpbEta0;
   NTuple::Item<double>  m_mpbEta1;
   NTuple::Item<double>  m_mppbEta0;
   NTuple::Item<double>  m_mppbEta1;


   NTuple::Array<double>  m_P_tot;
   NTuple::Array<double>  m_P_pp;
   NTuple::Array<double>  m_P_pm;
   NTuple::Array<double>  m_P_eta0;
   NTuple::Array<double>  m_P_eta1;
   NTuple::Array<double>  m_P_gam0_eta0;
   NTuple::Array<double>  m_P_gam1_eta0;
   NTuple::Array<double>  m_P_gam0_eta1;
   NTuple::Array<double>  m_P_gam1_eta1;

   NTuple::Item<double>  m_mPP;
   NTuple::Item<double>  m_mPM;
   NTuple::Item<double>  m_d3_pp;
   NTuple::Item<double>  m_d3_pm;

   NTuple::Item<double>  m_d3_eta0;
   NTuple::Item<double>  m_d3_eta1;
   NTuple::Item<double>  m_d3_gam0_eta0;
   NTuple::Item<double>  m_d3_gam1_eta0;
   NTuple::Item<double>  m_d3_gam0_eta1;
   NTuple::Item<double>  m_d3_gam1_eta1;

// rm 4gamma
//   NTuple::Item<double>  m_mg13_wc;
//   NTuple::Item<double>  m_mg14_wc;
//   NTuple::Item<double>  m_mg23_wc;
//   NTuple::Item<double>  m_mg24_wc;

   NTuple::Item<double>  m_mgamppbEta0_eta0;
   NTuple::Item<double>  m_mgamppbEta1_eta0;
   NTuple::Item<double>  m_mgamppbEta0_eta1;
   NTuple::Item<double>  m_mgamppbEta1_eta1;


   NTuple::Item<double>  m_mppb;
   NTuple::Item<double>  m_mppb2Eta;


///////////  2020-01-08
   NTuple::Item<double>  m_m2Eta;
   NTuple::Item<double>   m_mpEta;
   NTuple::Item<double>   m_mpbEta;
  

   NTuple::Item<double>  m_cos_pp;
   NTuple::Item<double>  m_cos_pm;
   NTuple::Item<double>  m_cos_Eta0;
   NTuple::Item<double>  m_cos_Eta1;
   NTuple::Item<double>  m_cos_gam1_Eta0;
   NTuple::Item<double>  m_cos_gam2_Eta0;
   NTuple::Item<double>  m_cos_gam1_Eta1;
   NTuple::Item<double>  m_cos_gam2_Eta1;
};
#endif
