#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "VertexFit/Helix.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "McTruth/McParticle.h"
#include "VertexFit/IVertexDbSvc.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TRandom.h" 
using CLHEP::Hep3Vector;
#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep2Vector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "ParticleID/ParticleID.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/SecondVertexFit.h"
#include <vector>
#include "JpsiToPiKKshortAlg/Comb.h"
#include "JpsiToPiKKshortAlg/JpsiToPiKKshort.h"
typedef std::vector<int> Vint;
typedef std::vector<Comb> VComb;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
const double mpi = 0.13497;
const double xmass[5] = 
{
	0.000511, 0.105658, 0.139570, 0.493677, 0.938272
}
;
bool debug = false;
const double mka = xmass[3];
int n_total_event,n_charge_cutted,n_charged_partical_cutted,Ncut2,n_particle_number_cutted,Ncut4,Ncut5,Ncut6,Ncut7,Ncut8,Ncut9,Ncut10,n_pip_pin_vtxfit_cutted,n_vtxfit_kshort_cutted,n_vtxfit_kstar_cutted,n_kmfit_cutted,n_mass_cutted,n_vtxfit_other_Kshort_ok,n_vtxfit_other_Kstar_ok;
////////////////////////////////////////////////////////////////////
JpsiToPiKKshort::JpsiToPiKKshort(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) 
{  
	declareProperty("CmsEnergy", m_ecms = 3.097 );
	declareProperty("Vr0cut", m_vr0cut=1.0);
	declareProperty("Vz0cut", m_vz0cut=10.0);
	declareProperty("EnergyThreshold", m_energyThreshold=0.02);
	declareProperty("EnergyThreshold_b", m_energyThreshold_b=0.025);
	declareProperty("EnergyThreshold_e", m_energyThreshold_e=0.05);
	declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
	declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
	declareProperty("GammaAngCut", m_gammaAngCut=20.0);
	declareProperty("idNo",m_idNo = 443);
	//Declare the properties  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode JpsiToPiKKshort::initialize()
{
	MsgStream log(msgSvc(), name());

	log << MSG::INFO << "in initialize()" << endmsg;

	StatusCode status;
	log << MSG::DEBUG << __LINE__ << " mctruth ntuple init " << endmsg;
	NTuplePtr nt_mctruth(ntupleSvc(), "FILE1/mctruth");
	if (nt_mctruth) {
		m_tuple_angular = nt_mctruth;
		log << MSG::DEBUG << __LINE__ << endmsg;
	}
	else
	{
		log << MSG::DEBUG << __LINE__ << endmsg;
		m_tuple_angular = ntupleSvc() -> book("FILE1/mctruth", CLID_ColumnWiseTuple,"ks N-tuple");

		if(m_tuple_angular)
		{
			status = m_tuple_angular->addItem("eta0_P4",      4,       mctrue_eta0_P4);
			status = m_tuple_angular->addItem("eta1_P4",      4,       mctrue_eta1_P4);
			status = m_tuple_angular->addItem("gam0_eta0_P4",      4,       mctrue_gam0_eta0_P4);
			status = m_tuple_angular->addItem("gam1_eta0_P4",      4,       mctrue_gam1_eta0_P4);
			status = m_tuple_angular->addItem("gam0_eta1_P4",      4,       mctrue_gam0_eta1_P4);
			status = m_tuple_angular->addItem("gam1_eta1_P4",      4,       mctrue_gam1_eta1_P4);
			status = m_tuple_angular->addItem("Ppp_P4",      4,       mctrue_Ppp_P4);
			status = m_tuple_angular->addItem("Ppm_P4",      4,       mctrue_Ppm_P4);

			status = m_tuple_angular->addItem("Peta0",        mctrue_Peta0);
			status = m_tuple_angular->addItem("Peta1",        mctrue_Peta1);
			status = m_tuple_angular->addItem("Pgam0_eta0",        mctrue_Pgam0_eta0);
			status = m_tuple_angular->addItem("Pgam1_eta0",        mctrue_Pgam1_eta0);
			status = m_tuple_angular->addItem("Pgam0_eta1",        mctrue_Pgam0_eta1);
			status = m_tuple_angular->addItem("Pgam1_eta1",        mctrue_Pgam1_eta1);
			status = m_tuple_angular->addItem("Ppp",        mctrue_Ppp);
			status = m_tuple_angular->addItem("Ppm",        mctrue_Ppm);


			status = m_tuple_angular->addItem("meta0",        mctrue_meta0);
			status = m_tuple_angular->addItem("meta1",        mctrue_meta1);
			status = m_tuple_angular->addItem("mgam0_eta0",        mctrue_mgam0_eta0);
			status = m_tuple_angular->addItem("mgam1_eta0",        mctrue_mgam1_eta0);
			status = m_tuple_angular->addItem("mgam0_eta1",        mctrue_mgam0_eta1);
			status = m_tuple_angular->addItem("mgam1_eta1",        mctrue_mgam1_eta1);
			status = m_tuple_angular->addItem("mpp",        mctrue_mpp);
			status = m_tuple_angular->addItem("mpm",        mctrue_mpm);

                        status = m_tuple_angular->addItem("mctruth_cos_gam1_Eta0",     m_mctruth_cos_gam1_Eta0);
                        status = m_tuple_angular->addItem("mctruth_cos_gam2_Eta0",     m_mctruth_cos_gam2_Eta0);
                        status = m_tuple_angular->addItem("mctruth_cos_gam1_Eta1",     m_mctruth_cos_gam1_Eta1);
                        status = m_tuple_angular->addItem("mctruth_cos_gam2_Eta1",     m_mctruth_cos_gam2_Eta1);

                        status = m_tuple_angular->addItem("mctruth_cos_Eta0",     m_mctruth_cos_Eta0);
                        status = m_tuple_angular->addItem("mctruth_cos_Eta1",     m_mctruth_cos_Eta1);
                        status = m_tuple_angular->addItem("mctruth_cos_pp",     m_mctruth_cos_pp);
                        status = m_tuple_angular->addItem("mctruth_cos_pm",     m_mctruth_cos_pm);
			
                        status = m_tuple_angular->addItem("m_mctruth_P_Pip_Kshort",	m_mctruth_P_Pip_Kshort	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Pin_Kshort",	m_mctruth_P_Pin_Kshort	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Pin_Kstar",	m_mctruth_P_Pin_Kstar	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Pip_Kstar",	m_mctruth_P_Pip_Kstar	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Kp_Kstar",	m_mctruth_P_Kp_Kstar	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Kn_Kstar",	m_mctruth_P_Kn_Kstar	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Kshort",		m_mctruth_P_Kshort	);
                        status = m_tuple_angular->addItem("m_mctruth_P_Kstar",		m_mctruth_P_Kstar	);
	
                        status = m_tuple_angular->addItem("mctruth_M_Pip_Kshort",	m_mctruth_M_Pip_Kshort	);
                        status = m_tuple_angular->addItem("mctruth_M_Pin_Kshort",	m_mctruth_M_Pin_Kshort	);
                        status = m_tuple_angular->addItem("mctruth_M_Pin_Kstar",	m_mctruth_M_Pin_Kstar	);
                        status = m_tuple_angular->addItem("mctruth_M_Pip_Kstar",	m_mctruth_M_Pip_Kstar	);
                        status = m_tuple_angular->addItem("mctruth_M_Kp_Kstar",		m_mctruth_M_Kp_Kstar	);
                        status = m_tuple_angular->addItem("mctruth_M_Kn_Kstar",		m_mctruth_M_Kn_Kstar	);
                        status = m_tuple_angular->addItem("mctruth_M_Kshort",		m_mctruth_M_Kshort	);
                        status = m_tuple_angular->addItem("mctruth_M_Kstar",		m_mctruth_M_Kstar	);
                        
			status = m_tuple_angular->addItem("mctruth_cos_Pip_Kshort",	m_mctruth_cos_Pip_Kshort);
			status = m_tuple_angular->addItem("mctruth_cos_Pin_Kshort",	m_mctruth_cos_Pin_Kshort);
			status = m_tuple_angular->addItem("mctruth_cos_Pin_Kstar",	m_mctruth_cos_Pin_Kstar	);
			status = m_tuple_angular->addItem("mctruth_cos_Pip_Kstar",	m_mctruth_cos_Pip_Kstar	);
			status = m_tuple_angular->addItem("mctruth_cos_Kp_Kstar",	m_mctruth_cos_Kp_Kstar	);
			status = m_tuple_angular->addItem("mctruth_cos_Kn_Kstar",	m_mctruth_cos_Kn_Kstar	);
			status = m_tuple_angular->addItem("mctruth_cos_Kshort",		m_mctruth_cos_Kshort	);
			status = m_tuple_angular->addItem("mctruth_cos_Kstar",		m_mctruth_cos_Kstar	);
			
			status = m_tuple_angular->addItem("mctruth_P4_Pip_Kshort",	4,m_mctruth_P4_Pip_Kshort	);
			status = m_tuple_angular->addItem("mctruth_P4_Pin_Kshort",	4,m_mctruth_P4_Pin_Kshort	);
			status = m_tuple_angular->addItem("mctruth_P4_Pin_Kstar",	4,m_mctruth_P4_Pin_Kstar	);
			status = m_tuple_angular->addItem("mctruth_P4_Pip_Kstar",	4,m_mctruth_P4_Pip_Kstar	);
			status = m_tuple_angular->addItem("mctruth_P4_Kp_Kstar",	4,m_mctruth_P4_Kp_Kstar		);
			status = m_tuple_angular->addItem("mctruth_P4_Kn_Kstar",	4,m_mctruth_P4_Kn_Kstar		);
			status = m_tuple_angular->addItem("mctruth_P4_Kshort",		4,m_mctruth_P4_Kshort		);
			status = m_tuple_angular->addItem("mctruth_P4_Kstar",		4,m_mctruth_P4_Kstar		);
		}
		else
		{
			log << MSG::FATAL << "Cannot book N-tuple: " << long(m_tuple_angular) << endmsg;

			return StatusCode::FAILURE;
		}
	}

	log << MSG::DEBUG << __LINE__ << " single ntuple init" << endmsg;
	NTuplePtr nt2(ntupleSvc(), "FILE1/chi2omg");
	if ( nt2 ) m_tuple = nt2;
	else 
	{
		m_tuple = ntupleSvc()->book ("FILE1/chi2omg", CLID_ColumnWiseTuple,"chi2omg");
		if ( m_tuple )    
		{
			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("cc",        m_cc, 0,500);
			status = m_tuple->addItem("vx0",       m_cc,  m_vx0);
			status = m_tuple->addItem("vy0",       m_cc,  m_vy0);
			status = m_tuple->addItem("vz0",       m_cc,  m_vz0);
			status = m_tuple->addItem("vr0",       m_cc,  m_vr0);
			status = m_tuple->addItem("rvz0",      m_cc,  m_rvz0);
			status = m_tuple->addItem("rvxy0",     m_cc,  m_rvxy0);
			status = m_tuple->addItem("costheta",  m_cc,  m_costheta);
			status = m_tuple->addItem("rvphi0",    m_cc,  m_rvphi0);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("nn",        m_nn, 0,200);
			status = m_tuple->addItem("time",      m_nn, m_time);
			status = m_tuple->addItem("costhetaN", m_nn, m_costhetaN);
			status = m_tuple->addItem("eraw",      m_nn, m_eraw);
			status = m_tuple->addItem("dang",      m_nn, m_dang);
			status = m_tuple->addItem("dang2",     m_dang2);
			status = m_tuple->addItem("dang_proton",     m_dang_proton);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("ntotCh",  m_ntotCh );
			status = m_tuple->addItem("ntotNe",  m_ntotNe );

			status = m_tuple->addItem("run",  m_run );
			status = m_tuple->addItem("rec",  m_rec );
			status = m_tuple->addItem("indexmc",   m_idxmc, 0, 100);
			status = m_tuple->addItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple->addItem("trkIdx",     m_idxmc, m_trkIdx);
			status = m_tuple->addItem("motheridx", m_idxmc, m_motheridx);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("nGam",      m_nGam);


			status = m_tuple->addItem("chi2",      m_chi2);

			status = m_tuple->addItem("mPi0_bkg0_6C",      m_mPi0_bkg0_6C);
			status = m_tuple->addItem("mPi0_bkg1_6C",      m_mPi0_bkg1_6C);
			status = m_tuple->addItem("mPI0_bkg_6C",      m_mPI0_bkg_6C);
			status = m_tuple->addItem("mETA_bkg_6C",      m_mETA_bkg_6C);



			status = m_tuple->addItem("chi2_6C",      m_chi2_6C);
			status = m_tuple->addItem("chi2_2pi0_6C",      m_chi2_2pi0_6C);
			status = m_tuple->addItem("chi2_pi0Eta_6C",      m_chi2_pi0Eta_6C);


			log << MSG::DEBUG << __LINE__ << endmsg;

			status = m_tuple->addItem("m12g_4C",      m_m12g_4C);
			status = m_tuple->addItem("m13g_4C",      m_m13g_4C);
			status = m_tuple->addItem("m14g_4C",      m_m14g_4C);
			status = m_tuple->addItem("m23g_4C",      m_m23g_4C);
			status = m_tuple->addItem("m24g_4C",      m_m24g_4C);
			status = m_tuple->addItem("m34g_4C",      m_m34g_4C);

			status = m_tuple->addItem("cosEta0",     m_cosEta0);
			status = m_tuple->addItem("cosEta1",     m_cosEta1);


			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("P_tot",  4,   m_P_tot);
			status = m_tuple->addItem("P_pp",  4,   m_P_pp);
			status = m_tuple->addItem("P_pm",  4,   m_P_pm);
			status = m_tuple->addItem("P_eta0",  4,   m_P_eta0);
			status = m_tuple->addItem("P_eta1",  4,   m_P_eta1);
			status = m_tuple->addItem("P_gam0_eta0",  4,   m_P_gam0_eta0);
			status = m_tuple->addItem("P_gam1_eta0",  4,   m_P_gam1_eta0);
			status = m_tuple->addItem("P_gam0_eta1",  4,   m_P_gam0_eta1);
			status = m_tuple->addItem("P_gam1_eta1",  4,   m_P_gam1_eta1);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("mPP",   m_mPP);
			status = m_tuple->addItem("mPM",   m_mPM);
			status = m_tuple->addItem("d3_pp",   m_d3_pp);
			status = m_tuple->addItem("d3_pm",   m_d3_pm);
			status = m_tuple->addItem("d3_gam0_eta0",   m_d3_gam0_eta0);
			status = m_tuple->addItem("d3_gam1_eta0",   m_d3_gam1_eta0);
			status = m_tuple->addItem("d3_gam0_eta1",   m_d3_gam0_eta1);
			status = m_tuple->addItem("d3_gam1_eta1",   m_d3_gam1_eta1);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("d3_eta0",   m_d3_eta0);
			status = m_tuple->addItem("d3_eta1",   m_d3_eta1);
/*
			status = m_tuple->addItem("mg13_wc",   m_mg13_wc);
			status = m_tuple->addItem("mg14_wc",   m_mg14_wc);
			status = m_tuple->addItem("mg23_wc",   m_mg23_wc);
			status = m_tuple->addItem("mg24_wc",   m_mg24_wc);
*/


			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("chi23g",    m_chi23g);
			status = m_tuple->addItem("chi25g",    m_chi25g);
			status = m_tuple->addItem("chi26g",    m_chi26g);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("eraw1",     m_eraw1);
			status = m_tuple->addItem("eraw2",     m_eraw2);
			status = m_tuple->addItem("eoverp1",     m_eoverp1);
			status = m_tuple->addItem("eoverp2",     m_eoverp2);
			status = m_tuple->addItem("ptot1",     m_ptot1);
			status = m_tuple->addItem("ptot2",     m_ptot2);
			status = m_tuple->addItem("depth1",     m_depth1);
			status = m_tuple->addItem("depth2",     m_depth2);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("npp",     m_npp);//add for PID
			status = m_tuple->addItem("npm",     m_npm);//add for PID
			//status = m_tuple->addItem("npim",     m_npim);//add for PID
			//status = m_tuple->addItem("npip",     m_npip);//add for PID

			status = m_tuple->addItem("mEta0_4C",     m_mEta0_4C);
			status = m_tuple->addItem("mEta1_4C",     m_mEta1_4C);
			status = m_tuple->addItem("m2Etarec_bkg",     m_m2Etarec_bkg);

			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("eEta0_4C",     m_eEta0_4C);
			status = m_tuple->addItem("eEta1_4C",     m_eEta1_4C);
			status = m_tuple->addItem("pEta0_4C",     m_pEta0_4C);
			status = m_tuple->addItem("pEta1_4C",     m_pEta1_4C);


			status = m_tuple->addItem("mgamrec1_bkg",     m_mgamrec1_bkg);
			status = m_tuple->addItem("mgamrec2_bkg",     m_mgamrec2_bkg);
			status = m_tuple->addItem("mgamrec3_bkg",     m_mgamrec3_bkg);
			status = m_tuple->addItem("mgamrec4_bkg",     m_mgamrec4_bkg);

			log << MSG::DEBUG << __LINE__ << endmsg;
			/////for bkg
			status = m_tuple->addItem("mEtarec0_bkg",      m_mEtarec0_bkg);
			status = m_tuple->addItem("mEtarec1_bkg",      m_mEtarec1_bkg);

			status = m_tuple->addItem("m2gamrec1_bkg",    m_m2gamrec1_bkg);
			status = m_tuple->addItem("m2gamrec2_bkg",    m_m2gamrec2_bkg);
			status = m_tuple->addItem("m2gamrec3_bkg",    m_m2gamrec3_bkg);
			status = m_tuple->addItem("m2gamrec4_bkg",    m_m2gamrec4_bkg);
			status = m_tuple->addItem("m2gamrec5_bkg",    m_m2gamrec5_bkg);
			status = m_tuple->addItem("m2gamrec6_bkg",    m_m2gamrec6_bkg);




			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("mpEta0",      m_mpEta0);
			status = m_tuple->addItem("mpEta1",      m_mpEta1);

			status = m_tuple->addItem("mpbEta0",      m_mpbEta0);
			status = m_tuple->addItem("mpbEta1",      m_mpbEta1);

			status = m_tuple->addItem("mppbEta0",      m_mppbEta0);
			status = m_tuple->addItem("mppbEta1",      m_mppbEta1);

			status = m_tuple->addItem("mgamppbEta0_eta1",    m_mgamppbEta0_eta1);
			status = m_tuple->addItem("mgamppbEta1_eta1",    m_mgamppbEta1_eta1);
			status = m_tuple->addItem("mgamppbEta0_eta0",    m_mgamppbEta0_eta0);
			status = m_tuple->addItem("mgamppbEta1_eta0",    m_mgamppbEta1_eta0);

			status = m_tuple->addItem("mppb",    m_mppb);
			status = m_tuple->addItem("mppb2Eta",    m_mppb2Eta);



			log << MSG::DEBUG << __LINE__ << endmsg;
			////////////////   2020-01-08
			status = m_tuple->addItem("m2Eta",    m_m2Eta);
			status = m_tuple->addItem("mpEta",    m_mpEta);
			status = m_tuple->addItem("mpbEta",    m_mpbEta);

			status = m_tuple->addItem("cos_pp",    m_cos_pp);
			status = m_tuple->addItem("cos_pm",    m_cos_pm);
			status = m_tuple->addItem("cos_Eta0",    m_cos_Eta0);
			status = m_tuple->addItem("cos_Eta1",    m_cos_Eta1);
			status = m_tuple->addItem("cos_gam1_Eta0",    m_cos_gam1_Eta0);
			status = m_tuple->addItem("cos_gam2_Eta0",    m_cos_gam2_Eta0);
			status = m_tuple->addItem("cos_gam1_Eta1",    m_cos_gam1_Eta1);
			status = m_tuple->addItem("cos_gam2_Eta1",    m_cos_gam2_Eta1);

			
			log << MSG::DEBUG << __LINE__ << endmsg;
			status = m_tuple->addItem("n_pip",m_npip);
			status = m_tuple->addItem("n_pin",m_npin);
			status = m_tuple->addItem("n_kp" ,m_nkp );
			status = m_tuple->addItem("n_kn" ,m_nkn );

			status = m_tuple->addItem("chisq_Kshort" ,m_chisq_Kshort);
			status = m_tuple->addItem("chisq_kmfit" ,m_chisq_kmfit);
			status = m_tuple->addItem("chisq_other_Kshort" ,m_chisq_other_Kshort);
			status = m_tuple->addItem("chisq_other_Kstar" ,m_chisq_other_Kstar);
			status = m_tuple->addItem("ratio_Kshort",m_ratio_Kshort);
			status = m_tuple->addItem("ratio_other", m_ratio_other);

			status = m_tuple->addItem("p_total",		4,m_p_total		);
			status = m_tuple->addItem("p_Pip_Kshort",	4,m_p_Pip_Kshort	);
			status = m_tuple->addItem("p_Pin_Kshort",	4,m_p_Pin_Kshort	);
			status = m_tuple->addItem("p_Pin_Kstar",	4,m_p_Pin_Kstar		);
			status = m_tuple->addItem("p_Pip_Kstar",	4,m_p_Pip_Kstar		);
			status = m_tuple->addItem("p_Kp_Kstar",		4,m_p_Kp_Kstar		);
			status = m_tuple->addItem("p_Kn_Kstar",		4,m_p_Kn_Kstar		);
			status = m_tuple->addItem("p_Kshort",		4,m_p_Kshort		);
			status = m_tuple->addItem("p_Kstar",		4,m_p_Kstar		);
			status = m_tuple->addItem("p_Pip_other_Kshort",	4,m_p_Pip_other_Kshort	);
			status = m_tuple->addItem("p_Pin_other_Kshort",	4,m_p_Pin_other_Kshort	);
			status = m_tuple->addItem("p_other_Kshort",	4,m_p_other_Kshort	);
			status = m_tuple->addItem("p_Kp_other_Kstar",	4,m_p_Kp_other_Kstar	);
			status = m_tuple->addItem("p_Pin_other_Kstar",	4,m_p_Pin_other_Kstar	);
			status = m_tuple->addItem("p_other_Kstar",	4,m_p_other_Kstar	);

			status = m_tuple->addItem("d3_Pip_Kshort",	m_d3_Pip_Kshort		);
			status = m_tuple->addItem("d3_Pin_Kshort",	m_d3_Pin_Kshort		);
			status = m_tuple->addItem("d3_Pin_Kstar",	m_d3_Pin_Kstar		);
			status = m_tuple->addItem("d3_Pip_Kstar",	m_d3_Pip_Kstar		);
			status = m_tuple->addItem("d3_Kp_Kstar",	m_d3_Kp_Kstar		);
			status = m_tuple->addItem("d3_Kn_Kstar",	m_d3_Kn_Kstar		);
			status = m_tuple->addItem("d3_Kshort",		m_d3_Kshort		);
			status = m_tuple->addItem("d3_Kstar",		m_d3_Kstar		);
			status = m_tuple->addItem("d3_Pip_other_Kshort",m_d3_Pip_other_Kshort	);
			status = m_tuple->addItem("d3_Pin_other_Kshort",m_d3_Pin_other_Kshort	);
			status = m_tuple->addItem("d3_other_Kshort",	m_d3_other_Kshort	);
			status = m_tuple->addItem("d3_Kp_other_Kstar",	m_d3_Kp_other_Kstar	);
			status = m_tuple->addItem("d3_Pin_other_Kstar",	m_d3_Pin_other_Kstar	);
			status = m_tuple->addItem("d3_other_Kstar",	m_d3_other_Kstar	);

			status = m_tuple->addItem("m_Pip_Kshort",	m_m_Pip_Kshort	);
			status = m_tuple->addItem("m_Pin_Kshort",	m_m_Pin_Kshort	);
			status = m_tuple->addItem("m_Pin_Kstar",	m_m_Pin_Kstar	);
			status = m_tuple->addItem("m_Pip_Kstar",	m_m_Pip_Kstar	);
			status = m_tuple->addItem("m_Kp_Kstar",		m_m_Kp_Kstar	);
			status = m_tuple->addItem("m_Kn_Kstar",		m_m_Kn_Kstar	);
			status = m_tuple->addItem("m_Kshort",		m_m_Kshort	);
			status = m_tuple->addItem("m_Kstar",		m_m_Kstar	);
			status = m_tuple->addItem("m_Pip_other_Kshort",		m_m_Pip_other_Kshort	);
			status = m_tuple->addItem("m_Pin_other_Kshort",		m_m_Pin_other_Kshort	);
			status = m_tuple->addItem("m_other_Kshort",		m_m_other_Kshort	);
			status = m_tuple->addItem("m_Kp_other_Kstar",		m_m_Kp_other_Kstar	);
			status = m_tuple->addItem("m_Pin_other_Kstar",		m_m_Pin_other_Kstar	);
			status = m_tuple->addItem("m_other_Kstar",		m_m_other_Kstar	);

			status = m_tuple->addItem("cos_Pip_Kshort",	m_cos_Pip_Kshort);
			status = m_tuple->addItem("cos_Pin_Kshort",	m_cos_Pin_Kshort);
			status = m_tuple->addItem("cos_Pin_Kstar",	m_cos_Pin_Kstar );
			status = m_tuple->addItem("cos_Pip_Kstar",	m_cos_Pip_Kstar );
			status = m_tuple->addItem("cos_Kp_Kstar",	m_cos_Kp_Kstar  );
			status = m_tuple->addItem("cos_Kn_Kstar",	m_cos_Kn_Kstar  );
			status = m_tuple->addItem("cos_Kshort",		m_cos_Kshort    );
			status = m_tuple->addItem("cos_Kstar",		m_cos_Kstar     );
			status = m_tuple->addItem("cos_Pip_other_Kshort",	m_cos_Pip_other_Kshort     );
			status = m_tuple->addItem("cos_Pin_other_Kshort",	m_cos_Pin_other_Kshort     );
			status = m_tuple->addItem("cos_other_Kshort",		m_cos_other_Kshort     );
			status = m_tuple->addItem("cos_Kp_other_Kstar",		m_cos_Kp_other_Kstar     );
			status = m_tuple->addItem("cos_Pin_other_Kstar",	m_cos_Pin_other_Kstar     );
			status = m_tuple->addItem("cos_other_Kstar",		m_cos_other_Kstar     );

			status = m_tuple->addItem("bkg_m_Pip_kshort",	m_bkg_m_Pip_Kshort);
			status = m_tuple->addItem("bkg_m_Pin_Kshort",	m_bkg_m_Pin_Kshort);
			status = m_tuple->addItem("bkg_m_Pin_Kstar",	m_bkg_m_Pin_Kstar );
			status = m_tuple->addItem("bkg_m_Pip_Kstar",	m_bkg_m_Pip_Kstar );
			status = m_tuple->addItem("bkg_m_Kp_Kstar",	m_bkg_m_Kp_Kstar  );
			status = m_tuple->addItem("bkg_m_Kn_Kstar",	m_bkg_m_Kn_Kstar  );
			status = m_tuple->addItem("bkg_m_Kshort",	m_bkg_m_Kshort    );
			status = m_tuple->addItem("bkg_m_Kstar",	m_bkg_m_Kstar     );
			status = m_tuple->addItem("bkg_m_Pin_Pin",	m_bkg_m_Pin_Pin   );
			status = m_tuple->addItem("bkg_m_Pip_Kp",	m_bkg_m_Pip_Kp    );
			status = m_tuple->addItem("bkg_m_Pin_Kp",	m_bkg_m_Pin_Kp    );
			status = m_tuple->addItem("bkg_m_Pip_Pin",	m_bkg_m_Pip_Pin   );
			status = m_tuple->addItem("bkg_m_Pip_other_Kshort",	m_bkg_m_Pip_other_Kshort     );
			status = m_tuple->addItem("bkg_m_Pin_other_Kshort",	m_bkg_m_Pin_other_Kshort     );
			status = m_tuple->addItem("bkg_m_other_Kshort",		m_bkg_m_other_Kshort     );
			status = m_tuple->addItem("bkg_m_Kp_other_Kstar",	m_bkg_m_Kp_other_Kstar     );
			status = m_tuple->addItem("bkg_m_Pin_other_Kstar",	m_bkg_m_Pin_other_Kstar     );
			status = m_tuple->addItem("bkg_m_other_Kstar",		m_bkg_m_other_Kstar     );

			// for debug use
			status = m_tuple->addItem("m_orig_Kshort",		m_m_orig_Kshort				);
			status = m_tuple->addItem("m_orig_Kstar",		m_m_orig_Kstar				);
			status = m_tuple->addItem("m_vtxadd_Kshort",		m_m_vtxadd_Kshort			);
			status = m_tuple->addItem("m_vtxadd_Kstar",		m_m_vtxadd_Kstar			);
			status = m_tuple->addItem("m_km_Kshort",		m_m_km_Kshort				);
			status = m_tuple->addItem("m_km_Kstar",			m_m_km_Kstar				);
			status = m_tuple->addItem("m_km_1resonance_Kshort",	m_m_km_1resonance_Kshort		);
			status = m_tuple->addItem("m_km_1resonance_Kstar",	m_m_km_1resonance_Kstar			);
			status = m_tuple->addItem("m_km_2resonance_Kshort",	m_m_km_2resonance_Kshort		);
			status = m_tuple->addItem("m_km_2resonance_Kstar",	m_m_km_2resonance_Kstar			);
			status = m_tuple->addItem("m_km_4p2resonance_Kshort",	m_m_km_4p2resonance_Kshort		);
			status = m_tuple->addItem("m_km_4p2resonance_Kstar",	m_m_km_4p2resonance_Kstar		);
			status = m_tuple->addItem("chisq_kmfit_1resonance",	m_chisq_kmfit_1resonance		);
			status = m_tuple->addItem("chisq_kmfit_2resonance",	m_chisq_kmfit_2resonance		);
			status = m_tuple->addItem("chisq_kmfit_4p2resonance",	m_chisq_kmfit_4p2resonance		);
		}
		else   
		{ 
			log << MSG::FATAL << "    Cannot book N-tuple:" << long(m_tuple) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	//
	//--------end of book--------
	//

	//--------------------------------------------------------------
	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode JpsiToPiKKshort::execute() 
{
	StatusCode sc = StatusCode::SUCCESS;

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	n_total_event++;

	SmartDataPtr<EvtRecTrackCol> evtRecTrackCol( eventSvc(), "/Event/EvtRec/EvtRecTrackCol" );
	if ( !evtRecTrackCol  ) {
		log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
		return StatusCode::FAILURE;

	}

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if (!eventHeader)
	{
		log << MSG::FATAL << "Could not find Event Header" << endreq;
		return StatusCode::FAILURE;
	}
	if(debug)cout<<__LINE__<<endl;
	log << MSG::DEBUG << __LINE__ << "" << endmsg;
	m_run = eventHeader->runNumber();
	m_rec = eventHeader->eventNumber();
	int runNo = m_run;
	int event = m_rec;
	int nE  = 0;
	int nnE  = 0;
	int gEa = 0;
	int gEb = 0;
	if(m_rec%10000==0)cout<<"Run   "<<m_run<<"     Event   "<<m_rec<<endl; 
	//MC information
// rm mctrue lorentz	
// 	HepLorentzVector  mctrue_pp01, mctrue_pm01, mctrue_gam1_Eta0, mctrue_gam2_Eta0,mctrue_gam1_Eta1,mctrue_gam2_Eta1,mctrue_Eta0, mctrue_Eta1;
	HepLorentzVector  mctruth_Pip_Kshort, mctruth_Pin_Kshort, mctruth_Pin_Kstar, mctruth_Kp_Kstar, mctruth_Kshort, mctruth_Kstar;

	if(debug)cout<<__LINE__<<endl;
	log << MSG::DEBUG << __LINE__ << "" << endmsg;
	if (eventHeader->runNumber()<0)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
				"/Event/MC/McParticleCol");
		int m_numParticle = 0;
		if (!mcParticleCol)
		{
			log << MSG::FATAL << "Could not retrieve McParticelCol" << endmsg;
			return StatusCode::FAILURE;
		}
		else
		{
			log << MSG::DEBUG << __LINE__ << " in mc event part" << endmsg;
			bool jpsiDecay = false;
			bool psipDecay = false;
			int rootIndex = -1;
			bool strange =false;
			Event::McParticleCol::iterator iter_mc;
			for(iter_mc=mcParticleCol->begin(); iter_mc != mcParticleCol->end();
					iter_mc++)
			{
				log << MSG::DEBUG << __LINE__ << " in mc event loop" << endmsg;
				if((*iter_mc)->primaryParticle()&&(*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()==11) { strange =true;}
				if ((*iter_mc)->primaryParticle()) {
					log << MSG::DEBUG << __LINE__ << endmsg;
					continue;
				}
				if (!(*iter_mc)->decayFromGenerator()) {
					log << MSG::DEBUG << __LINE__ << endmsg;
					continue;
				}
				if ((*iter_mc)->particleProperty()==m_idNo) {
					log << MSG::DEBUG << __LINE__ << endmsg;
					jpsiDecay = true;
					rootIndex = (*iter_mc)->trackIndex();
				}
				if (!jpsiDecay) {
					log << MSG::DEBUG << __LINE__ << endmsg;
					continue;
				}


				log << MSG::DEBUG << __LINE__ << endmsg;
				int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
				if(strange && ((*iter_mc)->mother()).particleProperty()!=m_idNo) mcidx--;//added
				int pdgid = (*iter_mc)->particleProperty();
				m_pdgid[m_numParticle] = pdgid;
				m_trkIdx[m_numParticle]=rootIndex;
				m_motheridx[m_numParticle] = mcidx;
				m_numParticle += 1;
				log << MSG::DEBUG << __LINE__ << endmsg;
				if(debug)cout<<__LINE__<<endl;
				///@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				double mcpx = (*iter_mc)->initialFourMomentum().px(); //最初px
				double mcpy = (*iter_mc)->initialFourMomentum().py(); //最初py
				double mcpz = (*iter_mc)->initialFourMomentum().pz(); //最初pz
				double mcen = (*iter_mc)->initialFourMomentum().e();  //最初能量
				HepLorentzVector p4mc;
				p4mc.setPx(mcpx);
				p4mc.setPy(mcpy);
				p4mc.setPz(mcpz);
				p4mc.setE(mcen);
				double mmm = p4mc.mag();
				double ppp = p4mc.rho();

				if(debug)cout<<__LINE__<<endl;
				log << MSG::DEBUG << __LINE__ << endmsg;

				int mctrue_id = (*iter_mc)->particleProperty();
				int mctrue_mothPID = ((*iter_mc)->mother()).particleProperty();
				int mctrue_moth_mothPID=((*iter_mc)->mother().mother()).particleProperty();


				log << MSG::DEBUG << __LINE__ << endmsg;
				if (mctrue_id == 313 && mctrue_mothPID == 443) {
					mctruth_Kstar		= p4mc;
					m_mctruth_P4_Kstar[0]	= mcpx;
					m_mctruth_P4_Kstar[1] 	= mcpy;
					m_mctruth_P4_Kstar[2] 	= mcpz;
					m_mctruth_P4_Kstar[3] 	= mcen;
					m_mctruth_P_Kstar	= ppp;
					m_mctruth_M_Kstar 	= mmm;
				}
				if (mctrue_id == -211 && mctrue_mothPID == 313) {
					mctruth_Pin_Kstar		= p4mc;
					m_mctruth_P4_Pin_Kstar[0]	= mcpx;
					m_mctruth_P4_Pin_Kstar[1] 	= mcpy;
					m_mctruth_P4_Pin_Kstar[2] 	= mcpz;
					m_mctruth_P4_Pin_Kstar[3] 	= mcen;
					m_mctruth_P_Pin_Kstar		= ppp;
					m_mctruth_M_Pin_Kstar		= mmm;
				}
				if (mctrue_id == 321 && mctrue_mothPID == 313) {
					mctruth_Kp_Kstar		= p4mc;
					m_mctruth_P4_Kp_Kstar[0]	= mcpx;
					m_mctruth_P4_Kp_Kstar[1] 	= mcpy;
					m_mctruth_P4_Kp_Kstar[2] 	= mcpz;
					m_mctruth_P4_Kp_Kstar[3] 	= mcen;
					m_mctruth_P_Kp_Kstar		= ppp;
					m_mctruth_M_Kp_Kstar 		= mmm;
				}
				if (mctrue_id == 310 && mctrue_mothPID == 443) {
					mctruth_Kshort			= p4mc;
					m_mctruth_P4_Kshort[0]		= mcpx;
					m_mctruth_P4_Kshort[1]		= mcpy;
					m_mctruth_P4_Kshort[2]		= mcpz;
					m_mctruth_P4_Kshort[3]		= mcen;
					m_mctruth_P_Kshort		= ppp;
					m_mctruth_M_Kshort 		= mmm;
				}
				if (mctrue_id == 211 && mctrue_mothPID == 310) {
					mctruth_Pip_Kshort		= p4mc;
					m_mctruth_P4_Pip_Kshort[0]	= mcpx;
					m_mctruth_P4_Pip_Kshort[1] 	= mcpy;
					m_mctruth_P4_Pip_Kshort[2] 	= mcpz;
					m_mctruth_P4_Pip_Kshort[3] 	= mcen;
					m_mctruth_P_Pip_Kshort		= ppp;
					m_mctruth_M_Pip_Kshort 		= mmm;
				}
				if (mctrue_id == -211 && mctrue_mothPID == 310) {
					mctruth_Pin_Kshort		= p4mc;
					m_mctruth_P4_Pin_Kshort[0]	= mcpx;
					m_mctruth_P4_Pin_Kshort[1] 	= mcpy;
					m_mctruth_P4_Pin_Kshort[2] 	= mcpz;
					m_mctruth_P4_Pin_Kshort[3] 	= mcen;
					m_mctruth_P_Pin_Kshort		= ppp;
					m_mctruth_M_Pin_Kshort 		= mmm;
				}

// rm two eta
/*
				//@@@@@@@@@@@@@@@@@  two eta
				if(debug)cout<<__LINE__<<endl;
				if(mctrue_id== 221 && mctrue_mothPID ==443&& nE==0  )   //  eta  A
				{
					mctrue_Eta0=p4mc;

					mctrue_eta0_P4[0]=  mcpx; mctrue_eta0_P4[1] = mcpy; mctrue_eta0_P4[2] = mcpz; mctrue_eta0_P4[3] = mcen; mctrue_Peta0 = ppp; mctrue_meta0 = mmm ;
					nE++;
				}
				else if(mctrue_id== 221 && mctrue_mothPID ==443  &&nE ==1 )   //  eta B
				{
					mctrue_Eta1=p4mc;

					mctrue_eta1_P4[0]=  mcpx; mctrue_eta1_P4[1] = mcpy; mctrue_eta1_P4[2] = mcpz; mctrue_eta1_P4[3] = mcen; mctrue_Peta1 = ppp; mctrue_meta1 = mmm ;
				}
				if(debug)cout<<__LINE__<<endl;


				//@@@@@@@@@@@@@@@@@ two gamma for eta A / B
				if(debug)cout<<__LINE__<<endl;

				if(mctrue_id== 22 && mctrue_mothPID ==221 &&gEa==0  && nnE==0 )   //  22  gamma   fro  eta A
				{
					mctrue_gam1_Eta0=p4mc;

					mctrue_gam0_eta0_P4[0]=  mcpx; mctrue_gam0_eta0_P4[1] = mcpy; mctrue_gam0_eta0_P4[2] = mcpz; mctrue_gam0_eta0_P4[3] = mcen; mctrue_Pgam0_eta0 = ppp; mctrue_mgam0_eta0 = mmm ;
					gEa++;
				}
				else if(mctrue_id== 22 && mctrue_mothPID ==221  &&gEa ==1 &&nnE==0 )   //  22  gamma   fro  eta A
				{
					mctrue_gam2_Eta0=p4mc;

					mctrue_gam1_eta0_P4[0]=  mcpx; mctrue_gam1_eta0_P4[1] = mcpy; mctrue_gam1_eta0_P4[2] = mcpz; mctrue_gam1_eta0_P4[3] = mcen; mctrue_Pgam1_eta0 = ppp; mctrue_mgam1_eta0 = mmm ;
					nnE++;
				}
				else if ( mctrue_id== 22 && mctrue_mothPID ==221 && gEb==0 &&nnE==1)   //  22  gamma   fro  eta B
				{
					mctrue_gam1_Eta1=p4mc;

					mctrue_gam0_eta1_P4[0]=  mcpx; mctrue_gam0_eta1_P4[1] = mcpy; mctrue_gam0_eta1_P4[2] = mcpz; mctrue_gam0_eta1_P4[3] = mcen; mctrue_Pgam0_eta1 = ppp; mctrue_mgam0_eta1 = mmm ;
					gEb++;
				}
				else if(mctrue_id== 22 && mctrue_mothPID ==221  &&gEb ==1 &&nnE==1 )   //  22  gamma   fro  eta B
				{
					mctrue_gam2_Eta1=p4mc;

					mctrue_gam1_eta1_P4[0]=  mcpx; mctrue_gam1_eta1_P4[1] = mcpy; mctrue_gam1_eta1_P4[2] = mcpz; mctrue_gam1_eta1_P4[3] = mcen; mctrue_Pgam1_eta1 = ppp; mctrue_mgam1_eta1 = mmm ;
				}
				if(debug)cout<<__LINE__<<endl;

				//@@@@@@@@@@@@@@@@@    p   pbar
				if(mctrue_id== 2212 && mctrue_mothPID ==443  )   // p+ 
				{                  
					mctrue_pp01=p4mc;

					mctrue_Ppp_P4[0]=  mcpx; mctrue_Ppp_P4[1] = mcpy; mctrue_Ppp_P4[2] = mcpz; mctrue_Ppp_P4[3] = mcen; mctrue_Ppp = ppp; mctrue_mpp = mmm ;
				}

				if(mctrue_id== -2212 && mctrue_mothPID ==443  )   // anti-p-
				{                  
					mctrue_pm01=p4mc;

					mctrue_Ppm_P4[0]=  mcpx; mctrue_Ppm_P4[1] = mcpy; mctrue_Ppm_P4[2] = mcpz; mctrue_Ppm_P4[3] = mcen; mctrue_Ppm = ppp; mctrue_mpm = mmm ;
				}
*/
				if(debug)cout<<__LINE__<<endl;
				log << MSG::DEBUG << __LINE__ << "" << endmsg;

                                HepLorentzVector ZXX(0.011*m_ecms,0.0,0.0,m_ecms);
                                Hep3Vector betaLab2 = ZXX.boostVector();
				
                                HepLorentzVector temp_Kshort     = mctruth_Kshort;
                                HepLorentzVector temp_Kstar      = mctruth_Kstar;
                                HepLorentzVector temp_Pip_Kshort = mctruth_Pip_Kshort;
                                HepLorentzVector temp_Pin_Kshort = mctruth_Pin_Kshort;
                                HepLorentzVector temp_Pin_Kstar  = mctruth_Pin_Kstar;
                                HepLorentzVector temp_Kp_Kstar   = mctruth_Kp_Kstar;
			
				m_mctruth_cos_Kshort	=(temp_Kshort    ).cosTheta();
				m_mctruth_cos_Kstar     =(temp_Kstar     ).cosTheta();
				m_mctruth_cos_Pip_Kshort=(temp_Pip_Kshort).cosTheta();
				m_mctruth_cos_Pin_Kshort=(temp_Pin_Kshort).cosTheta();
				m_mctruth_cos_Pin_Kstar =(temp_Pin_Kstar ).cosTheta();
				m_mctruth_cos_Kp_Kstar  =(temp_Kp_Kstar  ).cosTheta();


                                //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  0422 Add For MCtruth
                                
/* rm mctruth mc 
                                HepLorentzVector ZXX(0.011*m_ecms,0.0,0.0,m_ecms);
                                Hep3Vector betaLab2 = ZXX.boostVector();
                                HepLorentzVector temGam1_Eta0= mctrue_gam1_Eta0;
                                HepLorentzVector temGam2_Eta0= mctrue_gam2_Eta0;
                                HepLorentzVector temGam1_Eta1= mctrue_gam1_Eta1;
                                HepLorentzVector temGam2_Eta1= mctrue_gam2_Eta1;

                                HepLorentzVector tem_Eta0= mctrue_Eta0;
                                HepLorentzVector tem_Eta1= mctrue_Eta1;
                                HepLorentzVector tem_pp01= mctrue_pp01;
                                HepLorentzVector tem_pm01= mctrue_pm01;

                                m_mctruth_cos_gam1_Eta0=(temGam1_Eta0).cosTheta();
                                m_mctruth_cos_gam2_Eta0=(temGam2_Eta0).cosTheta();
                                m_mctruth_cos_gam1_Eta1=(temGam1_Eta1).cosTheta();
                                m_mctruth_cos_gam2_Eta1=(temGam2_Eta1).cosTheta();

                                m_mctruth_cos_Eta0=(tem_Eta0).cosTheta();
                                m_mctruth_cos_Eta1=(tem_Eta1).cosTheta();
                                m_mctruth_cos_pp=(tem_pp01).cosTheta();
                                m_mctruth_cos_pm=(tem_pm01).cosTheta();
*/



			}
//0422			m_tuple_angular->write();

		}
		m_idxmc = m_numParticle;
		if(debug)cout<<__LINE__<<endl;
		log << MSG::DEBUG << __LINE__ << "" << endmsg;

	}

	if(debug)cout<<__LINE__<<endl;
	log << MSG::DEBUG << __LINE__ << "" << endmsg;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),EventModel::EvtRec::EvtRecTrackCol);
	//	if(evtRecEvent->totalCharged()>10)return StatusCode::SUCCESS;
	//	if(evtRecEvent->totalNeutral()>20)return StatusCode::SUCCESS;
	m_ntotCh = evtRecEvent->totalCharged();
	m_ntotNe = evtRecEvent->totalNeutral();
	HepLorentzVector psip(0.011*m_ecms,0.0,0.0,m_ecms);

	Vint iGood,iChrgp,iChrgn;;
	iGood.clear();
	iChrgp.clear();
	iChrgn.clear();
	int nCharge = 0;


	Hep3Vector xorigin(0,0,0);
	HepSymMatrix xoriginEx(3, 0);
	IVertexDbSvc* vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid())
	{
		double* dbv = vtxsvc->PrimaryVertex();     //vertex[0] = vx; vertex[1]= vy; vertex[2] = vz;
		double* dbvsigma = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
		xoriginEx[0][0] = dbvsigma[0] * dbvsigma[0];
		xoriginEx[1][1] = dbvsigma[1] * dbvsigma[1];
		xoriginEx[2][2] = dbvsigma[2] * dbvsigma[2];
	} else {
		log << MSG::DEBUG << __LINE__ << " returned when set IP" << endmsg;
		return StatusCode::SUCCESS;
	}
	VertexParameter bs;
	bs.setVx(xorigin);
	bs.setEvx(xoriginEx);


	int cc=0;
	for(int i = 0;i < evtRecEvent->totalCharged(); i++)
	{
		log << MSG::DEBUG << __LINE__ << " in loop with index " << i << endmsg;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk =(*itTrk)->mdcTrack();
		if (!(*itTrk)->isMdcKalTrackValid()) continue;//MdcKalTrk
		double x0 = mdcTrk->x();
		double y0 = mdcTrk->y();
		double z0 = mdcTrk->z();
		double r0 = mdcTrk->r();
		double phi0 = mdcTrk->helix(1);
		double ptrk = mdcTrk->p();
		double chi2 = mdcTrk->chi2();
		double xv=xorigin.x();
		double yv=xorigin.y();
		double zv=xorigin.z();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		m_vx0[cc] = x0;
		m_vy0[cc] = y0;
		m_vz0[cc] = z0;
		m_vr0[cc] = Rxy;
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];       
		double Vctheta =cos(mdcTrk->theta());
		m_costheta[cc]=Vctheta;
		m_rvxy0[cc]=Rvxy0;
		m_rvz0[cc]=Rvz0;
		m_rvphi0[cc]=Rvphi0;
		cc++;
		if(fabs(Rvz0) >=m_vz0cut ) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		if(fabs(Vctheta)>0.93)continue;
		iGood.push_back(i);
		nCharge += mdcTrk->charge();		
		if(mdcTrk->charge()>0)
			iChrgp.push_back(i);
		else
			iChrgn.push_back(i);
		log << MSG::DEBUG << __LINE__ << " exit loop with index " << i << endmsg;
	}
	m_cc=cc;
	if(debug)cout<<__LINE__<<endl;
	log << MSG::DEBUG << __LINE__ << endmsg;
	int ntot=iGood.size();
	if((ntot != 4)||(nCharge != 0)) {
		log << MSG::DEBUG << __LINE__ << " failure with ntot " << ntot << " and nCharge " << nCharge << endmsg;
		return sc;
	}
	n_charge_cutted++;

	int nChp=iChrgp.size();
	int nChn=iChrgn.size();
	if((nChp != 2)||(nChn != 2)) {
		log << MSG::DEBUG << __LINE__ << " failure with nChp " << nChp << " and nChn " << nChn << endmsg;
		return sc;
	}
	n_charged_partical_cutted++;

	//------------------------------------------------------------------------------------------
	//emc
/*	double eraw[2],eoverp[2],ptrk[2],px_[2],py_[2],pz_[2],ptot[2];
	int goodshower[2];
	for(int i=0;i<2;i++){ goodshower[i]=0;eraw[i]=eoverp[i]=ptrk[i]=0.; px_[i]=py_[i]=pz_[i]=ptot[i]=0.; }
	for(int i_emc = 0; i_emc < ntot; i_emc++)
	{
		EvtRecTrackIterator  itTrk = evtRecTrackCol->begin() + iGood[i_emc];
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		ptrk[i_emc] = mdcTrk->p();
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		eraw[i_emc] = emcTrk->energy();
		eoverp[i_emc] = eraw[i_emc]/ptrk[i_emc];

		px_[i_emc] = mdcTrk->px();
		py_[i_emc] = mdcTrk->py();
		pz_[i_emc] = mdcTrk->pz();
		ptot[i_emc] = sqrt(px_[i_emc]*px_[i_emc] + py_[i_emc]*py_[i_emc]+ pz_[i_emc]*pz_[i_emc]);

		goodshower[i_emc]++;
	}

	//muc
	double depth[2],hits[2],layers[2],brlastlay[2],eclastlay[2];

	for(int i=0;i<ntot;i++){
		layers[i] = -2;
		EvtRecTrackIterator  itTrk = evtRecTrackCol->begin() + iGood[i];
		if((*itTrk)->isMucTrackValid()){
			RecMucTrack *evt_mu = (*itTrk)->mucTrack();
			depth[i] = (evt_mu)->depth();
			hits[i] = (evt_mu)->numHits();
			layers[i] = (evt_mu)->numLayers();
			brlastlay[i]= (evt_mu)->brLastLayer();
			eclastlay[i] = (evt_mu)->ecLastLayer();
		}
	}

	m_eraw1 = eraw[0];
	m_eraw2 = eraw[1];

	m_eoverp1 = eoverp[0];
	m_eoverp2 = eoverp[1];

	m_ptot1 = ptot[0];
	m_ptot2 = ptot[1];

	m_depth1 = depth[0];
	m_depth2 = depth[1];
*/
	//------------------------------------------------------------------------------------------

	Vint iGam;
	iGam.clear();
	int getTime=0;
	int nn = 0;
	double dang = 999.;
	for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++)
	{
		log << MSG::DEBUG << __LINE__ << " in loop with index i " << i << endmsg;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

		//double dang = 200.;
		for(int j = 0; j < evtRecEvent->totalCharged(); j++)
		{
			log << MSG::DEBUG << __LINE__ << "in loop with index j" << j << endmsg;
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			double angd = extpos.angle(emcpos);
			if(angd < dang)dang = angd;
			log << MSG::DEBUG << __LINE__ << " exit loop with index j" << j << endmsg;
		}

		double eraw = emcTrk->energy();
		int time= emcTrk->time();
		double theta = emcTrk->theta();
		double costheta=cos(theta);
		dang = dang * 180 / (CLHEP::pi);
		m_costhetaN[nn]=costheta;
		m_dang[nn] = dang;
		m_eraw[nn] = eraw;
		m_time[nn] =time;
		nn++;
		if(fabs(costheta)<0.8)
		{
			//m_energyThreshold = 0.025;
			m_energyThreshold = m_energyThreshold_b;
		}

		else if(fabs(costheta)< 0.92&&fabs(costheta)>0.86)
		{
			//m_energyThreshold = 0.05;
			m_energyThreshold = m_energyThreshold_e;
		}
		else
		{
			continue;
		}

		if(time>15||time<0)continue;
		if(eraw < m_energyThreshold) continue;
		iGam.push_back(i);
		log << MSG::DEBUG << __LINE__ << "exit loop with index i" << i << endmsg;
	}
	m_nn=nn;
	m_dang2 = dang;
	int nGam=iGam.size();
	//1009	if((nGam<4||(nGam>12))) return sc;
// rm gamma limit
	// if(( nGam<4) || ( nGam > 20 )) return sc;
	m_nGam = nGam;
	Ncut2++;


	if(debug)cout<<__LINE__<<endl;
	log << MSG::DEBUG << __LINE__ << "" << endmsg;
	Vp4 pGam;
	pGam.clear();
	for (int i = 0; i < nGam; i++)
	{
		log << MSG::DEBUG << __LINE__ << " loop with index " << i << endmsg;
		EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + iGam[i];
		RecEmcShower* emcTrk = (*itTrk) -> emcShower();

		double eraw = emcTrk -> energy();
		double phi  = emcTrk -> phi();
		double the  = emcTrk -> theta();

		HepLorentzVector ptrk;
		ptrk.setPx(eraw * sin(the) * cos(phi));
		ptrk.setPy(eraw * sin(the) * sin(phi));
		ptrk.setPz(eraw * cos(the));
		ptrk.setE(eraw);
		//ptrk = ptrk.boost(-boov);
		pGam.push_back(ptrk);
		log << MSG::DEBUG << __LINE__ << " exit loop with index" << i << endmsg;
	}
	if(debug)cout<<__LINE__<<endl;
	log << MSG::DEBUG << __LINE__ << "" << endmsg;

	Vint ikap;
	ikap.clear();
	Vp4 pkap;
	pkap.clear();

	//-----------------------Test fot PID ------------------------//


	double mpip = xmass[2], mpim = xmass[2], mpm = xmass[4], mpp = xmass[4];//pi+ pi- p- p质量
	double mkp = xmass[3]; // mass of kaon plus
	//  Vp4 ppip,ppm,ppin,ppp;//动量
	Vp4 ppm,ppin,ppp,ppip,pkp,pkn;//动量
	//  Vint ipip,ipm,ipin,ipp;//容器 存储事例信息？
	Vint ipm,ipin,ipp,ipip,ikp,ikn;//容器 存储事例信息？
	ppip.clear();
	ppm.clear();
	ppin.clear();
	ppp.clear();
	ipip.clear();
	ipm.clear();
	ipin.clear();
	ipp.clear();
	ppip.clear();
	pkp.clear();
	ikp.clear();
	pkn.clear();
	ikn.clear();
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nChp; i++)
	{
		log << MSG::DEBUG << __LINE__ << " in loop with index" << i << endmsg;
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgp[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon() |pid->onlyProton());    // seperater Pion/Kaon
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		double m_prob_pi = pid->probPion();
		double m_prob_K = pid->probKaon();
		double m_prob_p = pid->probProton();
		//if((pid->probPion()<pid->probKaon())||(pid->probPion()<pid->probProton())) continue;
		if(m_prob_p>m_prob_pi && m_prob_p>m_prob_K)//如果质子p的概率大于pi+,K+,则它就是质子
		{
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);//PID can set to e//这里是质子 proton 
			ipp.push_back(iChrgp[i]);


		}
		if(m_prob_pi>m_prob_K && m_prob_pi>m_prob_p) {
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
			RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
			ipip.push_back(iChrgp[i]);
		}
		if(m_prob_K>m_prob_pi && m_prob_K>m_prob_p) {
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
			RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
			ikp.push_back(iChrgp[i]);
		}
		log << MSG::DEBUG << __LINE__ << " exit loop with index " << i << endmsg;
	}
	for(int i = 0; i < nChn; i++)//chrn指  带负电径迹
	{
		log << MSG::DEBUG << __LINE__ << " in loop with index " << i << endmsg;
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgn[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon() |pid->onlyProton());    // seperater Pion/Kaon
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		double m_prob_p = pid->probProton();
		double m_prob_pi = pid->probPion();
		double m_prob_K = pid->probKaon();
		if(m_prob_p>m_prob_pi && m_prob_p>m_prob_K)//如果质子p的概率大于pi+,K+,则它就是质子
		{
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);//PID can set to e//这里是质子 proton 
			ipm.push_back(iChrgn[i]);

		}	
		if(m_prob_pi>m_prob_K && m_prob_pi>m_prob_p) {
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
			RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
			ipin.push_back(iChrgn[i]);
		}
		if(m_prob_K>m_prob_pi && m_prob_K>m_prob_p) {
			RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
			RecMdcKalTrack::setPidType(RecMdcKalTrack::kaon);
			ikn.push_back(iChrgn[i]);
		}
		log << MSG::DEBUG << __LINE__ << " exit loop with index " << i << endmsg;
	}
	int npp = ipp.size();
	int npm = ipm.size();
	int npip = ipip.size();
	int npin = ipin.size();
	int nkp  = ikp.size();
	int nkn  = ikn.size();

	
	
	m_npp=npp;//number of 质子
	m_npm=npm; //反质子数
	m_npip = npip;
        m_npin = npin;
        m_nkp  = nkp ;
        m_nkn  = nkn ;
// rm proton number limit
	log << MSG::DEBUG << __LINE__ << endmsg;
	// if(npp*npm !=1) return sc;//add  0721 
	if (npip*npin != 2 || (nkp+nkn) != 1) {
		log << MSG::DEBUG << __LINE__ << " failure with npip " << npip << " and npin " << npin << " and nkp " << nkp << " and nkn " << nkn << endmsg;
		return sc;
	}
	n_particle_number_cutted++;

	//-----------------PID END--------------------------------------------------------------------------------//    

	HepPoint3D vx(0.0, 0.0, 0.0);//初始化顶点
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx*bx;
	Evx[1][1] = by*by;
	Evx[2][2] = bz*bz;

	VertexParameter vxpar;
	vxpar.setVx(vx);
	vxpar.setEvx(Evx);

	log << MSG::DEBUG << __LINE__ << "" << endmsg;
	Vint temp_ipip, temp_ipin;
	temp_ipip.clear(); temp_ipin.clear();
	Vdouble temp_dchi, temp_dratio;
	temp_dchi.clear(); temp_dratio.clear();

	for (int i=0;i<ipip.size();i++) {
		log << MSG::DEBUG << __LINE__ << " in loop with index i " << i << endmsg;
		RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+ipip[i]))->mdcKalTrack();
		WTrackParameter wvPipTrk = WTrackParameter(xmass[2], pipTrk->getZHelixP(), pipTrk->getZErrorP());

		for (int j=0;j<ipin.size();j++) {
			log << MSG::DEBUG << __LINE__ << " in loop with index j" << j << endmsg;
			RecMdcKalTrack *pinTrk = (*(evtRecTrkCol->begin()+ipin[j]))->mdcKalTrack();
			WTrackParameter wvPinTrk = WTrackParameter(xmass[2], pinTrk->getZHelixP(), pinTrk->getZErrorP());
		
			VertexFit* vtxfit = VertexFit::instance();
			vtxfit->init();
			vtxfit->AddTrack(0, wvPipTrk);
			vtxfit->AddTrack(1, wvPinTrk);
			vtxfit->AddVertex(0, vxpar,0,1);
			if (vtxfit->Fit()) {
				vtxfit->BuildVirtualParticle(0);
				WTrackParameter wvKshortTrk = vtxfit->wVirtualTrack(0);
		
				SecondVertexFit* svtxfit = SecondVertexFit::instance();
				svtxfit->init();
				svtxfit->setPrimaryVertex(bs);
				svtxfit->AddTrack(0, wvKshortTrk);
				svtxfit->setVpar(vtxfit->vpar(0));
				if(svtxfit->Fit())
					temp_dratio.push_back(svtxfit->decayLength()/svtxfit->decayLengthError());
				else 
					temp_dratio.push_back(0);
			} else {	
				temp_dratio.push_back(0);
			}
			temp_ipip.push_back(ipip[i]);
			temp_ipin.push_back(ipin[j]);
			temp_dchi.push_back(vtxfit->chisq()); // When vtxfit->Fit() return false, the chisq would be 9999 in default. NO to use vtxfit->Fit(0), it will NO updata chisq then only 9999 you can get by chisq().
			log << MSG::DEBUG << __LINE__ << " exit loop with index j " << j << endmsg;
		}
		log << MSG::DEBUG << __LINE__ << " exit loop with index i "<< i << endmsg;
	}
	log << MSG::DEBUG << __LINE__ << "" << endmsg;
	
	int id_pip_kshort, id_pin_kshort, id_pin_kstar;
	if (ipin.size()!=2) {
		log << MSG::DEBUG << __LINE__ << " failure with Pi- number " << ipin.size() << endmsg;
		return sc;
	}
	if (temp_ipin.size()!=2) {
		log << MSG::DEBUG << __LINE__ << " failure with Pi+ number " << ipip.size() << endmsg;
		return sc;
	}
	// TODO pi0 eta0 highest E gamma
	// TODO K0 how to Ks
	// TODO compare kinematic fit with Kshort(Pi+ Pi-) and Pi(Kstar) K(Kstar)
	// TODO refactor ntuple about Kstar final states
	// TODO roofit
	// TODO topologic topoana.exe inMC_DpDm.card in /besfs5/groups/psipp/psippgroup/public/zhanghan/DpKsKlpi/04_all_newroot/topo
	n_pip_pin_vtxfit_cutted++;
	if (temp_dratio[0]>temp_dratio[1]) {
		id_pip_kshort = temp_ipip[0];
		id_pin_kshort = temp_ipin[0];
		id_pin_kstar  = temp_ipin[1];
		if (temp_dratio[0]<1)
			return sc;
		m_ratio_Kshort= temp_dratio[0];
		m_ratio_other = temp_dratio[1];
	} else {
		id_pip_kshort = temp_ipip[1];
		id_pin_kshort = temp_ipin[1];
		id_pin_kstar  = temp_ipin[0];
		if (temp_dratio[1]<1)
			return sc;
		m_ratio_Kshort= temp_dratio[1];
		m_ratio_other = temp_dratio[0];
	}

	log << MSG::DEBUG << __LINE__ << " id print: pip_kshort pin_kshort pin_kstar" << id_pip_kshort << " " << id_pin_kshort << " " << id_pin_kstar << endmsg;
	
	RecMdcKalTrack *pip_kshortTrk = (*(evtRecTrkCol->begin()+id_pip_kshort))->mdcKalTrack();
	RecMdcKalTrack *pin_kshortTrk = (*(evtRecTrkCol->begin()+id_pin_kshort))->mdcKalTrack();
	RecMdcKalTrack *pin_kstarTrk  = (*(evtRecTrkCol->begin()+id_pin_kstar))->mdcKalTrack();
	RecMdcKalTrack *kp_kstarTrk   = (*(evtRecTrkCol->begin()+ikp[0]))->mdcKalTrack();
	WTrackParameter wvpip_kshortTrk = WTrackParameter(xmass[2], pip_kshortTrk->getZHelixP(), pip_kshortTrk->getZErrorP());
	WTrackParameter wvpin_kshortTrk = WTrackParameter(xmass[2], pin_kshortTrk->getZHelixP(), pin_kshortTrk->getZErrorP());
	WTrackParameter wvpin_kstarTrk  = WTrackParameter(xmass[2], pin_kstarTrk->getZHelixP(),  pin_kstarTrk->getZErrorP());
	WTrackParameter wvkp_kstarTrk   = WTrackParameter(xmass[3], kp_kstarTrk->getZHelixP(),   kp_kstarTrk->getZErrorP());

	VertexFit* vtxfit_Kshort = VertexFit::instance();
	vtxfit_Kshort->init();
	vtxfit_Kshort->setChisqCut(5);
	vtxfit_Kshort->AddTrack(0,wvpip_kshortTrk);
	vtxfit_Kshort->AddTrack(1,wvpin_kshortTrk);
	vtxfit_Kshort->AddVertex(0,vxpar,0,1);
	log << MSG::DEBUG << __LINE__ << endmsg;
	if(!vtxfit_Kshort->Fit()) {
		log << MSG::DEBUG << __LINE__ << " vtxfit_kshort failure with chisq " << vtxfit_Kshort->chisq() << endmsg;
		return sc;
	}
	m_chisq_Kshort = vtxfit_Kshort->chisq();
	n_vtxfit_kshort_cutted++;
	log << MSG::DEBUG << __LINE__ << endmsg;
	WTrackParameter wpip_kshort = vtxfit_Kshort->wtrk(0);
	WTrackParameter wpin_kshort = vtxfit_Kshort->wtrk(1);
	vtxfit_Kshort->BuildVirtualParticle(0);
	WTrackParameter wkshort = vtxfit_Kshort->wVirtualTrack(0);
	
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	kmfit->init();
	kmfit->AddTrack(0,wvpip_kshortTrk);
	kmfit->AddTrack(1,wvpin_kshortTrk);
	kmfit->AddTrack(2,wvpin_kstarTrk);
	kmfit->AddTrack(3,wvkp_kstarTrk);
	kmfit->AddResonance(0,0.4876,0,1);
	log << MSG::DEBUG << __LINE__ << endmsg;
	if(!kmfit->Fit()) {
		log << MSG::DEBUG << __LINE__ << " kmfit failure" << endmsg;
		m_chisq_kmfit = 9999;
		// return sc;
	} else {
		n_kmfit_cutted++;
		m_chisq_kmfit = kmfit->chisq();
	}

	// for debug use
	KalmanKinematicFit *kmfit_1resonance = KalmanKinematicFit::instance();
	kmfit_1resonance->init();
	kmfit_1resonance->AddTrack(0,wvpip_kshortTrk);
	kmfit_1resonance->AddTrack(1,wvpin_kshortTrk);
	kmfit_1resonance->AddTrack(2,wvpin_kstarTrk);
	kmfit_1resonance->AddTrack(3,wvkp_kstarTrk);
	kmfit_1resonance->AddResonance(0,0.4876,0,1);
	log << MSG::DEBUG << __LINE__ << endmsg;
	if(!kmfit_1resonance->Fit()) {
		log << MSG::DEBUG << __LINE__ << " kmfit_1resonance failure" << endmsg;
		m_chisq_kmfit_1resonance = 9999;
		// return sc;
	} else {
		m_chisq_kmfit_1resonance = kmfit_1resonance->chisq();
	}

	// for debug use
	KalmanKinematicFit *kmfit_2resonance = KalmanKinematicFit::instance();
	kmfit_2resonance->init();
	kmfit_2resonance->AddTrack(0,wvpip_kshortTrk);
	kmfit_2resonance->AddTrack(1,wvpin_kshortTrk);
	kmfit_2resonance->AddTrack(2,wvpin_kstarTrk);
	kmfit_2resonance->AddTrack(3,wvkp_kstarTrk);
	kmfit_2resonance->AddResonance(0,0.4876,0,1);
	kmfit_2resonance->AddResonance(1,0.896,2,3);
	log << MSG::DEBUG << __LINE__ << endmsg;
	if(!kmfit_2resonance->Fit()) {
		log << MSG::DEBUG << __LINE__ << " kmfit_2resonance failure" << endmsg;
		m_chisq_kmfit_2resonance = 9999;
		// return sc;
	} else {
		m_chisq_kmfit_2resonance = kmfit_2resonance->chisq();
	}

	// for debug use
	KalmanKinematicFit *kmfit_4p2resonance = KalmanKinematicFit::instance();
	kmfit_4p2resonance->init();
	kmfit_4p2resonance->AddTrack(0,wvpip_kshortTrk);
	kmfit_4p2resonance->AddTrack(1,wvpin_kshortTrk);
	kmfit_4p2resonance->AddTrack(2,wvpin_kstarTrk);
	kmfit_4p2resonance->AddTrack(3,wvkp_kstarTrk);
	kmfit_4p2resonance->AddResonance(0,0.4876,0,1);
	kmfit_4p2resonance->AddResonance(1,0.896,2,3);
	kmfit_4p2resonance->AddFourMomentum(2,psip);
	log << MSG::DEBUG << __LINE__ << endmsg;
	if(!kmfit_4p2resonance->Fit()) {
		log << MSG::DEBUG << __LINE__ << " kmfit_4p2resonance failure" << endmsg;
		m_chisq_kmfit_4p2resonance = 9999;
		// return sc;
	} else {
		m_chisq_kmfit_4p2resonance = kmfit_4p2resonance->chisq();
	}

	// for other selecte
	VertexFit* vtxfit_other_Kshort = VertexFit::instance();
	vtxfit_other_Kshort->setChisqCut(9999);
	vtxfit_other_Kshort->init();
	vtxfit_other_Kshort->AddTrack(0,wvpip_kshortTrk);
	vtxfit_other_Kshort->AddTrack(1,wvpin_kstarTrk);
	vtxfit_other_Kshort->AddVertex(0,vxpar,0,1);
	if(vtxfit_other_Kshort->Fit())
		n_vtxfit_other_Kshort_ok++;
	vtxfit_other_Kshort->BuildVirtualParticle(0);
	WTrackParameter wpip_other_Kshort = vtxfit_other_Kshort->wtrk(0);
	WTrackParameter wpin_other_Kshort = vtxfit_other_Kshort->wtrk(1);
	WTrackParameter wother_Kshort = vtxfit_other_Kshort->wVirtualTrack(0);
	m_chisq_other_Kshort = vtxfit_other_Kshort->chisq();


	VertexFit* vtxfit_other_Kstar = VertexFit::instance();
	vtxfit_other_Kstar->setChisqCut(9999);
	vtxfit_other_Kstar->init();
	vtxfit_other_Kstar->AddTrack(0,wvkp_kstarTrk);
	vtxfit_other_Kstar->AddTrack(1,wvpin_kshortTrk);
	vtxfit_other_Kstar->AddVertex(0,vxpar,0,1);
	if(vtxfit_other_Kstar->Fit())
		n_vtxfit_other_Kstar_ok++;
	vtxfit_other_Kstar->BuildVirtualParticle(0);
	WTrackParameter wpip_other_Kstar = vtxfit_other_Kstar->wtrk(0);
	WTrackParameter wpin_other_Kstar = vtxfit_other_Kstar->wtrk(1);
	WTrackParameter wother_Kstar = vtxfit_other_Kstar->wVirtualTrack(0);
	m_chisq_other_Kstar = vtxfit_other_Kstar->chisq();
	

	log << MSG::DEBUG << __LINE__ << endmsg;
	// m_chi2 = kmfit->chisq();
	m_chi2 = 0;

	log << MSG::DEBUG << __LINE__ << " init lorentz vector from vwTrk";
	Hep3Vector betaLab = psip.boostVector();

	HepLorentzVector P_Pip_Kshort	= HepLorentzVector(wpip_kshort.w()[0],wpip_kshort.w()[1],wpip_kshort.w()[2],wpip_kshort.w()[3]).boost(-betaLab);
	HepLorentzVector P_Pin_Kshort	= HepLorentzVector(wpin_kshort.w()[0],wpin_kshort.w()[1],wpin_kshort.w()[2],wpin_kshort.w()[3]).boost(-betaLab);
	HepLorentzVector P_Pin_Kstar	= HepLorentzVector(wvpin_kstarTrk.w()[0],wvpin_kstarTrk.w()[1],wvpin_kstarTrk.w()[2],wvpin_kstarTrk.w()[3]).boost(-betaLab);
	HepLorentzVector P_Kp_Kstar	= HepLorentzVector(wvkp_kstarTrk.w()[0],wvkp_kstarTrk.w()[1],wvkp_kstarTrk.w()[2],wvkp_kstarTrk.w()[3]).boost(-betaLab);
	HepLorentzVector P_Kshort	= HepLorentzVector(wkshort.w()[0],wkshort.w()[1],wkshort.w()[2],wkshort.w()[3]).boost(-betaLab);
	HepLorentzVector P_Kstar	= P_Pin_Kstar + P_Kp_Kstar;

	HepLorentzVector P_Pip_other_Kshort	= HepLorentzVector(wpip_other_Kshort.w()[0],wpip_other_Kshort.w()[1],wpip_other_Kshort.w()[2],wpip_other_Kshort.w()[3]).boost(-betaLab);
	HepLorentzVector P_Pin_other_Kshort	= HepLorentzVector(wpin_other_Kshort.w()[0],wpin_other_Kshort.w()[1],wpin_other_Kshort.w()[2],wpin_other_Kshort.w()[3]).boost(-betaLab);
	HepLorentzVector P_other_Kshort	= HepLorentzVector(wother_Kshort.w()[0],wother_Kshort.w()[1],wother_Kshort.w()[2],wother_Kshort.w()[3]).boost(-betaLab);
	
	HepLorentzVector P_Kp_other_Kstar	= HepLorentzVector(wpip_other_Kstar.w()[0],wpip_other_Kstar.w()[1],wpip_other_Kstar.w()[2],wpip_other_Kstar.w()[3]).boost(-betaLab);
	HepLorentzVector P_Pin_other_Kstar	= HepLorentzVector(wpin_other_Kstar.w()[0],wpin_other_Kstar.w()[1],wpin_other_Kstar.w()[2],wpin_other_Kstar.w()[3]).boost(-betaLab);
	HepLorentzVector P_other_Kstar	= HepLorentzVector(wother_Kstar.w()[0],wother_Kstar.w()[1],wother_Kstar.w()[2],wother_Kstar.w()[3]).boost(-betaLab);


	// cut for mass
	if (fabs(P_Kstar.m()-0.892)>0.12) {
		log << MSG::DEBUG << __LINE__ << " return with Kstar mass is " << P_Kstar.m() << endmsg;
		return sc;
	}
	n_mass_cutted++;


	//HepLorentzVector P_Pip_Kshort = kmfit->pfit(0);
	//HepLorentzVector P_Pin_Kshort = kmfit->pfit(1);
	//HepLorentzVector P_Pin_Kstar = kmfit->pfit(2);
	//HepLorentzVector P_Kp_Kstar = kmfit->pfit(3);
	//HepLorentzVector P_Kshort = P_Pip_Kshort + P_Pin_Kshort;
	//HepLorentzVector P_Kstar  = P_Pin_Kstar  + P_Kp_Kstar;

	log << MSG::DEBUG << __LINE__ << " save data to ntuple" << endmsg;
	// HepLorentzVector ptotal=(kmfit->pfit(0)+kmfit->pfit(1)+kmfit->pfit(2)+kmfit->pfit(3));
	HepLorentzVector ptotal = (P_Pip_Kshort+P_Pin_Kshort+P_Pin_Kstar+P_Kp_Kstar);
	m_p_total[0] = ptotal.px();
	m_p_total[1] = ptotal.py();
	m_p_total[2] = ptotal.pz();
	m_p_total[3] = ptotal.e();

	
	m_p_Pip_Kshort[0] = P_Pip_Kshort.px();
	m_p_Pip_Kshort[1] = P_Pip_Kshort.py();
	m_p_Pip_Kshort[2] = P_Pip_Kshort.pz();
	m_p_Pip_Kshort[3] = P_Pip_Kshort.e();
	m_d3_Pip_Kshort=sqrt( P_Pip_Kshort.px()*P_Pip_Kshort.px()+P_Pip_Kshort.py()*P_Pip_Kshort.py()+P_Pip_Kshort.pz()*P_Pip_Kshort.pz());
	m_m_Pip_Kshort = P_Pip_Kshort.m();


	m_p_Pin_Kshort[0] = P_Pin_Kshort.px();
	m_p_Pin_Kshort[1] = P_Pin_Kshort.py();
	m_p_Pin_Kshort[2] = P_Pin_Kshort.pz();
	m_p_Pin_Kshort[3] = P_Pin_Kshort.e();
	m_d3_Pin_Kshort=sqrt( P_Pin_Kshort.px()*P_Pin_Kshort.px()+P_Pin_Kshort.py()*P_Pin_Kshort.py()+P_Pin_Kshort.pz()*P_Pin_Kshort.pz());
	m_m_Pin_Kshort = P_Pin_Kshort.m();

	
	m_p_Pin_Kstar[0] = P_Pin_Kstar.px();
	m_p_Pin_Kstar[1] = P_Pin_Kstar.py();
	m_p_Pin_Kstar[2] = P_Pin_Kstar.pz();
	m_p_Pin_Kstar[3] = P_Pin_Kstar.e();
	m_d3_Pin_Kstar=sqrt( P_Pin_Kstar.px()*P_Pin_Kstar.px()+P_Pin_Kstar.py()*P_Pin_Kstar.py()+P_Pin_Kstar.pz()*P_Pin_Kstar.pz());
	m_m_Pin_Kstar = P_Pin_Kstar.m();

	
	m_p_Kp_Kstar[0] = P_Kp_Kstar.px();
	m_p_Kp_Kstar[1] = P_Kp_Kstar.py();
	m_p_Kp_Kstar[2] = P_Kp_Kstar.pz();
	m_p_Kp_Kstar[3] = P_Kp_Kstar.e();
	m_d3_Kp_Kstar=sqrt( P_Kp_Kstar.px()*P_Kp_Kstar.px()+P_Kp_Kstar.py()*P_Kp_Kstar.py()+P_Kp_Kstar.pz()*P_Kp_Kstar.pz());
	m_m_Kp_Kstar = P_Kp_Kstar.m();


	m_p_Kshort[0] = P_Kshort.px();
	m_p_Kshort[1] = P_Kshort.py();
	m_p_Kshort[2] = P_Kshort.pz();
	m_p_Kshort[3] = P_Kshort.e();
	m_d3_Kshort=sqrt( P_Kshort.px()*P_Kshort.px()+P_Kshort.py()*P_Kshort.py()+P_Kshort.pz()*P_Kshort.pz());
	m_m_Kshort = P_Kshort.m();


	m_p_Kstar[0] = P_Kstar.px();
	m_p_Kstar[1] = P_Kstar.py();
	m_p_Kstar[2] = P_Kstar.pz();
	m_p_Kstar[3] = P_Kstar.e();
	m_d3_Kstar=sqrt( P_Kstar.px()*P_Kstar.px()+P_Kstar.py()*P_Kstar.py()+P_Kstar.pz()*P_Kstar.pz());
	m_m_Kstar = P_Kstar.m();


	m_p_other_Kshort[0] = P_other_Kshort.px();
	m_p_other_Kshort[1] = P_other_Kshort.py();
	m_p_other_Kshort[2] = P_other_Kshort.pz();
	m_p_other_Kshort[3] = P_other_Kshort.e();
	m_d3_other_Kshort=sqrt( P_other_Kshort.px()*P_other_Kshort.px()+P_other_Kshort.py()*P_other_Kshort.py()+P_other_Kshort.pz()*P_other_Kshort.pz());
	m_m_other_Kshort = P_other_Kshort.m();

	
	m_p_Pip_other_Kshort[0] = P_Pip_other_Kshort.px();
	m_p_Pip_other_Kshort[1] = P_Pip_other_Kshort.py();
	m_p_Pip_other_Kshort[2] = P_Pip_other_Kshort.pz();
	m_p_Pip_other_Kshort[3] = P_Pip_other_Kshort.e();
	m_d3_Pip_other_Kshort=sqrt( P_Pip_other_Kshort.px()*P_Pip_other_Kshort.px()+P_Pip_other_Kshort.py()*P_Pip_other_Kshort.py()+P_Pip_other_Kshort.pz()*P_Pip_other_Kshort.pz());
	m_m_Pip_other_Kshort = P_Pip_other_Kshort.m();
	

	m_p_Pin_other_Kshort[0] = P_Pin_other_Kshort.px();
	m_p_Pin_other_Kshort[1] = P_Pin_other_Kshort.py();
	m_p_Pin_other_Kshort[2] = P_Pin_other_Kshort.pz();
	m_p_Pin_other_Kshort[3] = P_Pin_other_Kshort.e();
	m_d3_Pin_other_Kshort=sqrt( P_Pin_other_Kshort.px()*P_Pin_other_Kshort.px()+P_Pin_other_Kshort.py()*P_Pin_other_Kshort.py()+P_Pin_other_Kshort.pz()*P_Pin_other_Kshort.pz());
	m_m_Pin_other_Kshort = P_Pin_other_Kshort.m();
	

	m_p_other_Kstar[0] = P_other_Kstar.px();
	m_p_other_Kstar[1] = P_other_Kstar.py();
	m_p_other_Kstar[2] = P_other_Kstar.pz();
	m_p_other_Kstar[3] = P_other_Kstar.e();
	m_d3_other_Kstar=sqrt( P_other_Kstar.px()*P_other_Kstar.px()+P_other_Kstar.py()*P_other_Kstar.py()+P_other_Kstar.pz()*P_other_Kstar.pz());
	m_m_other_Kstar = P_other_Kstar.m();

	
	m_p_Kp_other_Kstar[0] = P_Kp_other_Kstar.px();
	m_p_Kp_other_Kstar[1] = P_Kp_other_Kstar.py();
	m_p_Kp_other_Kstar[2] = P_Kp_other_Kstar.pz();
	m_p_Kp_other_Kstar[3] = P_Kp_other_Kstar.e();
	m_d3_Kp_other_Kstar=sqrt( P_Kp_other_Kstar.px()*P_Kp_other_Kstar.px()+P_Kp_other_Kstar.py()*P_Kp_other_Kstar.py()+P_Kp_other_Kstar.pz()*P_Kp_other_Kstar.pz());
	m_m_Kp_other_Kstar = P_Kp_other_Kstar.m();
	

	m_p_Pin_other_Kstar[0] = P_Pin_other_Kstar.px();
	m_p_Pin_other_Kstar[1] = P_Pin_other_Kstar.py();
	m_p_Pin_other_Kstar[2] = P_Pin_other_Kstar.pz();
	m_p_Pin_other_Kstar[3] = P_Pin_other_Kstar.e();
	m_d3_Pin_other_Kstar=sqrt( P_Pin_other_Kstar.px()*P_Pin_other_Kstar.px()+P_Pin_other_Kstar.py()*P_Pin_other_Kstar.py()+P_Pin_other_Kstar.pz()*P_Pin_other_Kstar.pz());
	m_m_Pin_other_Kstar = P_Pin_other_Kstar.m();
	

	m_bkg_m_Pip_Kshort	=(psip-P_Pip_Kshort		).m();
	m_bkg_m_Pin_Kshort	=(psip-P_Pin_Kshort		).m();
	m_bkg_m_Pin_Kstar	=(psip-P_Pin_Kstar		).m();
	m_bkg_m_Kp_Kstar	=(psip-P_Kp_Kstar		).m();
	m_bkg_m_Kshort		=(psip-P_Kshort			).m();
	m_bkg_m_Kstar		=(psip-P_Kstar			).m();
	m_bkg_m_Pin_Pin		=(psip-P_Pin_Kshort-P_Pin_Kstar	).m();
	m_bkg_m_Pip_Kp		=(psip-P_Pip_Kshort-P_Kp_Kstar	).m();
	m_bkg_m_Pin_Kp		=(psip-P_Pin_Kshort-P_Kp_Kstar	).m();
	m_bkg_m_Pip_Pin		=(psip-P_Pip_Kshort-P_Pin_Kstar	).m();
	m_bkg_m_other_Kshort	=(psip-P_other_Kshort		).m();
	m_bkg_m_Pip_other_Kshort=(psip-P_Pip_other_Kshort	).m();
	m_bkg_m_Pin_other_Kshort=(psip-P_Pin_other_Kshort	).m();
	m_bkg_m_other_Kstar	=(psip-P_other_Kstar		).m();
	m_bkg_m_Kp_other_Kstar	=(psip-P_Kp_other_Kstar		).m();
	m_bkg_m_Pin_other_Kstar	=(psip-P_Pin_other_Kstar	).m();
	log << MSG::DEBUG << __LINE__ << " save ok" << endmsg;


	log << MSG::DEBUG << __LINE__ << " save polarization angle to ntuple" << endmsg;

	m_cos_Pip_Kshort	=(P_Pip_Kshort	).cosTheta();
	m_cos_Pin_Kshort	=(P_Pin_Kshort	).cosTheta();
	m_cos_Pin_Kstar		=(P_Pin_Kstar	).cosTheta();
	m_cos_Kp_Kstar		=(P_Kp_Kstar	).cosTheta();
	m_cos_Kshort		=(P_Kshort	).cosTheta();
	m_cos_Kstar		=(P_Kstar	).cosTheta();
	m_cos_Pip_other_Kshort	=(P_Pip_other_Kshort).cosTheta();
	m_cos_Pin_other_Kshort	=(P_Pin_other_Kshort).cosTheta();
	m_cos_other_Kshort	=(P_other_Kshort).cosTheta();
	m_cos_Kp_other_Kstar	=(P_Kp_other_Kstar).cosTheta();
	m_cos_Pin_other_Kshort	=(P_Pin_other_Kshort).cosTheta();
	m_cos_other_Kshort	=(P_other_Kshort).cosTheta();


	// for debug use
	m_m_orig_Kshort			= (wvpip_kshortTrk.p()		+ wvpin_kshortTrk.p()).m();
	m_m_orig_Kstar			= (wvpin_kstarTrk.p()		+ wvkp_kstarTrk.p()).m();
	m_m_vtxadd_Kshort		= (wpip_kshort.p()		+ wpin_kshort.p()).m();
	m_m_vtxadd_Kstar		= (wvpin_kstarTrk.p()		+ wvkp_kstarTrk.p()).m(); // equal to m_m_orig_Kstar
	m_m_km_Kshort			= (kmfit->pfit(0)		+ kmfit->pfit(1)).m();
	m_m_km_Kstar			= (kmfit->pfit(2)		+ kmfit->pfit(3)).m();
	m_m_km_1resonance_Kshort	= (kmfit_1resonance->pfit(0)	+ kmfit_1resonance->pfit(1)).m();
	m_m_km_1resonance_Kstar		= (kmfit_1resonance->pfit(2)	+ kmfit_1resonance->pfit(3)).m();
	m_m_km_2resonance_Kshort	= (kmfit_2resonance->pfit(0)	+ kmfit_2resonance->pfit(1)).m();
	m_m_km_2resonance_Kstar		= (kmfit_2resonance->pfit(2)	+ kmfit_2resonance->pfit(3)).m();
	m_m_km_4p2resonance_Kshort	= (kmfit_4p2resonance->pfit(0)	+ kmfit_4p2resonance->pfit(1)).m();
	m_m_km_4p2resonance_Kstar	= (kmfit_4p2resonance->pfit(2)	+ kmfit_4p2resonance->pfit(3)).m();
	log << MSG::DEBUG << __LINE__ << " save ok" << endmsg;


	log << MSG::DEBUG << __LINE__ << " write ntuple" << endmsg;
	m_tuple->write();
	m_tuple_angular->write();
	log << MSG::DEBUG << __LINE__ << " write ok" << endmsg;
	log << MSG::INFO  << "excute() done" << endmsg;
	return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode JpsiToPiKKshort::finalize() 
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	log << MSG::INFO << "total event number:		" << n_total_event		<< endmsg;
	log << MSG::INFO << "ntot==4, nCharge==0:		" << n_charge_cutted		<< endmsg;
	log << MSG::INFO << "nChp==2, nChn==2:			" << n_charged_partical_cutted 	<< endmsg;
	log << MSG::INFO << "npip*npin==2&&(nkp+nkn)==1		" << n_particle_number_cutted	<< endmsg;
	log << MSG::INFO << "npin==2&&vtexfitted_npin==2	" << n_pip_pin_vtxfit_cutted	<< endmsg;
	log << MSG::INFO << "vertex fit for kshort success	" << n_vtxfit_kshort_cutted	<< endmsg;
	log << MSG::INFO << "vertex fit for kstar success	" << n_vtxfit_kstar_cutted	<< endmsg;
	log << MSG::INFO << "vertex fit for other Kshort success" << n_vtxfit_other_Kshort_ok	<< endmsg;
	log << MSG::INFO << "vertex fit for other Kstar success	" << n_vtxfit_other_Kstar_ok	<< endmsg;
	log << MSG::INFO << "kinematic fit success		" << n_kmfit_cutted		<< endmsg;
	log << MSG::INFO << "mass check succes			" << n_mass_cutted		<< endmsg;
	//log <<"pi0 list:			" << Ncut4  << endmsg;
	//log <<"VertexFit:			" << Ncut5  << endmsg;
	//log <<"6C:				" << Ncut6  << endmsg;
	//log <<"6C:				" << Ncut7  << endmsg;
	return StatusCode::SUCCESS;
}
// *****************************************************************
// ** A macro to create correlated Gaussian-distributed variables **
// *****************************************************************
