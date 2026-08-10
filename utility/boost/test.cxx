#include <Math/GenVector/Boost.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/Vector4Dfwd.h>
#include <ROOT/TProcessExecutor.hxx>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/TSeq.hxx>
#include <fstream>
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector> LorentzVector;
using ROOT::Math::Boost;
#include "../read/read.h"
#include "func.h"
const size_t      kNProcess       = 10;
const double      kCmsEnergy      = 3.097;
const std::string kTreeName       = "mctruth";
const std::string kFriendName     = "mctruth";
const std::string kP4Gamma        = "mctruth_P4_Gamma_Kshort";
const std::string kP4Pip          = "mctruth_P4_Pip_Kshort";
const std::string kP4Pin          = "mctruth_P4_Pin_Kshort";
const std::string kP4Pi           = "mctruth_P4_Pi_Kstar";
const std::string kP4K            = "mctruth_P4_K_Kstar";
const std::string kTag            = "mctruth_is_Gamma_PiPi";
const std::string kNumGamma       = "nn";
const std::string kDAng           = "dang";
const std::string kSelection      = "chisq_kmfit<100&&ratio_Kshort>2&&m_Kshort>0.4&&m_Kshort<0.6";
const std::string kOutD3GammaFrameKshort     = "d3_Gamma_frameKshort";
const std::string kOutEGammaFrameKshort      = "e_Gamma_frameKshort";
const std::string kOutAngGammaFrameKshort    = "ang_Gamma_Pi_frameKshort";
const std::string kOutD3GammaFrameLab        = "d3_Gamma_frameLab";
const std::string kOutEGammaFrameLab         = "e_Gamma_frameLab";
const std::string kOutAngGammaFrameLab       = "ang_Gamma_Pi_frameLab";
const std::string kOutD3GammaFramePiPi       = "d3_Gamma_framePiPi";
const std::string kOutEGammaFramePiPi        = "e_Gamma_framePiPi";
const std::string kOutAngGammaFramePiPi      = "ang_Gamma_Pi_framePiPi";
const std::string kOutAngGammaTrack          = "ang_Gamma_track";
const std::string kOutMRecoilPip             = "m_recoil_pip";
const std::string kOutMRecoilPipPim          = "m_recoil_pip_pim";
const std::string kOutMRecoilPi              = "m_recoil_pi";
const std::string kOutMRecoilK               = "m_recoil_k";
std::string out_filename(const std::string& input)
{
    return std::string(input+".friend");
}
bool write_into(const std::string &filename)
{
    auto tree = open_root(filename, kTreeName, "");
    TFile file_friend(out_filename(filename).c_str(),"recreate");
    auto tree_friend = new TTree(kTreeName.c_str(), kTreeName.c_str());

    // those would be write in to friend tree
    double d3_Gamma_frameKshort = 0;
    double e_Gamma_frameKshort  = 0;
    double ang_Gamma_Pi_frameKshort[2] = {0,0};
    double d3_Gamma_frameLab = 0;
    double e_Gamma_frameLab  = 0;
    double ang_Gamma_Pi_frameLab[2] = {0,0};
    double d3_Gamma_framePiPi = 0;
    double e_Gamma_framePiPi  = 0;
    double ang_Gamma_Pi_framePiPi[2] = {0,0};
    double ang_Gamma_track  = 0;
    double m_recoil_pip     = 0;
    double m_recoil_pip_pim = 0;
    double m_recoil_pi      = 0;
    double m_recoil_k       = 0;
    tree_friend->Branch(kOutD3GammaFrameKshort.c_str(),  &d3_Gamma_frameKshort,
        (kOutD3GammaFrameKshort  +"/D"   ).c_str());
    tree_friend->Branch(kOutEGammaFrameKshort.c_str(),   &e_Gamma_frameKshort,
        (kOutEGammaFrameKshort   +"/D"   ).c_str());
    tree_friend->Branch(kOutAngGammaFrameKshort.c_str(), &ang_Gamma_Pi_frameKshort,
        (kOutAngGammaFrameKshort +"[2]/D").c_str());
    tree_friend->Branch(kOutD3GammaFrameLab.c_str(),  &d3_Gamma_frameLab,
        (kOutD3GammaFrameLab  +"/D"   ).c_str());
    tree_friend->Branch(kOutEGammaFrameLab.c_str(),   &e_Gamma_frameLab,
        (kOutEGammaFrameLab   +"/D"   ).c_str());
    tree_friend->Branch(kOutAngGammaFrameLab.c_str(), &ang_Gamma_Pi_frameLab,
        (kOutAngGammaFrameLab +"[2]/D").c_str());
    tree_friend->Branch(kOutD3GammaFramePiPi.c_str(),  &d3_Gamma_framePiPi,
        (kOutD3GammaFramePiPi  +"/D"   ).c_str());
    tree_friend->Branch(kOutEGammaFramePiPi.c_str(),   &e_Gamma_framePiPi,
        (kOutEGammaFramePiPi   +"/D"   ).c_str());
    tree_friend->Branch(kOutAngGammaFramePiPi.c_str(), &ang_Gamma_Pi_framePiPi,
        (kOutAngGammaFramePiPi +"[2]/D").c_str());
    tree_friend->Branch(kOutAngGammaTrack.c_str(), &ang_Gamma_track,
        (kOutAngGammaTrack +"/D").c_str());
    tree_friend->Branch(kOutMRecoilPip.c_str(), &m_recoil_pip,
        (kOutMRecoilPip +"/D").c_str());
    tree_friend->Branch(kOutMRecoilPipPim.c_str(), &m_recoil_pip_pim,
        (kOutMRecoilPipPim +"/D").c_str());
    tree_friend->Branch(kOutMRecoilPi.c_str(), &m_recoil_pi,
        (kOutMRecoilPi +"/D").c_str());
    tree_friend->Branch(kOutMRecoilK.c_str(), &m_recoil_k,
        (kOutMRecoilK +"/D").c_str());

    // those would be read from tree for caculating
    double raw_gamma[4];
    double raw_pip[4];
    double raw_pin[4];
    double raw_pi[4]; // for Kstar decay
    double raw_k[4];
    double raw_total[4] = {0.011*kCmsEnergy, 0, 0, kCmsEnergy};
    int    is_gamma_pipi;
    int    num_gamma;
    double dang[128];
    tree->SetBranchAddress(kP4Gamma.c_str(), &raw_gamma);
    tree->SetBranchAddress(kP4Pip.c_str(),   &raw_pip);
    tree->SetBranchAddress(kP4Pin.c_str(),   &raw_pin);
    tree->SetBranchAddress(kP4Pi.c_str(),    &raw_pi);
    tree->SetBranchAddress(kP4K.c_str(),     &raw_k);
    tree->SetBranchAddress(kTag.c_str(),     &is_gamma_pipi);
    // tree->SetBranchAddress(kNumGamma.c_str(),&num_gamma);
    // tree->SetBranchAddress(kDAng.c_str(),    &dang);

    // caculating and writing friend
    for(auto i : ROOT::TSeqI(tree->GetEntries())) {
        tree->GetEntry(i);

        // boost
        if (is_gamma_pipi == 1) {
            GammaData in_pipigamma_frame = BoostGammaIntoPiPiGammaFrame(raw_gamma, raw_pip, raw_pin);
            d3_Gamma_frameKshort         = in_pipigamma_frame.d3_gamma;
            e_Gamma_frameKshort          = in_pipigamma_frame.e_gamma;
            ang_Gamma_Pi_frameKshort[0]  = in_pipigamma_frame.ang_gamma_pip;
            ang_Gamma_Pi_frameKshort[1]  = in_pipigamma_frame.ang_gamma_pin;
            GammaData in_lab_frame       = BoostGammaIntoLabFrame(raw_gamma, raw_pip, raw_pin, raw_total);
            d3_Gamma_frameLab            = in_lab_frame.d3_gamma;
            e_Gamma_frameLab             = in_lab_frame.e_gamma;
            ang_Gamma_Pi_frameLab[0]     = in_lab_frame.ang_gamma_pip;
            ang_Gamma_Pi_frameLab[1]     = in_lab_frame.ang_gamma_pin;
            GammaData in_pipi_frame      = BoostGammaIntoPiPiFrame(raw_gamma, raw_pip, raw_pin);
            d3_Gamma_framePiPi           = in_pipi_frame.d3_gamma;
            e_Gamma_framePiPi            = in_pipi_frame.e_gamma;
            ang_Gamma_Pi_framePiPi[0]    = in_pipi_frame.ang_gamma_pip;
            ang_Gamma_Pi_framePiPi[1]    = in_pipi_frame.ang_gamma_pin;
        } else {
            d3_Gamma_frameKshort = 0;
            e_Gamma_frameKshort  = 0;
            ang_Gamma_Pi_frameKshort[0] = 0;
            ang_Gamma_Pi_frameKshort[1] = 0;
            d3_Gamma_frameLab = 0;
            e_Gamma_frameLab  = 0;
            ang_Gamma_Pi_frameLab[0] = 0;
            ang_Gamma_Pi_frameLab[1] = 0;
            d3_Gamma_framePiPi = 0;
            e_Gamma_framePiPi  = 0;
            ang_Gamma_Pi_framePiPi[0] = 0;
            ang_Gamma_Pi_framePiPi[1] = 0;
        }

        // select minimal dang as AngGammaTrack
        // ang_Gamma_track = MinDAng(num_gamma, dang);
        
        // calculate recoil
        RecoilData recoil_pip = CalRecoilPip(
            raw_pip,raw_pin,raw_gamma,raw_pi,raw_k,raw_total
        );
        RecoilData recoil_pip_pim = CalRecoilPipPim(
            raw_pip,raw_pin,raw_gamma,raw_pi,raw_k,raw_total
        );
        RecoilData recoil_pi = CalRecoilPi(
            raw_pip,raw_pin,raw_gamma,raw_pi,raw_k,raw_total
        );
        RecoilData recoil_k = CalRecoilK(
            raw_pip,raw_pin,raw_gamma,raw_pi,raw_k,raw_total
        );
        m_recoil_pip     = recoil_pip.m;
        m_recoil_pip_pim = recoil_pip_pim.m;
        m_recoil_pi      = recoil_pi.m;
        m_recoil_k       = recoil_k.m;

        tree_friend->Fill();
    }
    file_friend.Write();
    file_friend.Close();
    return true;
}
bool process_file(const std::string &listname)
{
    std::vector<std::string> lines;
    std::ifstream infile(listname);
    std::string line;
    while(std::getline(infile,line)) {
        lines.push_back(line);
    }
    auto num_line = lines.size();
    size_t num_work = 0;
    size_t num_done = 0;
    auto workItem = [num_line, &num_work, &num_done](std::string line) -> int {
        std::clog << "[INFO] " << ++num_work << "/" << num_line << " writing friend for " << line << std::endl;
        write_into(line);
        std::clog << "[INFO] " << ++num_done << "/" << num_line << " writen friend " << out_filename(line) << std::endl;
        return 0;
    };
    ROOT::TProcessExecutor pool(kNProcess);
    pool.Map(workItem, lines);
    return true;
}
int main()
{
    // manual mkdir first
    process_file("boost.list");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/generator/IB/run/mc_job/selected.root");
    // process_file("../../inmc.list");
    // process_file("../../data.list");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/topo/selected.root");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/data/selected.root");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma_only_PiPi/run/topo/selected.root");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma_only_PiPiGamma/run/topo/selected.root");
    return 0;
}
