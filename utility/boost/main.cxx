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
const size_t      kNProcess       = 10;
const double      kCmsEnergy      = 3.097;
const std::string kTreeName       = "chi2omg";
const std::string kFriendName     = "chi2omg";
const std::string kP4Gamma        = "p_Gamma_Kshort";
const std::string kP4Pip          = "p_Pip_Kshort";
const std::string kP4Pin          = "p_Pin_Kshort";
const std::string kTag            = "is_Gamma_PiPi";
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

    // those would be read from tree for caculating
    double raw_gamma[4];
    double raw_pip[4];
    double raw_pin[4];
    int    is_gamma_pipi;
    int    num_gamma;
    double dang[128];
    tree->SetBranchAddress(kP4Gamma.c_str(), &raw_gamma);
    tree->SetBranchAddress(kP4Pip.c_str(),   &raw_pip);
    tree->SetBranchAddress(kP4Pin.c_str(),   &raw_pin);
    tree->SetBranchAddress(kTag.c_str(),     &is_gamma_pipi);
    tree->SetBranchAddress(kNumGamma.c_str(),&num_gamma);
    tree->SetBranchAddress(kDAng.c_str(),    dang);

    // caculating and writing friend
    for(auto i : ROOT::TSeqI(tree->GetEntries())) {
        tree->GetEntry(i);
        if (is_gamma_pipi == 1) {
            LorentzVector p4_gamma;
            LorentzVector p4_pip;
            LorentzVector p4_pin;
            p4_gamma.SetCoordinates(raw_gamma);
            p4_pip  .SetCoordinates(raw_pip);
            p4_pin  .SetCoordinates(raw_pin);
            auto kshort_frame = Boost(
                (p4_pip+p4_pin+p4_gamma).BoostToCM() );
            auto p4_gamma_frameKshort = kshort_frame (p4_gamma);
            auto p4_pip_frameKshort   = kshort_frame (p4_pip);
            auto p4_pin_frameKshort   = kshort_frame (p4_pin);
            d3_Gamma_frameKshort        = p4_gamma_frameKshort.P();
            e_Gamma_frameKshort         = p4_gamma_frameKshort.E();
            ang_Gamma_Pi_frameKshort[0] = p4_gamma_frameKshort.Vect().Unit().Dot(
                p4_pip_frameKshort.Vect().Unit());
            ang_Gamma_Pi_frameKshort[1] = p4_gamma_frameKshort.Vect().Unit().Dot(
                p4_pin_frameKshort.Vect().Unit());
            auto lab_frame = Boost( 
                LorentzVector(0.011*kCmsEnergy, 0, 0, kCmsEnergy).BoostToCM()
                ).Inverse();
            auto p4_gamma_frameLab = lab_frame (p4_gamma);
            auto p4_pip_frameLab   = lab_frame (p4_pip);
            auto p4_pin_frameLab   = lab_frame (p4_pin);
            d3_Gamma_frameLab        = p4_gamma_frameLab.P();
            e_Gamma_frameLab         = p4_gamma_frameLab.E();
            ang_Gamma_Pi_frameLab[0] = p4_gamma_frameLab.Vect().Unit().Dot(
                p4_pip_frameLab.Vect().Unit());
            ang_Gamma_Pi_frameLab[1] = p4_gamma_frameLab.Vect().Unit().Dot(
                p4_pin_frameLab.Vect().Unit());
            auto pipi_frame = Boost( (p4_pip+p4_pin).BoostToCM() );
            auto p4_gamma_framePiPi = pipi_frame (p4_gamma);
            auto p4_pip_framePiPi   = pipi_frame (p4_pip);
            auto p4_pin_framePiPi   = pipi_frame (p4_pin);
            d3_Gamma_framePiPi        = p4_gamma_framePiPi.P();
            e_Gamma_framePiPi         = p4_gamma_framePiPi.E();
            ang_Gamma_Pi_framePiPi[0] = p4_gamma_framePiPi.Vect().Unit().Dot(
                p4_pip_framePiPi.Vect().Unit());
            ang_Gamma_Pi_framePiPi[1] = p4_gamma_framePiPi.Vect().Unit().Dot(
                p4_pin_framePiPi.Vect().Unit());
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
        if (num_gamma <= 0) {ang_Gamma_track=-99;}
        else {
            ang_Gamma_track=3600;
            for (int i=0; i<num_gamma; i++) {
                if (ang_Gamma_track>num_gamma) {
                    ang_Gamma_track=num_gamma;
                }
            }
        }
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
    write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/generator/IB/run/mc_job/selected.root");
    // process_file("../../inmc.list");
    // process_file("../../data.list");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/topo/selected.root");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/data/selected.root");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma_only_PiPi/run/topo/selected.root");
    // write_into("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma_only_PiPiGamma/run/topo/selected.root");
    return 0;
}
