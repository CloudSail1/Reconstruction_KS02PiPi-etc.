#include "../read/read.h"
#include <TFile.h>
#include <memory>

const std::string kTreeName   = "mctruth";
const std::string kFriendName = "";

bool select(const std::string &file, const std::string &exp, const std::string &outname)
{
    auto old_chain = open_root(file, kTreeName, kFriendName);
    TFile new_file(outname.c_str(), "recreate");
    auto new_chain = old_chain->CopyTree(exp.c_str());
    INFO("Writing: "+outname);
    new_file.Write();
    INFO("Writed.");
    return true;
}
int main()
{
    select("select.list",
            "1",
        "./selected.root");
    // select("../../inmc.list",
        // "chisq_kmfit<=100&&ratio_Kshort>=2&&m_Kshort>=0.47&&m_Kshort<=0.51",
    // "/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/topo/selected.root");
    // select("../../data.list",
    //     "chisq_kmfit<=100&&ratio_Kshort>=2&&m_Kshort>=0.47&&m_Kshort<=0.51",
    // "/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/data/selected.root");
    // select("../../exmc_pipi.list",
    //     "chisq_kmfit<=100&&ratio_Kshort>=2&&m_Kshort>=0.47&&m_Kshort<=0.51",
    // "/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma_only_PiPi/run/topo/selected.root");
    // select("../../exmc_pipigamma.list",
    //     "chisq_kmfit<=100&&ratio_Kshort>=2&&m_Kshort>=0.47&&m_Kshort<=0.51",
    // "/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma_only_PiPiGamma/run/topo/selected.root");
    return 0;
}
