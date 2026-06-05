#include "draw.h"
#include <fstream>
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

const size_t      kNProcessDefault = 20;
const size_t      kNBin            = 100;
const std::string kNoDraw          = "0";
const std::string kDrawAll         = "1";
const std::string kTreeName        = "chi2omg";
const std::string kJsonFile        = "../../draw/params.json";
const bool        kNormalize       = false;

DrawParams json_to_draw_params(const json &j_defaults, const json &j_item) {
    DrawParams p;
    // Helper: get a field from item, fall back to defaults
    auto get_str = [&](const std::string &key) -> std::string {
        if (j_item.contains(key))
            return j_item[key].get<std::string>();
        if (j_defaults.contains(key))
            return j_defaults[key].get<std::string>();
        return "";
    };
    auto get_num = [&](const std::string &key) -> double {
        if (j_item.contains(key))
            return j_item[key].get<double>();
        if (j_defaults.contains(key))
            return j_defaults[key].get<double>();
        return 0.0;
    };

    p.lo        = get_num("lo");
    p.up        = get_num("up");
    p.target    = get_str("target");
    p.title     = get_str("title");
    p.outname   = get_str("outname");
    p.inmc_list = get_str("inmc_list");
    p.data_list = get_str("data_list");

    // exp_sig/exp_err/exp_bkg/exp_data: treat bool true/false into a c_str
    auto resolve_exp = [&](const std::string &key) -> std::string {
        if (j_item.contains(key)) {
            auto &v = j_item[key];
            if (v.is_boolean())
                return v.get<bool>() ? kDrawAll : kNoDraw;
            if (v.is_string())
                return v.get<std::string>();
        }
        if (j_defaults.contains(key)) {
            auto &v = j_defaults[key];
            if (v.is_boolean())
                return v.get<bool>() ? kDrawAll : kNoDraw;
            if (v.is_string())
                return v.get<std::string>();
        }
        return "";
    };

    p.exp_sig  = resolve_exp("exp_sig");
    p.exp_err  = resolve_exp("exp_err");
    p.exp_bkg  = resolve_exp("exp_bkg");
    p.exp_data = resolve_exp("exp_data");

    // n_bin: if not present, use kNBin (matches original default)
    if (j_item.contains("n_bin"))
        p.n_bin = j_item["n_bin"].get<size_t>();
    else if (j_defaults.contains("n_bin"))
        p.n_bin = j_defaults["n_bin"].get<size_t>();
    else
        p.n_bin = kNBin;

    // bool normalize
    if (j_item.contains("normalize"))
        p.normalize = j_item["normalize"].get<bool>();
    else if (j_defaults.contains("normalize"))
        p.normalize = j_defaults["normalize"].get<bool>();
    else
        p.normalize = kNormalize;

    return p;
}

int main() {
    // Read JSON file
    std::ifstream jf(kJsonFile);
    if (!jf.is_open()) {
        std::cerr << "PANIC: fail to open " << kJsonFile << std::endl;
        std::abort();
    }
    json j;
    jf >> j;
    jf.close();

    json j_defaults = j.value("defaults", json::object());
    json j_items    = j.value("items",    json::array());

    // Determine kNProcess from defaults
    size_t kNProcess = kNProcessDefault;
    if (j_defaults.contains("kNProcess"))
        kNProcess = j_defaults["kNProcess"].get<size_t>();

    // Build params_map from JSON items
    std::vector<DrawParams> params_map;
    for (const auto &j_item : j_items) {
        params_map.push_back(json_to_draw_params(j_defaults, j_item));
    }

    size_t num_total = params_map.size();
    size_t num_work  = 0;
    size_t num_done  = 0;
    auto workItem = [num_total, &num_work, &num_done](const DrawParams &params) -> int {
        std::clog << "[INFO] " << ++num_work << "/" << num_total << " drawing " <<params.outname << std::endl;
        draw(params);
        std::clog << "[INFO] " << ++num_done << "/" << num_total << " writen " <<params.outname << std::endl;
        return 0;
    };
    ROOT::TProcessExecutor workers(kNProcess);
    workers.Map(workItem, params_map);
    std::clog << "[INFO] all graph drawn." << std::endl;
    return 0;
}

void draw(const DrawParams &params)
{
    auto chain_inmc = open_root(params.inmc_list);
    auto chain_data = open_root(params.data_list);
    std::vector<TH1D*> th1d_inmc;
    th1d_inmc.push_back(new TH1D("th1d_sig", "th1d_sig", params.n_bin, params.lo, params.up));
    th1d_inmc.push_back(new TH1D("th1d_err", "th1d_err", params.n_bin, params.lo, params.up));
    th1d_inmc.push_back(new TH1D("th1d_bkg", "th1d_bkg", params.n_bin, params.lo, params.up));
    if (!params.inmc_list.empty()) {
    std::cout << params.target << std::endl;
        chain_inmc->Draw((params.target+">>"+"th1d_sig").c_str(), params.exp_sig.c_str());
        chain_inmc->Draw((params.target+">>"+"th1d_err").c_str(), params.exp_err.c_str());
        chain_inmc->Draw((params.target+">>"+"th1d_bkg").c_str(), params.exp_bkg.c_str());
    }
    auto th1d_data = new TH1D("th1d_data", params.title.c_str(), params.n_bin, params.lo, params.up);
    if (!params.data_list.empty())
        chain_data->Draw((params.target+">>"+"th1d_data").c_str(), params.exp_data.c_str());
    
    if (params.normalize) {
        double scale = th1d_data->Integral()
            /(th1d_inmc[0]->Integral()
                +th1d_inmc[1]->Integral()
                +th1d_inmc[2]->Integral()
            );
        th1d_inmc[0]->Scale(scale);
        th1d_inmc[1]->Scale(scale);
        th1d_inmc[2]->Scale(scale);
    }
    THStack st("stack", params.title.c_str());
    th1d_inmc[0]->SetFillColor(2);
    th1d_inmc[1]->SetFillColor(4);
    th1d_inmc[2]->SetFillColor(14);
    st.Add(th1d_inmc[2]);
    st.Add(th1d_inmc[1]);
    st.Add(th1d_inmc[0]);
    TCanvas canvas("canvas", "canvas", 1);
    if (!params.inmc_list.empty()&&!params.data_list.empty()) {
        if (st.GetMaximum()<th1d_data->GetMaximum())
            st.SetMaximum(th1d_data->GetMaximum());
    }
    if (!params.inmc_list.empty()) st.Draw("hist");
    if (!params.data_list.empty()) {
        th1d_data->SetMarkerStyle(20);
        th1d_data->Draw("same e1");
    }
    canvas.Print(params.outname.c_str());
}
TChain* open_root(const std::string &filename, const std::string &treename) {
    auto result = new TChain(treename.c_str());
    auto treefriend = new TChain(treename.c_str());
    result->AddFriend(treefriend, "friend");
    if (filename.empty()) return result;
    if ( is_suffix_root(filename) ) {
        result->Add(filename.c_str());
        treefriend->Add((filename+".friend").c_str());
        return result;
    }
    std::ifstream infile(filename);
    std::string line;
    if (!infile.is_open()) PANIC("fail to open "+filename);
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t start = line.find_first_not_of(" \t");
        size_t end   = line.find_last_not_of(" \t");
        if (start == end) continue;
        std::string trimmed = line.substr(start, end-start+1);
        if ( is_suffix_root(trimmed) ) {
            result->Add(trimmed.c_str());
            treefriend->Add((trimmed+".friend").c_str());
        }
    }
    return result;
}