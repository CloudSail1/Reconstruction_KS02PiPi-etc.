#include "read.h"

#include "TSystem.h"
#include "TFile.h"
#include <fstream>

std::string name_friend(const std::string &str) {
    return std::string(str+".friend");
}
std::unique_ptr<TChain> open_root(
    const std::string &filename,
    const std::string &treename,
    const std::string &friendname
) {
    bool is_read_friend = !friendname.empty();
    std::unique_ptr<TChain> result(new TChain(treename.c_str()));
    std::unique_ptr<TChain> treefriend(new TChain(treename.c_str()));
    if (is_read_friend)
        result->AddFriend(treefriend.get(), friendname.c_str());
    if (filename.empty()) return result; // TODO add WARN
    if ( is_suffix_root(filename) ) {
        if ( have_tree(filename, treename) ) {
            result->Add(filename.c_str());
            INFO("Adding to chain: "+filename);
        } else {
            WARN("Failed open Tree cause by "+
                (have_root(filename)
                    ?"no tree "+treename+"in file "+filename
                    :"failing to open root "+filename));
        }
        if ( is_read_friend ) {
            if ( have_friend(name_friend(filename)) ) {
                treefriend->Add((name_friend(filename)).c_str());
                INFO("Adding as friend: "+name_friend(filename));
            } else {
                WARN("friend no found with name: "+name_friend(filename));
            }
        }
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
            if ( have_tree(trimmed, treename) ) {
                result->Add(trimmed.c_str());
                INFO("Adding to chain: "+trimmed);
            } else {
                WARN("Failed open Tree cause by "+
                    (have_root(trimmed)
                        ?"no tree "+treename+"in file "+trimmed
                        :"failing to open root "+trimmed));
            }
            if ( is_read_friend ) {
                if ( have_friend(name_friend(trimmed)) ) {
                    treefriend->Add((name_friend(trimmed)).c_str());
                    INFO("Adding as friend: "+name_friend(trimmed));
                } else {
                    WARN("friend no found with name: "+name_friend(trimmed));
                }
            }
        }
    }
    INFO("TChain has readed: "+filename);
    if (result->GetEntries()<=0) PANIC("TChain is empty: "+filename);
    return result;
}
bool is_suffix_root(const std::string &str) {
    return str.size()>=5 && str.substr(str.size()-5)==".root";
}
bool have_root(const std::string &str) {
    if (gSystem->AccessPathName(str.c_str(),kReadPermission))
        return false;
    std::unique_ptr<TFile> file( TFile::Open(str.c_str()) );
    return file && !file->IsZombie();
}
bool have_tree(const std::string &str, const std::string &treename) {
    if (!have_root(str))
        return false;
    TFile file(str.c_str());
    auto tree = (TTree*)file.GetObjectChecked(treename.c_str(), "TTree");
    if (tree) return true;
    else return false;
}
bool have_friend(const std::string &str) {
    if (gSystem->AccessPathName(name_friend(str).c_str(),kReadPermission))
        return false;
    std::unique_ptr<TFile> file( TFile::Open(name_friend(str).c_str()) );
    return file && !file->IsZombie();
}
