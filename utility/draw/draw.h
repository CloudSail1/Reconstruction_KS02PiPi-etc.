#include <string>
#include <ROOT/TProcessExecutor.hxx>
#include "TChain.h"

extern const std::string kTreeName;
/*
extern const std::string kInmcList;
extern const std::string kDataList;
extern const size_t      kNBin;
extern const double      kLo;
extern const double      kUp;
extern const std::string kTarget;
extern const std::string kTitle;
extern const std::string kOutName;
extern const std::string kExpSig;
extern const std::string kExpErr;
extern const std::string kExpBkg;
extern const std::string kExpData;
*/

struct DrawParams{
    double       lo        ;
    double       up        ;
    std::string  target    ;
    std::string  title     ;
    std::string  outname   ;
    std::string  inmc_list ;
    std::string  data_list ;
    std::string  exp_sig   ;
    std::string  exp_err   ;
    std::string  exp_bkg   ;
    std::string  exp_data  ;
    size_t       n_bin     ;
    bool         normalize ;
};
/*
= kLo;
= kUp;
= kTarget;
= kTitle;
= kOutName;
= kInmcList;
= kDataList;
= kExpSig;
= kExpErr;
= kExpBkg;
= kExpData;
= kNBin;
*/

#define PANIC(msg) \
    do { \
        std::cerr << "PANIC: " << msg \
                  << " at " << __FILE__ << ":" << __LINE__ \
                  << " in function " << __func__ << std::endl; \
        std::abort(); \
    } while(0)

TChain* open_root(const std::string &filename, const std::string &treename=kTreeName);
void draw(const DrawParams &params);
bool is_suffix_root(const std::string &str) {
    return str.size()>=5 && str.substr(str.size()-5)==".root";
}
