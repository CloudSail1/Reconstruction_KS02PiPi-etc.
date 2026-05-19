#include "read.h"

int main() {
    std::cout << "hello world" << std::endl;
    auto chain = open_root("/scratchfs/bes/cloudsail/JpsiToPiKKshort/feat-PiPiandGamma/run/data/data_1.root", "chi2omg", "");
    chain->Print();
    return 0;
}