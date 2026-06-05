#include "../func.h"
#include <iostream>
double MinDAng(int size, double* dang) {
    if (size<=0) {return -99;}
    if (size>128) { // warning: see readme for bug message
        std::clog
            << "[WARN] dang size larger than pre-allocate size, some dang may be discard."
            << std::endl;
    }
    double min_dang = 9999;
    for (int i=0; i<size&&i<128; ++i) {
        if (min_dang > dang[i] && dang[i]!=0) {
            min_dang = dang[i];
        }
    }
    return min_dang;
}