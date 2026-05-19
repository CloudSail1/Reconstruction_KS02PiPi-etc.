#ifndef READ_H
#define READ_H
#include <string>
#include <memory>
#include "TChain.h"
#include <iostream>

std::unique_ptr<TChain> open_root(const std::string &filename,
    const std::string &treename,
    const std::string &friendname);
std::string name_friend(const std::string &str);
bool is_suffix_root(const std::string &str);
bool have_root(const std::string &str);
bool have_tree(const std::string &str, const std::string &treename);
bool have_friend(const std::string &str);
#endif
#ifndef PANIC
#define PANIC(msg) \
    do { \
        std::cerr << "PANIC " << msg \
                  << " at " << __FILE__ << ":" << __LINE__ \
                  << ": " << __func__ << std::endl; \
        std::abort(); \
    } while(0)
#endif
#ifndef INFO
#define INFO(msg) \
    do { \
        std::clog << "[INFO] " << msg \
                  << " at " << __FILE__ << ":" << __LINE__ \
                  << ": " << __func__ << std::endl; \
    } while(0)
#endif
#ifndef WARN
#define WARN(msg) \
    do { \
        std::clog << "[WARN] " << msg \
                  << " at " << __FILE__ << ":" << __LINE__ \
                  << ": " << __func__ << std::endl; \
    } while(0)
#endif