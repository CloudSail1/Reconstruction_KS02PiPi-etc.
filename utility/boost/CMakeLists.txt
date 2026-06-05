cmake_minimum_required(VERSION 3.15...4.1)

project(boost LANGUAGES CXX)

# Finding the ROOT package
## [find_package]
find_package(ROOT 6.16 CONFIG REQUIRED)
## [find_package]

# 收集func文件夹下的所有.cxx文件
file(GLOB FUNC_SOURCES func/*.cxx)

# Adding an executable program and linking to needed ROOT libraries
## [add_and_link]
add_executable(boost ../read/read.cxx test.cxx ${FUNC_SOURCES})
target_link_libraries(boost PUBLIC ROOT::Core ROOT::Physics ROOT::Tree ROOT::Thread)
#add_library(read_test read.cxx)
## [add_and_link]