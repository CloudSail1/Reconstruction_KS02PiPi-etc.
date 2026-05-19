cmake_minimum_required(VERSION 3.15...4.1)

project(draw LANGUAGES CXX)

# Finding the ROOT package
## [find_package]
find_package(ROOT 6.16 CONFIG REQUIRED)
## [find_package]

# Adding an executable program and linking to needed ROOT libraries
## [add_and_link]
add_executable(draw main.cxx)
target_link_libraries(draw PUBLIC ROOT::Core ROOT::Physics ROOT::Tree ROOT::Thread ROOT::Gpad)
#add_library(read_test read.cxx)
## [add_and_link]