cmake_minimum_required(VERSION 3.30)
project(QuantitativeFuturesStrategy)

set(CMAKE_CXX_STANDARD 20)

find_package(CURL REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

add_executable(FuturesStrategyImproved main.cpp)

target_link_libraries(FuturesStrategyImproved PRIVATE CURL::libcurl nlohmann_json::nlohmann_json)
