cmake_minimum_required(VERSION 3.16)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(SigmaCM LANGUAGES CXX)
endif()

if(NOT TARGET ROOT::Core OR NOT TARGET ROOT::Minuit2 OR NOT TARGET ROOT::Gpad)
  find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist Gpad MathCore Minuit2)
endif()

add_library(sigmacm
  SigmaCMConfig.cpp
  SigmaCMExtractor.cpp
  SigmaCMPlotOutput.cpp
  SigmaCMInput.cpp
  SigmaCMResultIO.cpp
)
target_include_directories(sigmacm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sigmacm PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::MathCore ROOT::Minuit2)
target_compile_features(sigmacm PUBLIC cxx_std_17)

if(TARGET Clas12Ana AND TARGET Clas12Debug AND TARGET PkgConfig::hipo4)
  target_compile_definitions(sigmacm PUBLIC SIGMACM_WITH_HIPO)
  target_link_libraries(sigmacm PUBLIC
    PkgConfig::hipo4
    -lEG
    -lClas12Banks
    -lClas12Root
    -L${CLAS12ROOT}/lib
    -L${CLAS12ROOT}/ccdb/lib
    Clas12Ana
    Clas12Debug
  )
endif()

add_executable(sigmacm_extract run_extract.cpp)
target_link_libraries(sigmacm_extract PRIVATE sigmacm)

add_executable(sigmacm_make_skim run_make_skim.cpp)
target_link_libraries(sigmacm_make_skim PRIVATE sigmacm)

foreach(app
  run_cut_toys
  run_gcf_toys
  run_combined_toys
  run_closure
  run_fit_range_scan
  run_diagnostics
  run_profile_scan
)
  add_executable(sigmacm_${app} ${app}.cpp)
  target_link_libraries(sigmacm_${app} PRIVATE sigmacm)
endforeach()

add_executable(sigmacm_unit_tests test_sigmacm.cpp)
target_link_libraries(sigmacm_unit_tests PRIVATE sigmacm)

add_custom_target(sigmacm_tests
  COMMAND sigmacm_unit_tests
  DEPENDS sigmacm_unit_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
