cmake_minimum_required(VERSION 3.16)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(SigmaCM LANGUAGES CXX)
endif()

if(NOT TARGET ROOT::Core OR NOT TARGET ROOT::Minuit2)
  find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist MathCore Minuit2)
endif()

add_library(sigmacm
  SigmaCMConfig.cpp
  SigmaCMExtractor.cpp
  SigmaCMSkimIO.cpp
  SigmaCMResultIO.cpp
)
target_include_directories(sigmacm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sigmacm PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::MathCore ROOT::Minuit2)
target_compile_features(sigmacm PUBLIC cxx_std_17)

foreach(app
  run_nominal
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
