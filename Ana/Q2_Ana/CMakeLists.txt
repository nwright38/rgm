add_library(many_plots many_plots.cpp)
target_link_libraries(many_plots ${ROOT_LIBRARIES})
foreach(fnameSrc KMiss.cpp Analysis_Note_Q2.cpp Analysis_Note_Q2_comp.cpp Analysis_Note_Q2_manyTheory.cpp KMiss_comparison_plot.cpp p_LUND.cpp Electron_DeltaE.cpp Radiation_Offset.cpp Central_Detector_Proton.cpp iso_p_LUND.cpp Main_Figs.cpp Main_sigmaCM.cpp Main_sigmaCM_Hists.cpp GCF_check.cpp Main_Figs_Sys_Err.cpp MMiss.cpp For_Florian.cpp Main_ForJulian.cpp ratio_ForJulian.cpp skim_elastic.cpp edp_lowxB.cpp pd_new.cpp)
  message(STATUS ${fnameSRC})
  string (REPLACE ".cpp" "" fnameExe ${fnameSrc})
  add_executable(${fnameExe} ${fnameSrc})
  target_link_libraries(${fnameExe} ${ROOT_LIBRARIES} PkgConfig::hipo4 -lEG -lClas12Banks -lClas12Root -L${CLAS12ROOT}/lib -L${CLAS12ROOT}/ccdb/lib Clas12Ana Clas12Debug reweighter many_plots)
endforeach()
# Analysis_Note_Q2_manyTheory.cpp
