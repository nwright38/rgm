#include "src_skim_to_plot_common.C"

void build_src_skim_stat_sys(
    const char* inputFile = "~/data/RGM_DATA/c12_sim_skim_noMmiss.root",
    const char* outputFile = "~/data/RGM_DATA/c12_sim_stat_sys_plots.root",
    const char* treeName = "srcTree",
    int nToys = 100,
    unsigned int toySeed = 1234,
  const char* leadMode = "both",
  int massNumber = 0) {
  src_skim_plot::CutConfig cuts = src_skim_plot::DefaultCutConfig();

  // Edit cuts here.
  // Example:
  // cuts.usePcmMinusPrelMax = true;
  // cuts.pcmMinusPrel_max = 0.0;          // pCM - pRel < 0  (equivalent to pCM < pRel)
  // cuts.pcmMinusPrel_max_sigma = 0.02;   // toy spread for that threshold

  src_skim_plot::runPipeline(
      inputFile,
      outputFile,
      true,
      nToys,
      toySeed,
      treeName,
      leadMode,
      cuts,
      massNumber);
}
