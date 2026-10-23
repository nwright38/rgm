// Standard GCF genT -> generator-level SRC kinematics (73 skim-compatible branches).
#include "convert_events2N_to_srcTree.C"

void convert_genT_to_srcTree(
    const char *inputFileName,
    const char *outputFileName,
    int targetA = 12,
    double eBeam = -1.0,
    bool requireEpp = false,
    bool applyMCSmearing = true,
    bool applyBasicSrcCuts = false,
    const char *inputTreeName = "genT",
    double minXB = 1.2,
    double minQ2 = 1.5,
    double minLeadP = 0.8,
    double minRecP = 0.3,
    double minKMiss = 0.3,
    double maxKMiss = 1.0,
    double minMMiss = 0.0,
    double maxMMiss = 1.1) {
  convert_events2N_to_srcTree(inputFileName, outputFileName, inputTreeName,
      true, targetA, eBeam, requireEpp, applyMCSmearing, applyBasicSrcCuts,
      minXB, minQ2, minLeadP, minRecP, minKMiss, maxKMiss,
      minMMiss, maxMMiss, true);
}
