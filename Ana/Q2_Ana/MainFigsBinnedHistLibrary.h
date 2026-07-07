#ifndef MAINFIGSBINNEDHISTLIBRARY_HH
#define MAINFIGSBINNEDHISTLIBRARY_HH

#include <functional>
#include <string>
#include <vector>

#include "BinnedHistStore.h"

template <typename EventT>
struct OneDHistDef {
  std::string name;
  Selection selection;
  std::function<bool(const EventT&)> passFn;
  std::function<double(const EventT&)> valueFn;
  int nbins;
  double xmin;
  double xmax;
  std::vector<double> valueEdges;
};

template <typename EventT>
void appendOneDHist(std::vector<FillTask<EventT>>& tasks, const OneDHistDef<EventT>& def) {
  tasks.push_back({def.name, def.selection, def.passFn, {}, {},
                   def.valueFn, def.nbins, def.xmin, def.xmax, def.valueEdges});
}

template <typename EventT>
void appendOneDHistList(std::vector<FillTask<EventT>>& tasks,
                        const std::vector<OneDHistDef<EventT>>& defs) {
  for (const auto& def : defs) appendOneDHist(tasks, def);
}

#endif
