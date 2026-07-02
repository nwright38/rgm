#ifndef Q2_REWEIGHT_H
#define Q2_REWEIGHT_H

#include <iostream>
#include <memory>
#include <string>

#include <TFile.h>
#include <TH1.h>

class Q2Reweight {
public:
  Q2Reweight() = default;

  explicit Q2Reweight(const std::string& fileName) {
    load(fileName);
  }

  bool load(const std::string& fileName) {
    file_.reset(TFile::Open(fileName.c_str(), "READ"));
    if (!file_ || file_->IsZombie()) {
      std::cerr << "Could not open Q2 reweight file: " << fileName << std::endl;
      file_.reset();
      return false;
    }

    TH1* h = nullptr;
    file_->GetObject("q2_reweighting/h_Q2_reweight", h);
    if (!h) file_->GetObject("h_Q2_reweight", h);
    if (!h) {
      std::cerr << "Could not find h_Q2_reweight in " << fileName << std::endl;
      file_.reset();
      return false;
    }

    weights_.reset(static_cast<TH1*>(h->Clone("loaded_h_Q2_reweight")));
    weights_->SetDirectory(nullptr);
    enabled_ = true;
    return true;
  }

  bool enabled() const { return enabled_ && weights_; }

  double weight(double q2) const {
    if (!enabled()) return 1.0;
    const int bin = weights_->FindFixBin(q2);
    if (bin < 1 || bin > weights_->GetNbinsX()) return 0.0;
    return weights_->GetBinContent(bin);
  }

private:
  std::unique_ptr<TFile> file_;
  std::unique_ptr<TH1> weights_;
  bool enabled_ = false;
};

#endif
