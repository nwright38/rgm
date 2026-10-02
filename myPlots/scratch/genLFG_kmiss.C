// =====================================================================================
//  genLFG_kmiss.C
//
//  Mean-field (e,e'p) event generator using a Local Fermi Gas (LFG) spectral function,
//  with detector smearing, writing truth + reconstructed pMiss and kMiss to a ROOT tree.
//
//  Physics
//  -------
//   * Struck proton position r is sampled from the proton density rho_p(r) r^2
//     (modified harmonic oscillator for A<=16, Woods-Saxon otherwise).
//   * Initial momentum is uniform inside the local Fermi sphere,
//       k_F(r) = hbar*c * (3 pi^2 rho_p(r))^(1/3).
//   * Removal energy  E_m = S_p + [k_F(r)^2 - p^2]/(2 m_N)  (hole depth below the local
//     Fermi surface; set cfg.useHoleEnergy=false for a constant E_m = S_p).
//   * Residual A-1 system: M*_{A-1} = M_A - m_p + E_m, recoiling with -p_i.
//   * Electron angles are sampled uniformly in (cos theta_e, phi_e); E' is then fixed by
//     exact 4-momentum conservation (unique root, solved by bisection).
//   * PWIA event weight [nb]:  w = Z * DeltaOmega_e * sigma_cc1 * |dE'/dE_m|
//     (de Forest CC1 off-shell e-p cross section, Kelly 2004 proton form factors).
//     Sum(w)/N estimates the cross section in nb for the generated electron phase space.
//     Histograms of resolution *should be filled with this weight* so that the kinematic
//     mix (Q2, x_B, proton angle) looks like data.
//
//  Missing-momentum variables (paper definitions)
//  ----------------------------------------------
//     p_miss^mu = p_lead^mu - q^mu,   E_miss = E_lead - omega
//     p_miss^-  = E_miss - p_miss,z        (z along q-vector)
//     k_miss^2  = m_N^2 (p_miss,perp^2 + m_N^2) / [p_miss^- (2 m_N - p_miss^-)] - m_N^2
//
//  Detector smearing (EDIT the Config block below for your setup)
//  -------------------------------------------------------------------
//     sigma_p/p   = sqrt( a^2 + (b*p)^2 )       (FD: total p;  CD: transverse p_T)
//     sigma_theta = sqrt( c^2 + (d/p)^2 )       [rad]
//     sigma_phi   = sqrt( e^2 + (f/p)^2 )       [rad]
//   Separate parameter sets for: electron in FD, proton in CD, proton in FD.
//   Region is decided from the TRUE proton polar angle.
//
//  Usage
//  -----
//     root -l -b -q 'genLFG_kmiss.C+(1000000,"lfg_C12.root","C12",5.986)'
//     root -l -b -q 'genLFG_kmiss.C+(1000000,"lfg_C12.root","C12",5.986,12345,1.5,1.2)'
//   Targets: "He4","C12","O16","Ar40","Ca40","Ca48","Sn120"
//
//   (Testing without ROOT:  g++ -O2 -DLFG_NO_ROOT -x c++ genLFG_kmiss.C -o gen
//    writes a CSV instead of a tree.)
// =====================================================================================

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#include <algorithm>

#ifndef LFG_NO_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

namespace lfg {

// ---------------------------------------------------------------- constants (GeV, fm)
const double mP     = 0.93827208816;
const double mE     = 0.00051099895;
const double hbarc  = 0.1973269804;       // GeV fm
const double alphaEM= 1.0 / 137.035999084;
const double amu    = 0.93149410242;      // GeV
const double GeV2nb = 389379.37;          // 1 GeV^-2 in nb
const double muP    = 2.79284734;
const double deg    = M_PI / 180.0;

// ---------------------------------------------------------------- tiny 3-vector
struct V3 {
  double x = 0, y = 0, z = 0;
  V3() {}
  V3(double a, double b, double c) : x(a), y(b), z(c) {}
  V3 operator+(const V3 &o) const { return V3(x + o.x, y + o.y, z + o.z); }
  V3 operator-(const V3 &o) const { return V3(x - o.x, y - o.y, z - o.z); }
  V3 operator*(double s) const { return V3(x * s, y * s, z * s); }
  double Dot(const V3 &o) const { return x * o.x + y * o.y + z * o.z; }
  V3 Cross(const V3 &o) const { return V3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x); }
  double Mag2() const { return Dot(*this); }
  double Mag() const { return std::sqrt(Mag2()); }
  V3 Unit() const { double m = Mag(); return m > 0 ? (*this) * (1.0 / m) : V3(); }
  double Theta() const { return std::atan2(std::sqrt(x * x + y * y), z); }
  double Phi() const { return std::atan2(y, x); }
  double Perp() const { return std::sqrt(x * x + y * y); }
};
inline V3 Sph(double p, double th, double ph) {
  return V3(p * std::sin(th) * std::cos(ph), p * std::sin(th) * std::sin(ph), p * std::cos(th));
}

// ---------------------------------------------------------------- resolution model
struct Res {
  double a, b;   // sigma_p/p   = sqrt(a^2 + (b p)^2)     (p in GeV)
  double c, d;   // sigma_theta = sqrt(c^2 + (d/p)^2)     (rad, rad*GeV)
  double e, f;   // sigma_phi   = sqrt(e^2 + (f/p)^2)     (rad, rad*GeV)
  bool   usePT;  // true: momentum resolution applies to p_T (central/solenoid tracker)
};

// ---------------------------------------------------------------- configuration
struct Config {
  // beam / target
  double Ebeam  = 5.986;      // GeV
  int    A = 12, Z = 6;
  double massU  = 12.0;       // atomic mass [u]
  double Sp     = 0.015957;   // proton separation energy [GeV]
  int    densModel = 1;       // 0 = Woods-Saxon, 1 = modified harmonic oscillator
  double hoA = 1.687, hoAlpha = 1.067;          // MHO: rho ~ (1+alpha (r/a)^2) exp(-(r/a)^2)
  double wsR = 0, wsA = 0.54;                    // WS : rho ~ 1/(1+exp((r-R)/a)), R set from A if 0
  bool   useHoleEnergy = true;

  // generated electron phase space
  double thetaEmin = 5.0 * deg, thetaEmax = 35.0 * deg;
  // truth-level generation cuts (set <=0 to disable)
  double Q2min = 0.0;   // GeV^2
  double xBmin = 0.0;

  // acceptance regions (true proton polar angle) -- flags only, no events are cut
  double thetaFDmin = 5.0 * deg,  thetaFDmax = 35.0 * deg;
  double thetaCDmin = 35.0 * deg, thetaCDmax = 125.0 * deg;

  // ===== DETECTOR RESOLUTIONS: EDIT THESE ==========================================
  // Placeholder CLAS12-like values. Replace with your detector numbers.
  //                  a       b      c       d       e       f      usePT
  Res eFD      = {0.0070, 0.0,   0.0010, 0.0,   0.0030, 0.0,   false};
  Res pFD      = {0.0070, 0.0,   0.0010, 0.0,   0.0030, 0.0,   false};
  Res pCD      = {0.0500, 0.0,   0.0150, 0.0,   0.0050, 0.0,   true };
  // =================================================================================
};

inline bool SetTarget(Config &c, const std::string &t) {
  c.densModel = 0; c.wsR = 0; c.wsA = 0.54;
  if (t == "He4")        { c.A = 4;   c.Z = 2;  c.massU = 4.00260325413; c.Sp = 0.019814; c.densModel = 1; c.hoA = 1.372; c.hoAlpha = 0.0; }
  else if (t == "C12")   { c.A = 12;  c.Z = 6;  c.massU = 12.0;          c.Sp = 0.015957; c.densModel = 1; c.hoA = 1.687; c.hoAlpha = 1.067; }
  else if (t == "O16")   { c.A = 16;  c.Z = 8;  c.massU = 15.9949146196; c.Sp = 0.012127; c.densModel = 1; c.hoA = 1.833; c.hoAlpha = 1.544; }
  else if (t == "Ar40")  { c.A = 40;  c.Z = 18; c.massU = 39.9623831237; c.Sp = 0.012529; }
  else if (t == "Ca40")  { c.A = 40;  c.Z = 20; c.massU = 39.962590863;  c.Sp = 0.008328; }
  else if (t == "Ca48")  { c.A = 48;  c.Z = 20; c.massU = 47.95252276;   c.Sp = 0.015806; }
  else if (t == "Sn120") { c.A = 120; c.Z = 50; c.massU = 119.90220163;  c.Sp = 0.010690; }
  else return false;
  return true;
}

// ---------------------------------------------------------------- nuclear density + LFG
struct Density {
  std::vector<double> r, cdf, kF;   // r grid [fm], CDF of rho_p r^2, local proton k_F [GeV]
  void Build(const Config &c) {
    const int n = 4000;
    double R = (c.wsR > 0) ? c.wsR : 1.12 * std::cbrt(c.A) - 0.86 / std::cbrt(c.A);
    double rmax = (c.densModel == 1) ? 6.0 * c.hoA : R + 15.0 * c.wsA;
    r.resize(n + 1); cdf.resize(n + 1); kF.resize(n + 1);
    std::vector<double> rho(n + 1);
    double dr = rmax / n, norm = 0;
    for (int i = 0; i <= n; i++) {
      r[i] = i * dr;
      double x = r[i];
      if (c.densModel == 1) { double u = x / c.hoA; rho[i] = (1 + c.hoAlpha * u * u) * std::exp(-u * u); }
      else                  rho[i] = 1.0 / (1.0 + std::exp((x - R) / c.wsA));
      if (i > 0) norm += 0.5 * (rho[i] * x * x + rho[i - 1] * r[i - 1] * r[i - 1]) * dr;
    }
    norm *= 4 * M_PI;                               // integral of shape
    double scaleP = double(c.Z) / norm;             // rho_p normalised to Z protons  [fm^-3]
    cdf[0] = 0;
    for (int i = 0; i <= n; i++) {
      double rhoP = scaleP * rho[i];
      kF[i] = hbarc * std::cbrt(3 * M_PI * M_PI * rhoP);
      if (i > 0) cdf[i] = cdf[i - 1] + 0.5 * (rho[i] * r[i] * r[i] + rho[i - 1] * r[i - 1] * r[i - 1]) * dr;
    }
    for (auto &v : cdf) v /= cdf[n];
  }
  // returns r, sets kFr
  double Sample(double u, double &kFr) const {
    int i = int(std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin());
    if (i <= 0) { kFr = kF[0]; return 0; }
    double f = (u - cdf[i - 1]) / std::max(1e-300, cdf[i] - cdf[i - 1]);
    kFr = kF[i - 1] + f * (kF[i] - kF[i - 1]);
    return r[i - 1] + f * (r[i] - r[i - 1]);
  }
};

// ---------------------------------------------------------------- form factors & CC1
inline void KellyFF(double Q2, double &GE, double &GM) {
  double t = Q2 / (4 * mP * mP);
  GE = (1 - 0.24 * t) / (1 + 10.98 * t + 12.82 * t * t + 21.97 * t * t * t);
  GM = muP * (1 + 0.12 * t) / (1 + 10.97 * t + 18.86 * t * t + 6.55 * t * t * t);
}

// de Forest CC1 off-shell e-p cross section, GeV^-2 / sr  (no recoil factor)
inline double SigmaCC1(double Eb, const V3 &kIn, const V3 &kOut, double Ep_e,
                       const V3 &pLead, double ELead, const V3 &pm) {
  V3 q = kIn - kOut;
  double omega = Eb - Ep_e, q2 = q.Mag2(), qm = std::sqrt(q2);
  double Q2 = q2 - omega * omega;
  double thE = std::acos(std::max(-1.0, std::min(1.0, kIn.Unit().Dot(kOut.Unit()))));
  double s2 = std::sin(thE / 2), c2 = std::cos(thE / 2), tan2 = (s2 * s2) / (c2 * c2);
  double sigM = alphaEM * alphaEM * c2 * c2 / (4 * Eb * Eb * s2 * s2 * s2 * s2);

  double Ebar = std::sqrt(mP * mP + pm.Mag2());
  double wbar = ELead - Ebar;
  double Qb2  = q2 - wbar * wbar;
  double GE, GM; KellyFF(Q2, GE, GM);
  double tau = Q2 / (4 * mP * mP);
  double F1 = (GE + tau * GM) / (1 + tau);
  double kF2 = (GM - GE) / (1 + tau);                       // kappa*F2
  double Aff = F1 * F1 + Qb2 / (4 * mP * mP) * kF2 * kF2;
  double Bff = (F1 + kF2) * (F1 + kF2);

  double pp = pLead.Mag();
  double cosG = pLead.Unit().Dot(q.Unit());
  double sinG = std::sqrt(std::max(0.0, 1 - cosG * cosG));
  V3 n1 = q.Cross(kIn), n2 = q.Cross(pLead);
  double cosPhi = (n1.Mag() > 0 && n2.Mag() > 0) ? n1.Unit().Dot(n2.Unit()) : 1.0;

  double EE = Ebar * ELead;
  double wC = ((Ebar + ELead) * (Ebar + ELead) * Aff - q2 * Bff) / (4 * EE);
  double wT = Qb2 * Bff / (2 * EE);
  double wS = pp * pp * sinG * sinG * Aff / EE;
  double wI = -pp * sinG * (Ebar + ELead) * Aff / EE;
  double lam = Q2 / q2;
  double s = lam * lam * wC + (lam / 2 + tan2) * wT + (lam * cosPhi * cosPhi + tan2) * wS
           + lam * std::sqrt(lam + tan2) * cosPhi * wI;
  (void)qm;
  return sigM * std::max(0.0, s);
}

// ---------------------------------------------------------------- missing-momentum variables
struct Miss {
  double Q2, xB, omega, q, pMiss, Emiss, pPerp, pMinus, kMiss;
  bool   kValid;
};
inline Miss ComputeMiss(double Eb, const V3 &kOut, const V3 &pLead, double mN) {
  Miss m;
  V3 kIn(0, 0, Eb);
  double Ep = std::sqrt(kOut.Mag2() + mE * mE);
  V3 q = kIn - kOut;
  m.omega = Eb - Ep;
  m.q = q.Mag();
  m.Q2 = q.Mag2() - m.omega * m.omega;
  m.xB = m.Q2 / (2 * mP * m.omega);
  V3 pm = pLead - q;
  double ELead = std::sqrt(pLead.Mag2() + mP * mP);
  m.pMiss = pm.Mag();
  m.Emiss = ELead - m.omega;
  V3 qu = q.Unit();
  double pz = pm.Dot(qu);
  m.pPerp = (pm - qu * pz).Mag();
  m.pMinus = m.Emiss - pz;
  double den = m.pMinus * (2 * mN - m.pMinus);
  if (den > 0) {
    double k2 = mN * mN * (m.pPerp * m.pPerp + mN * mN) / den - mN * mN;
    m.kMiss = std::sqrt(std::max(0.0, k2));
    m.kValid = true;
  } else { m.kMiss = -1; m.kValid = false; }
  return m;
}

// ---------------------------------------------------------------- smearing
inline V3 Smear(const V3 &p, const Res &r, std::mt19937_64 &g) {
  std::normal_distribution<double> N(0, 1);
  double pm = p.Mag(), th = p.Theta(), ph = p.Phi();
  double sp  = std::sqrt(r.a * r.a + (r.b * pm) * (r.b * pm));
  double sth = std::sqrt(r.c * r.c + (r.d / pm) * (r.d / pm));
  double sph = std::sqrt(r.e * r.e + (r.f / pm) * (r.f / pm));
  double thS = th + sth * N(g), phS = ph + sph * N(g);
  if (r.usePT) {
    double pT = pm * std::sin(th);
    double pTs = pT * (1 + std::sqrt(r.a * r.a + (r.b * pT) * (r.b * pT)) * N(g));
    double ps = pTs / std::max(1e-6, std::sin(thS));
    return Sph(ps, thS, phS);
  }
  double ps = pm * (1 + sp * N(g));
  return Sph(ps, thS, phS);
}

} // namespace lfg

// =====================================================================================
//  main generator
// =====================================================================================
#ifdef LFG_NO_ROOT
typedef long long Long64_t;
#endif

void genLFG_kmiss(Long64_t nEvents = 1000000, const char *outName = "lfg_kmiss.root",
                  const char *target = "C12", double Ebeam = 5.986, unsigned seed = 12345,
                  double Q2min = 1., double xBmin = 1.0)
{
  using namespace lfg;
  Config cfg;
  cfg.Ebeam = Ebeam;
  cfg.Q2min = Q2min;
  cfg.xBmin = xBmin;
  if (!SetTarget(cfg, target)) { printf("Unknown target %s\n", target); return; }

  const double mN = mP;                                      // m_N used in the k_miss formula
  const double MA = cfg.massU * amu - cfg.Z * mE;            // nuclear mass
  Density dens; dens.Build(cfg);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);

  const V3 kIn(0, 0, cfg.Ebeam);
  const double cmin = std::cos(cfg.thetaEmax), cmax = std::cos(cfg.thetaEmin);
  const double dOmega = 2 * M_PI * (cmax - cmin);

  // --------------------------------------------------------------- output
  double r_fm, kF, Em, weight;
  double Q2_t, xB_t, omega_t, q_t, pMiss_t, kMiss_t, Emiss_t, pPerp_t, pMinus_t;
  double Q2_r, xB_r, omega_r, q_r, pMiss_r, kMiss_r, Emiss_r, pPerp_r, pMinus_r;
  double pe_t, the_t, phe_t, pp_t, thp_t, php_t;
  double pe_r, the_r, phe_r, pp_r, thp_r, php_r;
  int    region, kValid_t, kValid_r;

#ifndef LFG_NO_ROOT
  TFile *fout = TFile::Open(outName, "RECREATE");
  TTree *T = new TTree("T", Form("LFG mean-field (e,e'p), %s, E=%.3f GeV", target, Ebeam));
  // requested variables
  T->Branch("kMiss_truth", &kMiss_t);  T->Branch("kMiss", &kMiss_r);
  T->Branch("pMiss_truth", &pMiss_t);  T->Branch("pMiss", &pMiss_r);
  // everything else useful for cuts / diagnostics
  T->Branch("weight", &weight);        T->Branch("region", &region);   // 1=CD, 2=FD, 0=none
  T->Branch("kValid_truth", &kValid_t); T->Branch("kValid", &kValid_r);
  T->Branch("r_fm", &r_fm); T->Branch("kF", &kF); T->Branch("Em_true", &Em);
  T->Branch("Q2_truth", &Q2_t);   T->Branch("Q2", &Q2_r);
  T->Branch("xB_truth", &xB_t);   T->Branch("xB", &xB_r);
  T->Branch("omega_truth", &omega_t); T->Branch("omega", &omega_r);
  T->Branch("q_truth", &q_t);     T->Branch("q", &q_r);
  T->Branch("Emiss_truth", &Emiss_t); T->Branch("Emiss", &Emiss_r);
  T->Branch("pMissPerp_truth", &pPerp_t); T->Branch("pMissPerp", &pPerp_r);
  T->Branch("pMinus_truth", &pMinus_t);   T->Branch("pMinus", &pMinus_r);
  T->Branch("pe_truth", &pe_t);  T->Branch("thetae_truth", &the_t);  T->Branch("phie_truth", &phe_t);
  T->Branch("pe", &pe_r);        T->Branch("thetae", &the_r);        T->Branch("phie", &phe_r);
  T->Branch("pp_truth", &pp_t);  T->Branch("thetap_truth", &thp_t);  T->Branch("phip_truth", &php_t);
  T->Branch("pp", &pp_r);        T->Branch("thetap", &thp_r);        T->Branch("phip", &php_r);
#else
  FILE *fout = fopen(outName, "w");
  fprintf(fout, "weight,region,kMiss_truth,kMiss,pMiss_truth,pMiss,Q2_truth,xB_truth,Emiss_truth,"
                "pp_truth,thetap_truth,pe_truth,thetae_truth,Em_true,kF\n");
#endif

  Long64_t nFilled = 0, nTry = 0;
  while (nFilled < nEvents) {
    nTry++;
    // ---- initial proton from LFG
    r_fm = dens.Sample(U(rng), kF);
    double p0 = kF * std::cbrt(U(rng));
    double ct = 2 * U(rng) - 1, phi0 = 2 * M_PI * U(rng);
    V3 pInit = Sph(p0, std::acos(ct), phi0);
    Em = cfg.Sp + (cfg.useHoleEnergy ? (kF * kF - p0 * p0) / (2 * mP) : 0.0);
    double Mstar = MA - mP + Em;
    double Erec  = std::sqrt(Mstar * Mstar + p0 * p0);

    // ---- electron direction
    double cte = cmin + (cmax - cmin) * U(rng), phe = 2 * M_PI * U(rng);
    V3 nE = Sph(1.0, std::acos(cte), phe);

    // ---- solve energy conservation for E'
    auto f = [&](double Epr) {
      V3 kOut = nE * std::sqrt(std::max(0.0, Epr * Epr - mE * mE));
      V3 pL = kIn - kOut + pInit;
      return cfg.Ebeam - Epr + MA - std::sqrt(pL.Mag2() + mP * mP) - Erec;
    };
    double lo = 2 * mE, hi = cfg.Ebeam;
    if (f(lo) <= 0) continue;                       // kinematically forbidden
    for (int it = 0; it < 80; it++) { double mid = 0.5 * (lo + hi); (f(mid) > 0 ? lo : hi) = mid; }
    double Epr = 0.5 * (lo + hi);
    V3 kOut = nE * std::sqrt(Epr * Epr - mE * mE);
    V3 pLead = kIn - kOut + pInit;
    double ELead = std::sqrt(pLead.Mag2() + mP * mP);

    // ---- weight: Z * dOmega * sigma_cc1 * |dE'/dE_m|
    double dEpdEpr = -(pLead.Dot(nE)) / ELead * (Epr / std::sqrt(Epr * Epr - mE * mE));
    double dEmdEpr = -(1.0 + dEpdEpr) * Erec / Mstar;
    double sig = SigmaCC1(cfg.Ebeam, kIn, kOut, Epr, pLead, ELead, pInit);
    weight = cfg.Z * dOmega * sig / std::fabs(dEmdEpr) * GeV2nb;

    // ---- truth variables
    Miss mt = ComputeMiss(cfg.Ebeam, kOut, pLead, mN);
    if (cfg.Q2min > 0.0 && mt.Q2 < cfg.Q2min) continue;
    if (cfg.xBmin > 0.0 && mt.xB < cfg.xBmin) continue;
    thp_t = pLead.Theta();
    region = (thp_t >= cfg.thetaCDmin && thp_t <= cfg.thetaCDmax) ? 1
           : (thp_t >= cfg.thetaFDmin && thp_t <  cfg.thetaFDmax) ? 2 : 0;

    // ---- smear
    V3 kOutR  = Smear(kOut, cfg.eFD, rng);
    V3 pLeadR = Smear(pLead, region == 2 ? cfg.pFD : cfg.pCD, rng);
    Miss mr = ComputeMiss(cfg.Ebeam, kOutR, pLeadR, mN);

    // ---- fill
    Q2_t = mt.Q2; xB_t = mt.xB; omega_t = mt.omega; q_t = mt.q; pMiss_t = mt.pMiss; kMiss_t = mt.kMiss;
    Emiss_t = mt.Emiss; pPerp_t = mt.pPerp; pMinus_t = mt.pMinus; kValid_t = mt.kValid;
    Q2_r = mr.Q2; xB_r = mr.xB; omega_r = mr.omega; q_r = mr.q; pMiss_r = mr.pMiss; kMiss_r = mr.kMiss;
    Emiss_r = mr.Emiss; pPerp_r = mr.pPerp; pMinus_r = mr.pMinus; kValid_r = mr.kValid;
    pe_t = kOut.Mag();  the_t = kOut.Theta();  phe_t = kOut.Phi();
    pe_r = kOutR.Mag(); the_r = kOutR.Theta(); phe_r = kOutR.Phi();
    pp_t = pLead.Mag(); php_t = pLead.Phi();
    pp_r = pLeadR.Mag(); thp_r = pLeadR.Theta(); php_r = pLeadR.Phi();

#ifndef LFG_NO_ROOT
    T->Fill();
#else
    fprintf(fout, "%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", weight, region, kMiss_t, kMiss_r,
            pMiss_t, pMiss_r, Q2_t, xB_t, Emiss_t, pp_t, thp_t, pe_t, the_t, Em, kF);
#endif
    nFilled++;
    if (nFilled % 200000 == 0) printf("  %lld events\n", (long long)nFilled);
  }

    printf("Generated %lld events (%lld trials) for %s at E=%.3f GeV with Q2>=%.3f GeV^2 and xB>=%.3f -> %s\n",
      (long long)nFilled, (long long)nTry, target, cfg.Ebeam, cfg.Q2min, cfg.xBmin, outName);
#ifndef LFG_NO_ROOT
  T->Write();
  fout->Close();
#else
  fclose(fout);
#endif
}

#ifdef LFG_NO_ROOT
#include <cstdlib>
int main(int argc, char **argv) {
  long long n = argc > 1 ? atoll(argv[1]) : 200000;
  genLFG_kmiss(n, argc > 2 ? argv[2] : "lfg.csv", argc > 3 ? argv[3] : "C12",
               argc > 4 ? atof(argv[4]) : 5.986,
               argc > 5 ? (unsigned)strtoul(argv[5], nullptr, 10) : 12345,
               argc > 6 ? atof(argv[6]) : 0.0,
               argc > 7 ? atof(argv[7]) : 0.0);
  return 0;
}
#endif