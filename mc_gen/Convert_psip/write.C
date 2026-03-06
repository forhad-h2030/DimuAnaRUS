#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TLorentzVector.h>
#include <vector>
#include <string>
#include <algorithm>

double MUON_MASS = 0.1056; // GeV/c^2

template <typename T>
T FirstByCharge(const std::vector<T>& values,
                const std::vector<int>& charges,
                bool wantPositive,
                T def = T{})
{
  const int n = std::min(values.size(), charges.size());
  for (int i = 0; i < n; ++i) {
    if (wantPositive && charges[i] > 0) return values[i];
    if (!wantPositive && charges[i] < 0) return values[i];
  }
  return def;
}

void write() {
  //std::string infile   = "../MC_DY_Pythia8_Target_July18.root";
  //std::string infile   = "../MC_DY_Pythia8_Mass1_8_legacy.root";
  //std::string infile   = "../MC_DY_Pythia8_Mass0.5_2.5.root";
  std::string infile   = "../MC_Psip_Pythia8_Target_Feb18.root";
  std::string treename = "tree";

  ROOT::RDataFrame df(treename, infile);

  auto df2 = df
    // muon momentum components from your tree (guard empty vectors if needed)
    .Define("rec_dimu_mu_pos_px", "rec_dimuon_px_pos_tgt.size()>0 ? rec_dimuon_px_pos_tgt[0] : 0.0")
    .Define("rec_dimu_mu_pos_py", "rec_dimuon_py_pos_tgt.size()>0 ? rec_dimuon_py_pos_tgt[0] : 0.0")
    .Define("rec_dimu_mu_pos_pz", "rec_dimuon_pz_pos_tgt.size()>0 ? rec_dimuon_pz_pos_tgt[0] : 0.0")
    .Define("rec_dimu_mu_neg_px", "rec_dimuon_px_neg_tgt.size()>0 ? rec_dimuon_px_neg_tgt[0] : 0.0")
    .Define("rec_dimu_mu_neg_py", "rec_dimuon_py_neg_tgt.size()>0 ? rec_dimuon_py_neg_tgt[0] : 0.0")
    .Define("rec_dimu_mu_neg_pz", "rec_dimuon_pz_neg_tgt.size()>0 ? rec_dimuon_pz_neg_tgt[0] : 0.0")

    // track info from ST1 using charge selection
    .Define("rec_track_pos_x_st1",
            [](const std::vector<double>& x_st1, const std::vector<int>& q) {
              return FirstByCharge(x_st1, q, true, 0.0);
            }, {"rec_track_x_st1", "rec_track_charge"})
    .Define("rec_track_neg_x_st1",
            [](const std::vector<double>& x_st1, const std::vector<int>& q) {
              return FirstByCharge(x_st1, q, false, 0.0);
            }, {"rec_track_x_st1", "rec_track_charge"})

    .Define("rec_track_pos_px_st1",
            [](const std::vector<double>& px_st1, const std::vector<int>& q) {
              return FirstByCharge(px_st1, q, true, 0.0);
            }, {"rec_track_px_st1", "rec_track_charge"})
    .Define("rec_track_neg_px_st1",
            [](const std::vector<double>& px_st1, const std::vector<int>& q) {
              return FirstByCharge(px_st1, q, false, 0.0);
            }, {"rec_track_px_st1", "rec_track_charge"})
    .Define("rec_track_pos_vx",
            [](const std::vector<double>& x, const std::vector<int>& q) {
              return FirstByCharge(x, q, true, 0.0);
            }, {"rec_track_vx", "rec_track_charge"})
    .Define("rec_track_pos_vy",
            [](const std::vector<double>& y, const std::vector<int>& q) {
              return FirstByCharge(y, q, true, 0.0);
            }, {"rec_track_vy", "rec_track_charge"})
    .Define("rec_track_pos_vz",
            [](const std::vector<double>& z, const std::vector<int>& q) {
              return FirstByCharge(z, q, true, 0.0);
            }, {"rec_track_vz", "rec_track_charge"})

    .Define("rec_track_neg_vx",
            [](const std::vector<double>& x, const std::vector<int>& q) {
              return FirstByCharge(x, q, false, 0.0);
            }, {"rec_track_vx", "rec_track_charge"})
    .Define("rec_track_neg_vy",
            [](const std::vector<double>& y, const std::vector<int>& q) {
              return FirstByCharge(y, q, false, 0.0);
            }, {"rec_track_vy", "rec_track_charge"})
    .Define("rec_track_neg_vz",
            [](const std::vector<double>& z, const std::vector<int>& q) {
              return FirstByCharge(z, q, false, 0.0);
            }, {"rec_track_vz", "rec_track_charge"})

    // invariant mass of the two muons
    .Define("mass",
            [](double px1, double py1, double pz1,
               double px2, double py2, double pz2) {
              TLorentzVector mu_pos, mu_neg;
              mu_pos.SetXYZM(px1, py1, pz1, MUON_MASS);
              mu_neg.SetXYZM(px2, py2, pz2, MUON_MASS);
              return (mu_pos + mu_neg).M();
            },
            {"rec_dimu_mu_pos_px","rec_dimu_mu_pos_py","rec_dimu_mu_pos_pz",
             "rec_dimu_mu_neg_px","rec_dimu_mu_neg_py","rec_dimu_mu_neg_pz"});

  std::vector<std::string> outBranches = {
    "rec_dimu_mu_pos_px","rec_dimu_mu_pos_py","rec_dimu_mu_pos_pz",
    "rec_dimu_mu_neg_px","rec_dimu_mu_neg_py","rec_dimu_mu_neg_pz",
    "rec_track_pos_x_st1","rec_track_neg_x_st1",
    "rec_track_pos_px_st1","rec_track_neg_px_st1",
    "rec_track_pos_vx","rec_track_pos_vy","rec_track_pos_vz",
    "rec_track_neg_vx","rec_track_neg_vy","rec_track_neg_vz",
    "mass"
  };

  auto df3 = df2.Filter("mass > 0.0 && mass < 8.0");
  df3.Snapshot("tree", "MC_Psip_Pythia8_Target_Feb18.root", outBranches);
}

