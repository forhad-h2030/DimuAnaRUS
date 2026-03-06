#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TLorentzVector.h>
#include <vector>
#include <string>

constexpr double MUON_MASS = 0.1056; // GeV/c^2
template <typename T>
T GetByCharge(const std::vector<T> &values,
              const std::vector<int> &charges,
              bool wantPositive)
{
    int idx_pos = (charges[0] > 0) ? 0 : 1;
    int idx_neg = (charges[0] > 0) ? 1 : 0;
    int idx = wantPositive ? idx_pos : idx_neg;
    return (idx < static_cast<int>(values.size())) ? values[idx] : T{};
}

void write() {
    std::string infile = "../MC_Non_Jpsi_Oct23_Target.root";
    std::string treename = "tree";

    ROOT::RDataFrame df(treename, infile);

    auto df2 = df
        // muon momentum components from your tree
        .Define("mu_pos_px", "rec_dimuon_px_pos_tgt[0]")
        .Define("mu_pos_py", "rec_dimuon_py_pos_tgt[0]")
        .Define("mu_pos_pz", "rec_dimuon_pz_pos_tgt[0]")
        .Define("mu_neg_px", "rec_dimuon_px_neg_tgt[0]")
        .Define("mu_neg_py", "rec_dimuon_py_neg_tgt[0]")
        .Define("mu_neg_pz", "rec_dimuon_pz_neg_tgt[0]")

        // track info with charge filtering
        .Define("rec_track_pos_x_st1",
                [](const std::vector<double>& rec_track_vx,
                   const std::vector<int>& rec_track_charge) {
                    return GetByCharge(rec_track_vx, rec_track_charge, true);
                }, {"rec_track_vx","rec_track_charge"})
        .Define("rec_track_neg_x_st1",
                [](const std::vector<double>& rec_track_vx,
                   const std::vector<int>& rec_track_charge) {
                    return GetByCharge(rec_track_vx, rec_track_charge, false);
                }, {"rec_track_vx","rec_track_charge"})
        .Define("rec_track_pos_px_st1",
                [](const std::vector<double>& rec_track_px,
                   const std::vector<int>& rec_track_charge) {
                    return GetByCharge(rec_track_px, rec_track_charge, true);
                }, {"rec_track_px","rec_track_charge"})
        .Define("rec_track_neg_px_st1",
                [](const std::vector<double>& rec_track_px,
                   const std::vector<int>& rec_track_charge) {
                    return GetByCharge(rec_track_px, rec_track_charge, false);
                }, {"rec_track_px","rec_track_charge"})

        // invariant mass of the two muons
        .Define("mass",
                [](double px1, double py1, double pz1,
                   double px2, double py2, double pz2) {
                    TLorentzVector track_pos, track_neg;
                    track_pos.SetXYZM(px1, py1, pz1, MUON_MASS);
                    track_neg.SetXYZM(px2, py2, pz2, MUON_MASS);
                    return (track_pos + track_neg).M();
                },
                {"mu_pos_px","mu_pos_py","mu_pos_pz",
                 "mu_neg_px","mu_neg_py","mu_neg_pz"});

    std::vector<std::string> outBranches = {
        "mu_pos_px","mu_pos_py","mu_pos_pz",
        "mu_neg_px","mu_neg_py","mu_neg_pz",
        "rec_track_pos_x_st1","rec_track_neg_x_st1",
        "rec_track_pos_px_st1","rec_track_neg_px_st1",
        "mass" // include the new mass branch
    };

    auto df3 = df2.Filter("mass>1.0 && mass < 8.0");

    df3.Snapshot("tree", "comb_target_aug19.root", outBranches);
}

