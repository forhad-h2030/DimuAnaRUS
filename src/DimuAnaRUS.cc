#include <fstream>
#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQHitVector.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQDimuonVector.h>
#include <interface_main/SQTrackVector_v1.h>
#include <ktracker/FastTracklet.h>
#include <ktracker/SRecEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <UtilAna/UtilTrack.h>
#include <UtilAna/UtilDimuon.h>
#include <UtilAna/UtilTrigger.h>
#include <UtilAna/TrigRoadset.h>
#include "DimuAnaRUS.h"
using namespace std;

/// DY pT reweight, cf. GetPtReWeight() in e906-gmc-reco/inc/Common.h (target dependence dropped)
double DimuAnaRUS::GetPtReWeight(double pT, double xF, double mass)
{
  const double P = 120.0;
  const double M = 0.938;
  const double S = 2*M*M + 2*M*sqrt(P*P+M*M);
  const double p1_old = 2.8;

  double tau = mass*mass/S;
  double pTsqMax = S/4 * (1-tau)*(1-tau)*(1-xF*xF);
  if (pTsqMax < 0.09) return 1.0;

  const double p1_new = 2.411 - 0.7658 * fabs(xF);

  double f_old = pow(1+pT*pT/p1_old/p1_old, 6) * p1_old*p1_old * (1 - pow(1 + pTsqMax/p1_old/p1_old, -5));
  double f_new = pow(1+pT*pT/p1_new/p1_new, 6) * p1_new*p1_new * (1 - pow(1 + pTsqMax/p1_new/p1_new, -5));
  return f_old / f_new;
}

DimuAnaRUS::DimuAnaRUS(const std::string& name)
  : SubsysReco  (name),
    m_file(0),
    m_tree(0),
    m_tree_name("tree"),
    m_file_name("output.root"),
    m_evt(0),
    m_sp_map(0),
    m_hit_vec(0),
    m_sq_trk_vec(0),
    m_sq_dim_vec(0),
    sqhit_flag(false),
    true_mode(false),
    reco_dimu_mode(false),
    reco_mode(false),
    data_trig_mode(false),
    mc_trig_mode(true)
{
  ;
}

DimuAnaRUS::~DimuAnaRUS()
{
  ;
}

int DimuAnaRUS::Init(PHCompositeNode* startNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int DimuAnaRUS::InitRun(PHCompositeNode* startNode)
{
	m_file = new TFile(m_file_name.c_str(), "RECREATE");

	if (!m_file || m_file->IsZombie()) {
		std::cerr << "Error: Could not create file " << m_file_name << std::endl;
		exit(1);
	} else {
		std::cout << "File " << m_file->GetName() << " opened successfully." << std::endl;
	}

	m_tree = new TTree(m_tree_name.c_str(), "Tree for storing events");
	if (!m_tree) {
		std::cerr << "Error: Could not create tree " << m_tree_name << std::endl;
		exit(1);
	} else {
		std::cout << "Tree " << m_tree->GetName() << " created successfully." << std::endl;
	}

	m_tree->Branch("eventID", &eventID, "eventID/I");
	m_tree->Branch("runID", &runID, "runID/I");
	m_tree->Branch("spillID", &spillID, "spillID/I");
	m_tree->Branch("weight", &weight, "weight/D");
	m_tree->Branch("evtPassTrigEmu", &evtPassTrigEmu, "evtPassTrigEmu/I");

    if(sqhit_flag){
        m_tree->Branch("hitID", &hitID);
        m_tree->Branch("hitTrackID", &hitTrackID);
        m_tree->Branch("processID", &processID);
        m_tree->Branch("detectorID", &detectorID);
        m_tree->Branch("elementID", &elementID);
        m_tree->Branch("tdcTime", &tdcTime);
        m_tree->Branch("driftDistance", &driftDistance);
    }

	if (true_mode) {
        m_tree->Branch("true_dimuon_mass", &true_dimuon_mass, "true_dimuon_mass/D");
        m_tree->Branch("true_dimuon_pT", &true_dimuon_pT, "true_dimuon_pT/D");
        m_tree->Branch("true_dimuon_xF", &true_dimuon_xF, "true_dimuon_xF/D");
        m_tree->Branch("ptReWeight", &ptReWeight, "ptReWeight/D");
        m_tree->Branch("gCharge", &gCharge);
        m_tree->Branch("gTrackID", &gTrackID);
        m_tree->Branch("gvx", &gvx);
        m_tree->Branch("gvy", &gvy);
        m_tree->Branch("gvz", &gvz);
        m_tree->Branch("gpx", &gpx);
        m_tree->Branch("gpy", &gpy);
        m_tree->Branch("gpz", &gpz);

        m_tree->Branch("gx_st1", &gx_st1);
        m_tree->Branch("gy_st1", &gy_st1);
        m_tree->Branch("gz_st1", &gz_st1);
        m_tree->Branch("gpx_st1", &gpx_st1);
        m_tree->Branch("gpy_st1", &gpy_st1);
        m_tree->Branch("gpz_st1", &gpz_st1);

        m_tree->Branch("gx_st3", &gx_st3);
        m_tree->Branch("gy_st3", &gy_st3);
        m_tree->Branch("gz_st3", &gz_st3);
        m_tree->Branch("gpx_st3", &gpx_st3);
        m_tree->Branch("gpy_st3", &gpy_st3);
        m_tree->Branch("gpz_st3", &gpz_st3);
}

	if (reco_mode==true){


        if (reco_dimu_mode == true) {
        m_tree->Branch("rec_dimuon_id", &rec_dimuon_id);
        m_tree->Branch("rec_dimuon_true_id", &rec_dimuon_true_id);
        m_tree->Branch("rec_dimuon_track_id_pos", &rec_dimuon_track_id_pos);
        m_tree->Branch("rec_dimuon_track_id_neg", &rec_dimuon_track_id_neg);

        m_tree->Branch("rec_dimuon_x", &rec_dimuon_x);
        m_tree->Branch("rec_dimuon_y", &rec_dimuon_y);
        m_tree->Branch("rec_dimuon_z", &rec_dimuon_z);

        m_tree->Branch("rec_dimuon_px_pos", &rec_dimuon_px_pos);
        m_tree->Branch("rec_dimuon_py_pos", &rec_dimuon_py_pos);
        m_tree->Branch("rec_dimuon_pz_pos", &rec_dimuon_pz_pos);

        m_tree->Branch("rec_dimuon_px_neg", &rec_dimuon_px_neg);
        m_tree->Branch("rec_dimuon_py_neg", &rec_dimuon_py_neg);
        m_tree->Branch("rec_dimuon_pz_neg", &rec_dimuon_pz_neg);

        m_tree->Branch("rec_dimuon_px_pos_tgt", &rec_dimuon_px_pos_tgt);
        m_tree->Branch("rec_dimuon_py_pos_tgt", &rec_dimuon_py_pos_tgt);
        m_tree->Branch("rec_dimuon_pz_pos_tgt", &rec_dimuon_pz_pos_tgt);

        m_tree->Branch("rec_dimuon_px_neg_tgt", &rec_dimuon_px_neg_tgt);
        m_tree->Branch("rec_dimuon_py_neg_tgt", &rec_dimuon_py_neg_tgt);
        m_tree->Branch("rec_dimuon_pz_neg_tgt", &rec_dimuon_pz_neg_tgt);

        m_tree->Branch("rec_dimuon_mass_tgt", &rec_dimuon_mass_tgt);
        m_tree->Branch("rec_dimuon_pass_trig", &rec_dimuon_pass_trig);

        // track vertex positions of the two muons
        m_tree->Branch("rec_dimuon_x_pos_vtx", &rec_dimuon_x_pos_vtx);
        m_tree->Branch("rec_dimuon_y_pos_vtx", &rec_dimuon_y_pos_vtx);
        m_tree->Branch("rec_dimuon_z_pos_vtx", &rec_dimuon_z_pos_vtx);

        // station-1 x/y of the two muons
        m_tree->Branch("rec_dimuon_x_pos_st1", &rec_dimuon_x_pos_st1);
        m_tree->Branch("rec_dimuon_y_pos_st1", &rec_dimuon_y_pos_st1);
        m_tree->Branch("rec_dimuon_x_neg_st1", &rec_dimuon_x_neg_st1);
        m_tree->Branch("rec_dimuon_y_neg_st1", &rec_dimuon_y_neg_st1);

        m_tree->Branch("rec_dimuon_x_neg_vtx", &rec_dimuon_x_neg_vtx);
        m_tree->Branch("rec_dimuon_y_neg_vtx", &rec_dimuon_y_neg_vtx);
        m_tree->Branch("rec_dimuon_z_neg_vtx", &rec_dimuon_z_neg_vtx);

        }
	}

	// The MC event node gives the event weight; it is optional unless true_mode needs it
	m_evt_true = findNode::getClass<SQMCEvent>(startNode, "SQMCEvent");

	if (true_mode) {
		m_vec_trk = findNode::getClass<SQTrackVector>(startNode, "SQTruthTrackVector");
		m_evt_true = findNode::getClass<SQMCEvent     >(startNode, "SQMCEvent");
		m_true_vec_dim = findNode::getClass<SQDimuonVector>(startNode, "SQTruthDimuonVector");
		if (!m_vec_trk || !m_evt_true) {
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	if (reco_mode) {
		m_sq_trk_vec = findNode::getClass<SQTrackVector>(startNode, "SQRecTrackVector");
		if(reco_dimu_mode==true) m_sq_dim_vec = findNode::getClass<SQDimuonVector>(startNode, "SQRecDimuonVector_PM");

		if (!m_sq_trk_vec) {
			std::cerr << "Error: m_sq_trk_vec is null!" << std::endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}

    if(reco_dimu_mode==true){
        if (!m_sq_dim_vec) {
            std::cerr << "Error: m_sq_dim_vec is null!" << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
    }
	}
	m_evt = findNode::getClass<SQEvent>(startNode, "SQEvent");
    if(sqhit_flag){
        m_hit_vec = findNode::getClass<SQHitVector>(startNode, "SQHitVector");
        trackletVec = findNode::getClass<TrackletVector>(startNode, "TrackletVector");
        if (!m_evt || !m_hit_vec || !trackletVec) {
            return Fun4AllReturnCodes::ABORTEVENT;
        }
   }
	if(reco_mode ==true && data_trig_mode ==true){
		//cout << "inside the data roadset mode: "<<endl;
		SQRun* sq_run = findNode::getClass<SQRun>(startNode, "SQRun");
		if (!sq_run) std::cout << "Error: SQRun  is null!" << std::endl;
		if (!sq_run) return Fun4AllReturnCodes::ABORTEVENT;
		int LBtop = sq_run->get_v1495_id(2);
		int LBbot = sq_run->get_v1495_id(3);
		int ret = m_rs.LoadConfig(LBtop, LBbot);
		if (ret != 0) {
			cout << "!!WARNING!!  OnlMonTrigEP::InitRunOnlMon():  roadset.LoadConfig returned " << ret << ".\n";
		}
		cout <<"Roadset " << m_rs.str(1) << endl;

	}

	if(reco_mode ==true && mc_trig_mode ==true){
		int ret = m_rs.LoadConfig(131);
		if (ret != 0) {
			cout << "!!WARNING!!  OnlMonTrigEP::InitRunOnlMon():  roadset.LoadConfig returned " << ret << ".\n";
		}
		cout <<"MC trigger roadset 131 loaded: " << m_rs.str(1) << endl;
	}
	return Fun4AllReturnCodes::EVENT_OK;
}


int DimuAnaRUS::process_event(PHCompositeNode* startNode)
{

	int proc_id=0;
	int SourceFlag=0;
	int true_proc = m_evt_true ? m_evt_true->get_process_id() : 0;
	//cout <<" true process id "<< true_proc <<endl;
 
	   eventID = m_evt->get_event_id();
	   weight = m_evt_true ? m_evt_true->get_weight() : 1.0;
	   runID = m_evt->get_run_id();
	   spillID = m_evt->get_spill_id();

	   evtPassTrigEmu = 0;


	const double muon_mass = 0.10566;
	TLorentzVector mu1;
	TLorentzVector mu2;
	if(true_mode == true){
		ResetTrueBranches();
		ResetHitBranches();
		if (m_true_vec_dim && m_true_vec_dim->size() > 0) {
			double x1, x2, costh, phi;
			UtilDimuon::CalcVar(m_true_vec_dim->at(0), true_dimuon_mass, true_dimuon_pT, x1, x2, true_dimuon_xF, costh, phi);
			ptReWeight = GetPtReWeight(true_dimuon_pT, true_dimuon_xF, true_dimuon_mass);
		}
		for (unsigned int ii = 0; ii < m_vec_trk->size(); ii++) {

			//cout << "inside the true track loop: "<<endl;
			SQTrack* trk = m_vec_trk->at(ii);
			gCharge.push_back(trk->get_charge());
			gTrackID.push_back(trk->get_track_id());
			gvx.push_back(trk->get_pos_vtx().X());
			gvy.push_back(trk->get_pos_vtx().Y());
			gvz.push_back(trk->get_pos_vtx().Z());
			gpx.push_back(trk->get_mom_vtx().Px());
			gpy.push_back(trk->get_mom_vtx().Py());
			gpz.push_back(trk->get_mom_vtx().Pz());

			gx_st1.push_back(trk->get_pos_st1().X());
			gy_st1.push_back(trk->get_pos_st1().Y());
			gz_st1.push_back(trk->get_pos_st1().Z());
			gpx_st1.push_back(trk->get_mom_st1().Px());
			gpy_st1.push_back(trk->get_mom_st1().Py());
			gpz_st1.push_back(trk->get_mom_st1().Pz());

			gx_st3.push_back(trk->get_pos_st3().X());
			gy_st3.push_back(trk->get_pos_st3().Y());
			gz_st3.push_back(trk->get_pos_st3().Z());
			gpx_st3.push_back(trk->get_mom_st3().Px());
			gpy_st3.push_back(trk->get_mom_st3().Py());
			gpz_st3.push_back(trk->get_mom_st3().Pz());

            //cout << "true track id: "<< trk->get_track_id() << " charge: "<<trk->get_charge()<<  endl;

			/*
			// Only calculate dimuon if we have exactly 2 tracks AND the true_proc is between 400 and 412
			if (m_vec_trk->size() == 2 && true_proc >= 400 && true_proc <= 412) {
				if (ii == 0) {
					mu1.SetXYZM(trk->get_mom_vtx().Px(),
							trk->get_mom_vtx().Py(),
							trk->get_mom_vtx().Pz(),
							muon_mass);
				}

				if (ii == 1) {
					mu2.SetXYZM(trk->get_mom_vtx().Px(),
							trk->get_mom_vtx().Py(),
							trk->get_mom_vtx().Pz(),
							muon_mass);

					TLorentzVector dimuon = mu1 + mu2;
					double dimuon_mass = dimuon.M();
					std::cout << "Dimuon invariant mass: " << dimuon_mass << " GeV" << std::endl;
					std::cout << "true_proc " << true_proc << std::endl;

					// Assign proc_id by mass
					if (std::abs(dimuon_mass - 3.097) < 0.05) {
						proc_id = 12;  // J/ψ → μ+μ−
					}
					else if (std::abs(dimuon_mass - 3.686) < 0.05) {
						proc_id = 13;  // ψ′ → μ+μ−
					}
				}
			}			
			*/

			//if (true_proc == 223)proc_id = 11;
			//if (true_proc == 0)proc_id = 14;

			int SourceFlag=0;
			double z = trk->get_pos_vtx().Z();
			if (z <= -296. && z >= -304.) SourceFlag = 1; // target
			else if (z >= 0. && z < 250) SourceFlag = 2;   // dump
			else if (z > -296. && z < 0.) SourceFlag = 3;  // gap
			else if (z < -304) SourceFlag = 0;              // upstream

            if(sqhit_flag){
                if (m_hit_vec) {
                    for (int ihit = 0; ihit < m_hit_vec->size(); ++ihit) {
                        SQHit* hit = m_hit_vec->at(ihit);
                        //if (!hit->is_in_time()) continue;
                        proc_id=12;
                        if(hit->get_track_id() != trk->get_track_id()) continue;
                        int processID_;
                        if(trk->get_charge() >0) processID_ =proc_id;
                        else processID_ = proc_id +10 ;
                        int sourceFlag_= SourceFlag;
                        unsigned int encodedValue = EncodeProcess(processID_, sourceFlag_);
                        //cout << "charge: "<< trk->get_charge() <<endl;
                        //cout << "encodedValue: "<< encodedValue<< endl;
                        //cout << "DecodeSourceFlag: "<< DecodeSourceFlag(encodedValue) << endl;
                        //cout << "DecodeProcessID: "<< DecodeProcessID(encodedValue) << endl;
                        //cout << "det Name: "<< hit->get_detector_id() << " Hit Id: "<< hit->get_hit_id () << endl;
                        processID.push_back(encodedValue);
                        hitID.push_back(hit->get_hit_id());
                        hitTrackID.push_back(hit->get_track_id());
                        detectorID.push_back(hit->get_detector_id());
                        elementID.push_back(hit->get_element_id());
                        tdcTime.push_back(hit->get_tdc_time());
                        driftDistance.push_back(hit->get_drift_distance());
                    }   
                }
            }
		}
	}



	//cout << "reco mode: "<< reco_mode << endl;
	//std::cout << "Size of m_sq_dim_vec: " << m_sq_dim_vec->size() << std::endl;

	if(reco_mode == true){
		if(reco_dimu_mode==true){
			ResetRecoDimuBranches();
			for (auto it = m_sq_dim_vec->begin(); it != m_sq_dim_vec->end(); it++) {
				SRecDimuon& sdim = dynamic_cast<SRecDimuon&>(**it);
				int trk_id_pos = sdim.get_track_id_pos();
				int trk_id_neg = sdim.get_track_id_neg();
				SRecTrack& trk_pos = dynamic_cast<SRecTrack&>(*(m_sq_trk_vec->at(trk_id_pos)));
				SRecTrack& trk_neg = dynamic_cast<SRecTrack&>(*(m_sq_trk_vec->at(trk_id_neg)));

				// --- Cut: z vertex > -690 cm (both track vertices) ---
				if (trk_pos.get_pos_vtx().Z() < -690.) continue;
				if (trk_neg.get_pos_vtx().Z() < -690.) continue;

				// --- Cut: |y_st1| > 3 cm (both tracks) ---
				if (fabs(trk_pos.get_pos_st1().Y()) < 3.) continue;
				if (fabs(trk_neg.get_pos_st1().Y()) < 3.) continue;

				// --- Cut: chi2 — all hypotheses valid, target must be the best ---
				double chisq_target_pos   = trk_pos.getChisqTarget();
				double chisq_dump_pos     = trk_pos.get_chisq_dump();
				double chisq_upstream_pos = trk_pos.get_chisq_upstream();
				double chisq_target_neg   = trk_neg.getChisqTarget();
				double chisq_dump_neg     = trk_neg.get_chisq_dump();
				double chisq_upstream_neg = trk_neg.get_chisq_upstream();
				bool pass_chisq_cut =
					(chisq_target_pos   >= 0 &&
					 chisq_dump_pos     >= 0 &&
					 chisq_upstream_pos >= 0 &&
					 (chisq_target_pos - chisq_dump_pos)     <= 0 &&
					 (chisq_target_pos - chisq_upstream_pos) <= 0 &&
					 chisq_target_neg   >= 0 &&
					 chisq_dump_neg     >= 0 &&
					 chisq_upstream_neg >= 0 &&
					 (chisq_target_neg - chisq_dump_neg)     <= 0 &&
					 (chisq_target_neg - chisq_upstream_neg) <= 0);
				if (!pass_chisq_cut) continue;

				// --- Target hypothesis & mass > 0 ---
				sdim.calcVariables(1); // 1 = target
				if ((sdim.p_pos_target + sdim.p_neg_target).M() <= 0.) continue;

				// Trigger emulation: matched roads in opposite top/bottom halves (flag only, no rejection)
				bool dimu_pass_trig = true; // no emulation requested -> nothing to reject
				if (reco_mode && (data_trig_mode || mc_trig_mode)) {
					std::vector<int> list_road_pos = UtilTrack::FindMatchedRoads(&trk_pos);
					std::vector<int> list_road_neg = UtilTrack::FindMatchedRoads(&trk_neg);
					std::vector<int> list_road_pt = m_rs.PosTop()->FindRoadIDs(list_road_pos);
					std::vector<int> list_road_pb = m_rs.PosBot()->FindRoadIDs(list_road_pos);
					std::vector<int> list_road_nt = m_rs.NegTop()->FindRoadIDs(list_road_neg);
					std::vector<int> list_road_nb = m_rs.NegBot()->FindRoadIDs(list_road_neg);
					bool pass_trig = (list_road_pt.size()>0 && list_road_nb.size()>0) ||
					                 (list_road_pb.size()>0 && list_road_nt.size()>0);
					static int n_trig_print = 0;
					if (n_trig_print < 20) { // sanity check on the first dimuons only
						++n_trig_print;
						cout << "TrigEmu: evt " << m_evt->get_event_id()
						     << " matched roads pos/neg = " << list_road_pos.size() << "/" << list_road_neg.size()
						     << "  PT/PB/NT/NB = " << list_road_pt.size() << "/" << list_road_pb.size()
						     << "/" << list_road_nt.size() << "/" << list_road_nb.size()
						     << "  pass = " << pass_trig << endl;
					}
					dimu_pass_trig = pass_trig; // flag only; selection is applied offline
				}
				if (dimu_pass_trig) evtPassTrigEmu = 1;


				//--------
				//vtx
				rec_dimuon_px_pos_vtx.push_back(trk_pos.get_mom_vtx().Px());
				rec_dimuon_py_pos_vtx.push_back(trk_pos.get_mom_vtx().Py());
				rec_dimuon_pz_pos_vtx.push_back(trk_pos.get_mom_vtx().Pz());

				rec_dimuon_x_pos_vtx.push_back(trk_pos.get_pos_vtx().X());
				rec_dimuon_y_pos_vtx.push_back(trk_pos.get_pos_vtx().Y());
				rec_dimuon_z_pos_vtx.push_back(trk_pos.get_pos_vtx().Z());

				rec_dimuon_px_neg_vtx.push_back(trk_neg.get_mom_vtx().Px());
				rec_dimuon_py_neg_vtx.push_back(trk_neg.get_mom_vtx().Py());
				rec_dimuon_pz_neg_vtx.push_back(trk_neg.get_mom_vtx().Pz());

				rec_dimuon_x_neg_vtx.push_back(trk_neg.get_pos_vtx().X());
				rec_dimuon_y_neg_vtx.push_back(trk_neg.get_pos_vtx().Y());
				rec_dimuon_z_neg_vtx.push_back(trk_neg.get_pos_vtx().Z());


				// Station 1 - positive trk
				rec_dimuon_px_pos_st1.push_back(trk_pos.get_mom_st1().Px());
				rec_dimuon_py_pos_st1.push_back(trk_pos.get_mom_st1().Py());
				rec_dimuon_pz_pos_st1.push_back(trk_pos.get_mom_st1().Pz());

				rec_dimuon_x_pos_st1.push_back(trk_pos.get_pos_st1().X());
				rec_dimuon_y_pos_st1.push_back(trk_pos.get_pos_st1().Y());
				rec_dimuon_z_pos_st1.push_back(trk_pos.get_pos_st1().Z());

				// Station 3 - positive trk
				rec_dimuon_px_pos_st3.push_back(trk_pos.get_mom_st3().Px());
				rec_dimuon_py_pos_st3.push_back(trk_pos.get_mom_st3().Py());
				rec_dimuon_pz_pos_st3.push_back(trk_pos.get_mom_st3().Pz());

				rec_dimuon_x_pos_st3.push_back(trk_pos.get_pos_st3().X());
				rec_dimuon_y_pos_st3.push_back(trk_pos.get_pos_st3().Y());
				rec_dimuon_z_pos_st3.push_back(trk_pos.get_pos_st3().Z());


				// Station 1 - negative trk
				rec_dimuon_px_neg_st1.push_back(trk_neg.get_mom_st1().Px());
				rec_dimuon_py_neg_st1.push_back(trk_neg.get_mom_st1().Py());
				rec_dimuon_pz_neg_st1.push_back(trk_neg.get_mom_st1().Pz());

				rec_dimuon_x_neg_st1.push_back(trk_neg.get_pos_st1().X());
				rec_dimuon_y_neg_st1.push_back(trk_neg.get_pos_st1().Y());
				rec_dimuon_z_neg_st1.push_back(trk_neg.get_pos_st1().Z());

				// Station 3 - negative trk
				rec_dimuon_px_neg_st3.push_back(trk_neg.get_mom_st3().Px());
				rec_dimuon_py_neg_st3.push_back(trk_neg.get_mom_st3().Py());
				rec_dimuon_pz_neg_st3.push_back(trk_neg.get_mom_st3().Pz());

				rec_dimuon_x_neg_st3.push_back(trk_neg.get_pos_st3().X());
				rec_dimuon_y_neg_st3.push_back(trk_neg.get_pos_st3().Y());
				rec_dimuon_z_neg_st3.push_back(trk_neg.get_pos_st3().Z());


				rec_dimuon_id.push_back(sdim.get_dimuon_id());
				rec_dimuon_true_id.push_back(sdim.get_rec_dimuon_id());
				rec_dimuon_track_id_pos.push_back(trk_id_pos);
				rec_dimuon_track_id_neg.push_back(trk_id_neg);
				rec_dimuon_px_pos.push_back(sdim.get_mom_pos().Px());
				rec_dimuon_py_pos.push_back(sdim.get_mom_pos().Py());
				rec_dimuon_pz_pos.push_back(sdim.get_mom_pos().Pz());
				rec_dimuon_px_neg.push_back(sdim.get_mom_neg().Px());
				rec_dimuon_py_neg.push_back(sdim.get_mom_neg().Py());
				rec_dimuon_pz_neg.push_back(sdim.get_mom_neg().Pz());
				rec_dimuon_x.push_back(sdim.get_pos().X());
				rec_dimuon_y.push_back(sdim.get_pos().Y());
				rec_dimuon_z.push_back(sdim.get_pos().Z());
				// ===== Target hypothesis =====
				sdim.calcVariables(1); // 1 = target
				rec_dimuon_px_pos_tgt.push_back(sdim.p_pos_target.Px());
				rec_dimuon_py_pos_tgt.push_back(sdim.p_pos_target.Py());
				rec_dimuon_pz_pos_tgt.push_back(sdim.p_pos_target.Pz());
				rec_dimuon_px_neg_tgt.push_back(sdim.p_neg_target.Px());
				rec_dimuon_py_neg_tgt.push_back(sdim.p_neg_target.Py());
				rec_dimuon_pz_neg_tgt.push_back(sdim.p_neg_target.Pz());

				TLorentzVector mom_tgt = sdim.p_pos_target + sdim.p_neg_target;
				rec_dimuon_mass_tgt.push_back(mom_tgt.M());
				rec_dimuon_pass_trig.push_back(dimu_pass_trig ? 1 : 0);
			}
		}
	}
	m_tree->Fill();
	return Fun4AllReturnCodes::EVENT_OK;
}

int DimuAnaRUS::End(PHCompositeNode* startNode)
{
  m_file->cd();
  m_file->Write();
  m_file->Close();  
  return Fun4AllReturnCodes::EVENT_OK;
}


void DimuAnaRUS::ResetHitBranches() {
 	hitID.clear();
	processID.clear();
 	hitTrackID.clear();
	detectorID.clear();
	elementID.clear();
	tdcTime.clear();
	driftDistance.clear();
}
void DimuAnaRUS::ResetTrueBranches() {
	true_dimuon_mass = true_dimuon_pT = true_dimuon_xF = -999;
	ptReWeight = 1.0;
	gCharge.clear(); gTrackID.clear();
	gvx.clear(); gvy.clear(); gvz.clear();
	gpx.clear(); gpy.clear(); gpz.clear();
	gx_st1.clear(); gy_st1.clear(); gz_st1.clear();
	gpx_st1.clear(); gpy_st1.clear(); gpz_st1.clear();
	gx_st3.clear(); gy_st3.clear(); gz_st3.clear();
	gpx_st3.clear(); gpy_st3.clear(); gpz_st3.clear();
}

void DimuAnaRUS::ResetRecoDimuBranches() {
    rec_dimuon_id.clear(); rec_dimuon_true_id.clear(); rec_dimuon_track_id_pos.clear(); rec_dimuon_track_id_neg.clear();
    rec_dimuon_x.clear(); rec_dimuon_y.clear(); rec_dimuon_z.clear();
    rec_dimuon_px_pos.clear(); rec_dimuon_py_pos.clear(); rec_dimuon_pz_pos.clear();
    rec_dimuon_px_neg.clear(); rec_dimuon_py_neg.clear(); rec_dimuon_pz_neg.clear();
    rec_dimuon_px_pos_tgt.clear(); rec_dimuon_py_pos_tgt.clear(); rec_dimuon_pz_pos_tgt.clear();
    rec_dimuon_px_neg_tgt.clear(); rec_dimuon_py_neg_tgt.clear(); rec_dimuon_pz_neg_tgt.clear();
    rec_dimuon_mass_tgt.clear();
    rec_dimuon_pass_trig.clear();

    rec_dimuon_x_pos_vtx.clear();
    rec_dimuon_y_pos_vtx.clear();
    rec_dimuon_z_pos_vtx.clear();
    rec_dimuon_px_pos_vtx.clear();
    rec_dimuon_py_pos_vtx.clear();
    rec_dimuon_pz_pos_vtx.clear();

    rec_dimuon_x_neg_vtx.clear();
    rec_dimuon_y_neg_vtx.clear();
    rec_dimuon_z_neg_vtx.clear();
    rec_dimuon_px_neg_vtx.clear();
    rec_dimuon_py_neg_vtx.clear();
    rec_dimuon_pz_neg_vtx.clear();

    rec_dimuon_x_pos_st1.clear();
    rec_dimuon_y_pos_st1.clear();
    rec_dimuon_z_pos_st1.clear();
    rec_dimuon_x_neg_st1.clear();
    rec_dimuon_y_neg_st1.clear();
    rec_dimuon_z_neg_st1.clear();
    rec_dimuon_x_pos_st3.clear();
    rec_dimuon_y_pos_st3.clear();
    rec_dimuon_z_pos_st3.clear();
    rec_dimuon_x_neg_st3.clear();
    rec_dimuon_y_neg_st3.clear();
    rec_dimuon_z_neg_st3.clear();
    rec_dimuon_px_pos_st1.clear();
    rec_dimuon_py_pos_st1.clear();
    rec_dimuon_pz_pos_st1.clear();
    rec_dimuon_px_neg_st1.clear();
    rec_dimuon_py_neg_st1.clear();
    rec_dimuon_pz_neg_st1.clear();
    rec_dimuon_px_pos_st3.clear();
    rec_dimuon_py_pos_st3.clear();
    rec_dimuon_pz_pos_st3.clear();
    rec_dimuon_px_neg_st3.clear();
    rec_dimuon_py_neg_st3.clear();
    rec_dimuon_pz_neg_st3.clear();
}

unsigned int DimuAnaRUS::EncodeProcess(int processID, int sourceFlag) {
        return ( (sourceFlag & 0x3) << 5 ) |  // 2 bits for sourceFlag (1-3)
	                (processID & 0x1F);            // 5 bits for processID (11-24)
			               }
int DimuAnaRUS::DecodeSourceFlag(unsigned int encoded) { 
	return (encoded >> 5) & 0x3;
}
int DimuAnaRUS::DecodeProcessID(unsigned int encoded) {
	return encoded & 0x1F;
}
