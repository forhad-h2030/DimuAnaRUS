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
#include <geom_svc/GeomSvc.h>
#include "DimuAnaRUS.h"
using namespace std;

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
    saveDimuonOnly(false),
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
	m_tree->Branch("eventID", &eventID, "eventID/I");
	m_tree->Branch("rfID", &rfID, "rfID/I");
	m_tree->Branch("turnID", &turnID, "turnID/I");
	m_tree->Branch("rfIntensity", rfIntensity, "rfIntensity[33]/I");
	m_tree->Branch("fpgaTrigger", fpgaTrigger, "fpgaTrigger[5]/I");
	m_tree->Branch("nimTrigger", nimTrigger, "nimTrigger[5]/I");

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
        m_tree->Branch("rec_track_id",             &rec_track_id);
        m_tree->Branch("rec_track_charge",         &rec_track_charge);
        m_tree->Branch("rec_track_vx",             &rec_track_vx);
        m_tree->Branch("rec_track_vy",             &rec_track_vy);
        m_tree->Branch("rec_track_vz",             &rec_track_vz);
        m_tree->Branch("rec_track_px",             &rec_track_px);
        m_tree->Branch("rec_track_py",             &rec_track_py);
        m_tree->Branch("rec_track_pz",             &rec_track_pz);
        m_tree->Branch("rec_track_x_st1",          &rec_track_x_st1);
        m_tree->Branch("rec_track_y_st1",          &rec_track_y_st1);
        m_tree->Branch("rec_track_z_st1",          &rec_track_z_st1);
        m_tree->Branch("rec_track_px_st1",         &rec_track_px_st1);
        m_tree->Branch("rec_track_py_st1",         &rec_track_py_st1);
        m_tree->Branch("rec_track_pz_st1",         &rec_track_pz_st1);
        m_tree->Branch("rec_track_x_st3",          &rec_track_x_st3);
        m_tree->Branch("rec_track_y_st3",          &rec_track_y_st3);
        m_tree->Branch("rec_track_z_st3",          &rec_track_z_st3);
        m_tree->Branch("rec_track_px_st3",         &rec_track_px_st3);
        m_tree->Branch("rec_track_py_st3",         &rec_track_py_st3);
        m_tree->Branch("rec_track_pz_st3",         &rec_track_pz_st3);
        m_tree->Branch("rec_track_num_hits",       &rec_track_num_hits);
        m_tree->Branch("rec_track_chisq",          &rec_track_chisq);
        m_tree->Branch("rec_track_chisq_tgt",   &rec_track_chisq_tgt);
        m_tree->Branch("rec_track_chisq_dump",     &rec_track_chisq_dump);
        m_tree->Branch("rec_track_chisq_upstream", &rec_track_chisq_upstream);
        m_tree->Branch("rec_track_x_tgt",          &rec_track_x_tgt);
        m_tree->Branch("rec_track_y_tgt",          &rec_track_y_tgt);
        m_tree->Branch("rec_track_z_tgt",          &rec_track_z_tgt);
        m_tree->Branch("rec_track_px_tgt",         &rec_track_px_tgt);
        m_tree->Branch("rec_track_py_tgt",         &rec_track_py_tgt);
        m_tree->Branch("rec_track_pz_tgt",         &rec_track_pz_tgt);
        m_tree->Branch("rec_track_x_dump",         &rec_track_x_dump);
        m_tree->Branch("rec_track_y_dump",         &rec_track_y_dump);
        m_tree->Branch("rec_track_z_dump",         &rec_track_z_dump);
        m_tree->Branch("rec_track_px_dump",        &rec_track_px_dump);
        m_tree->Branch("rec_track_py_dump",        &rec_track_py_dump);
        m_tree->Branch("rec_track_pz_dump",        &rec_track_pz_dump);
        m_tree->Branch("rec_hit_ids",              &rec_hit_ids);
        m_tree->Branch("rec_track_hit_x",              &rec_track_hit_x);
        m_tree->Branch("rec_track_hit_y",              &rec_track_hit_y);


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

        m_tree->Branch("rec_dimuon_px_pos_dump", &rec_dimuon_px_pos_dump);
        m_tree->Branch("rec_dimuon_py_pos_dump", &rec_dimuon_py_pos_dump);
        m_tree->Branch("rec_dimuon_pz_pos_dump", &rec_dimuon_pz_pos_dump);

        m_tree->Branch("rec_dimuon_px_neg_dump", &rec_dimuon_px_neg_dump);
        m_tree->Branch("rec_dimuon_py_neg_dump", &rec_dimuon_py_neg_dump);
        m_tree->Branch("rec_dimuon_pz_neg_dump", &rec_dimuon_pz_neg_dump);


        m_tree->Branch("rec_dimuon_mass_tgt", &rec_dimuon_mass_tgt);

        }
	}

	if (true_mode) {
		m_vec_trk = findNode::getClass<SQTrackVector>(startNode, "SQTruthTrackVector");
		m_evt_true = findNode::getClass<SQMCEvent     >(startNode, "SQMCEvent");
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
		//cout <<"Roadset " << m_rs.str(1) << endl;
	}
	return Fun4AllReturnCodes::EVENT_OK;
}


int DimuAnaRUS::process_event(PHCompositeNode* startNode)
{

	int proc_id=0;
	int SourceFlag=0;
	int true_proc = m_evt_true->get_process_id();
	//cout <<" true process id "<< true_proc <<endl;
	/*   //if (! m_evt->get_trigger(SQEvent::MATRIX1)) {
	   if (! m_evt->get_trigger(SQEvent::NIM4)) {
	   return Fun4AllReturnCodes::EVENT_OK;
	   }
	*/
	//if (saveDimuonOnly && m_sq_dim_vec->empty()) {
	//	return 0;  // Skip this event if no dimuons are present and the mode is enabled
	//}
 
	   eventID = m_evt->get_event_id();
	   runID = m_evt->get_run_id();
	   spillID = m_evt->get_spill_id();
	   rfID = m_evt->get_qie_rf_id();
	   turnID = m_evt->get_qie_turn_id();

	   fpgaTrigger[0] = m_evt->get_trigger(SQEvent::MATRIX1);
	   fpgaTrigger[1] = m_evt->get_trigger(SQEvent::MATRIX2);
	   fpgaTrigger[2] = m_evt->get_trigger(SQEvent::MATRIX3);
	   fpgaTrigger[3] = m_evt->get_trigger(SQEvent::MATRIX4);
	   fpgaTrigger[4] = m_evt->get_trigger(SQEvent::MATRIX5);

	   nimTrigger[0] = m_evt->get_trigger(SQEvent::NIM1);
	   nimTrigger[1] = m_evt->get_trigger(SQEvent::NIM2);
	   nimTrigger[2] = m_evt->get_trigger(SQEvent::NIM3);
	   nimTrigger[3] = m_evt->get_trigger(SQEvent::NIM4);
	   nimTrigger[4] = m_evt->get_trigger(SQEvent::NIM5);
	   for (int i = -16; i <= 16; ++i) {
	   rfIntensity[i + 16] = m_evt->get_qie_rf_intensity(i);
	   }

	const double muon_mass = 0.10566;
	TLorentzVector mu1;
	TLorentzVector mu2;
	if(true_mode == true){
		ResetTrueBranches();
		ResetHitBranches();
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

int index = -1;

	if(reco_mode == true){
		ResetRecoBranches();

		for (auto it = m_sq_trk_vec->begin(); it != m_sq_trk_vec->end(); ++it) {
			index+=1;
			SRecTrack* trk = dynamic_cast<SRecTrack*>(*it);

			//cout << "track id: "<< trk->get_rec_track_id ()  <<" charge: " << trk->get_charge() << "index: "<<endl;

			rec_track_id.push_back(index);  // where i_trk is the index in the reco track loop
			// Basic track info
			rec_track_charge.push_back(trk->get_charge());
			rec_track_vx.push_back(trk->get_pos_vtx().X());
			rec_track_vy.push_back(trk->get_pos_vtx().Y());
			rec_track_vz.push_back(trk->get_pos_vtx().Z());
			rec_track_px.push_back(trk->get_mom_vtx().Px());
			rec_track_py.push_back(trk->get_mom_vtx().Py());
			rec_track_pz.push_back(trk->get_mom_vtx().Pz());

			// Station 1
			rec_track_x_st1.push_back(trk->get_pos_st1().X());
			rec_track_y_st1.push_back(trk->get_pos_st1().Y());
			rec_track_z_st1.push_back(trk->get_pos_st1().Z());
			rec_track_px_st1.push_back(trk->get_mom_st1().Px());
			rec_track_py_st1.push_back(trk->get_mom_st1().Py());
			rec_track_pz_st1.push_back(trk->get_mom_st1().Pz());

			// Station 3
			rec_track_x_st3.push_back(trk->get_pos_st3().X());
			rec_track_y_st3.push_back(trk->get_pos_st3().Y());
			rec_track_z_st3.push_back(trk->get_pos_st3().Z());
			rec_track_px_st3.push_back(trk->get_mom_st3().Px());
			rec_track_py_st3.push_back(trk->get_mom_st3().Py());
			rec_track_pz_st3.push_back(trk->get_mom_st3().Pz());

			// Chi-squared
			rec_track_chisq.push_back(trk->get_chisq());
			rec_track_chisq_tgt.push_back(trk->getChisqTarget());
			rec_track_chisq_dump.push_back(trk->get_chisq_dump());
			rec_track_chisq_upstream.push_back(trk->get_chisq_upstream());

			// Number of hits
			rec_track_num_hits.push_back(trk->get_num_hits());

			// Target region
			rec_track_x_tgt.push_back(trk->get_pos_target().X());
			rec_track_y_tgt.push_back(trk->get_pos_target().Y());
			rec_track_z_tgt.push_back(trk->get_pos_target().Z());
			rec_track_px_tgt.push_back(trk->get_mom_target().Px());
			rec_track_py_tgt.push_back(trk->get_mom_target().Py());
			rec_track_pz_tgt.push_back(trk->get_mom_target().Pz());

			// Dump region
			rec_track_x_dump.push_back(trk->get_pos_dump().X());
			rec_track_y_dump.push_back(trk->get_pos_dump().Y());
			rec_track_z_dump.push_back(trk->get_pos_dump().Z());
			rec_track_px_dump.push_back(trk->get_mom_dump().Px());
			rec_track_py_dump.push_back(trk->get_mom_dump().Py());
			rec_track_pz_dump.push_back(trk->get_mom_dump().Pz());

			// Now fill hit IDs for this track
			std::vector<int> hit_ids;
			std::vector<double> hit_ids_pos_x;
			std::vector<double> hit_ids_pos_y;
			for (int i_hit = 0; i_hit < trk->get_num_hits(); ++i_hit) {
				//cout <<"charge: "<< trk->get_charge()<< "  hit id of the tracks: " << trk->get_hit_id(i_hit) << endl;
				hit_ids.push_back(trk->get_hit_id(i_hit));
				auto [det_id, det_z_pos] = GetDetElemIDFromHitID(trk->get_hit_id(i_hit));
				Double_t x = 0.0, y = 0.0;
				trk->getExpPositionFast(static_cast<Double_t>(det_z_pos), x, y);
				//cout << "x, y pos: (" << x << ", " << y << ")" << endl;
				hit_ids_pos_x.push_back(x);
				hit_ids_pos_y.push_back(y);
			}
			rec_hit_ids.push_back(hit_ids);
			rec_track_hit_x.push_back(hit_ids_pos_x);
			rec_track_hit_y.push_back(hit_ids_pos_y);
		}

		if(reco_dimu_mode==true){
			ResetRecoDimuBranches();
			for (auto it = m_sq_dim_vec->begin(); it != m_sq_dim_vec->end(); it++) {
				SRecDimuon& sdim = dynamic_cast<SRecDimuon&>(**it);
				int trk_id_pos = sdim.get_track_id_pos();
				int trk_id_neg = sdim.get_track_id_neg();
				SRecTrack& trk_pos = dynamic_cast<SRecTrack&>(*(m_sq_trk_vec->at(trk_id_pos)));
				SRecTrack& trk_neg = dynamic_cast<SRecTrack&>(*(m_sq_trk_vec->at(trk_id_neg)));


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
				// ===== Dump hypothesis =====
				sdim.calcVariables(2); // 2 = dump
				rec_dimuon_px_pos_dump.push_back(sdim.p_pos_dump.Px());
				rec_dimuon_py_pos_dump.push_back(sdim.p_pos_dump.Py());
				rec_dimuon_pz_pos_dump.push_back(sdim.p_pos_dump.Pz());
				rec_dimuon_px_neg_dump.push_back(sdim.p_neg_dump.Px());
				rec_dimuon_py_neg_dump.push_back(sdim.p_neg_dump.Py());
				rec_dimuon_pz_neg_dump.push_back(sdim.p_neg_dump.Pz());

				TLorentzVector mom_tgt = sdim.p_pos_target + sdim.p_neg_target;
				rec_dimuon_mass_tgt.push_back(mom_tgt.M());
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
	gCharge.clear(); gTrackID.clear();
	gvx.clear(); gvy.clear(); gvz.clear();
	gpx.clear(); gpy.clear(); gpz.clear();
	gx_st1.clear(); gy_st1.clear(); gz_st1.clear();
	gpx_st1.clear(); gpy_st1.clear(); gpz_st1.clear();
	gx_st3.clear(); gy_st3.clear(); gz_st3.clear();
	gpx_st3.clear(); gpy_st3.clear(); gpz_st3.clear();
}

void DimuAnaRUS::ResetRecoBranches() {
    rec_track_id.clear();
    rec_track_charge.clear();
    rec_track_vx.clear(); rec_track_vy.clear(); rec_track_vz.clear();
    rec_track_px.clear(); rec_track_py.clear(); rec_track_pz.clear();
    rec_track_x_st1.clear(); rec_track_y_st1.clear(); rec_track_z_st1.clear();
    rec_track_px_st1.clear(); rec_track_py_st1.clear(); rec_track_pz_st1.clear();
    rec_track_x_st3.clear(); rec_track_y_st3.clear(); rec_track_z_st3.clear();
    rec_track_px_st3.clear(); rec_track_py_st3.clear(); rec_track_pz_st3.clear();
    rec_track_chisq.clear();
    rec_track_chisq_tgt.clear();
    rec_track_chisq_dump.clear();
    rec_track_chisq_upstream.clear();
    rec_track_num_hits.clear();
    rec_track_x_tgt.clear(); rec_track_y_tgt.clear(); rec_track_z_tgt.clear();
    rec_track_px_tgt.clear(); rec_track_py_tgt.clear(); rec_track_pz_tgt.clear();
    rec_track_x_dump.clear(); rec_track_y_dump.clear(); rec_track_z_dump.clear();
    rec_track_px_dump.clear(); rec_track_py_dump.clear(); rec_track_pz_dump.clear();
    rec_hit_ids.clear();
    rec_track_hit_x.clear();
    rec_track_hit_y.clear();
}

void DimuAnaRUS::ResetRecoDimuBranches() {
    rec_dimuon_id.clear(); rec_dimuon_true_id.clear(); rec_dimuon_track_id_pos.clear(); rec_dimuon_track_id_neg.clear();
    rec_dimuon_x.clear(); rec_dimuon_y.clear(); rec_dimuon_z.clear();
    rec_dimuon_px_pos.clear(); rec_dimuon_py_pos.clear(); rec_dimuon_pz_pos.clear();
    rec_dimuon_px_neg.clear(); rec_dimuon_py_neg.clear(); rec_dimuon_pz_neg.clear();
    rec_dimuon_px_pos_tgt.clear(); rec_dimuon_py_pos_tgt.clear(); rec_dimuon_pz_pos_tgt.clear();
    rec_dimuon_px_neg_tgt.clear(); rec_dimuon_py_neg_tgt.clear(); rec_dimuon_pz_neg_tgt.clear();
    rec_dimuon_px_pos_dump.clear(); rec_dimuon_py_pos_dump.clear(); rec_dimuon_pz_pos_dump.clear();
    rec_dimuon_px_neg_dump.clear(); rec_dimuon_py_neg_dump.clear(); rec_dimuon_pz_neg_dump.clear();
    rec_dimuon_mass_tgt.clear();

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

std::pair<int, int> DimuAnaRUS::GetDetElemIDFromHitID(int hit_id) const {
    if (!m_hit_vec) return {-1, -1}; 

    for (const auto* hit : *m_hit_vec) {
        if (hit->get_hit_id() == hit_id) {
            GeomSvc::UseDbSvc(true);
            GeomSvc* geom = GeomSvc::instance();
            //cout << "det: "<< hit->get_detector_id() <<" detector z pos: "<< geom->getDetectorZ0(geom->getDetectorName(hit->get_detector_id())) <<endl;
            return {hit->get_detector_id(), geom->getDetectorZ0(geom->getDetectorName(hit->get_detector_id()))};
        }
    }

    return {-1, -1}; // Not found
}

