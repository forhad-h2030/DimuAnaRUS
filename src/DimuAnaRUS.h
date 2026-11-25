#ifndef _DIMU_ANA_RUS_H_
#define _DIMU_ANA_RUS_H_
#include <fun4all/SubsysReco.h>
#include <UtilAna/TrigRoadset.h>


class TFile;
class TTree;
class SQRun;
class SQSpillMap;
class SQEvent;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;
class SQMCEvent;
class TrackletVector;

class DimuAnaRUS: public SubsysReco {
	SQEvent* m_evt;
	SQMCEvent* m_evt_true;
	SQSpillMap* m_sp_map;
	SQHitVector* m_hit_vec;
	SQTrackVector * m_vec_trk;
	SQTrackVector*  m_sq_trk_vec;
    TrackletVector* trackletVec;
	SQDimuonVector* m_sq_dim_vec;
	SQDimuonVector* m_true_vec_dim;
	SQRun* sq_run;

	std::string m_file_name;
	std::string m_tree_name;
	TFile*      m_file;
	TTree*      m_tree;
	public:

	UtilTrigger::TrigRoadset m_rs;
	bool true_mode;
	bool reco_mode;
	bool data_trig_mode;
	bool mc_trig_mode;
	bool true_dimu_mode;
	int proc_id;
	int SourceFlag;
    bool sqhit_flag;
    bool reco_dimu_mode;

	void SetMCTrueMode(bool enable) { true_mode = enable; }
	void SetMCDimuonMode(bool enable) { true_dimu_mode = enable; }
	void SetRecoMode(bool enable) { reco_mode = enable; }
	void SetRecoDimuMode(bool enable) { reco_dimu_mode = enable; }
	void SetDataTriggerEmu(bool enable) { 
		data_trig_mode = enable; 
		if (enable) mc_trig_mode = false; // the two modes are mutually exclusive; disabling one leaves both off
	}   

	void SetMCTriggerEmu(bool enable) { 
		mc_trig_mode = enable; 
		if (enable) data_trig_mode = false; // the two modes are mutually exclusive; disabling one leaves both off
	}   

	void SetSourceFlag(int flag) { SourceFlag = flag; }
	void SetProcessId(int proc_id_) { proc_id = proc_id_; }




	unsigned int EncodeProcess(int processID, int sourceFlag);
	static int DecodeSourceFlag(unsigned int encoded);
	static int DecodeProcessID(unsigned int encoded);

	DimuAnaRUS(const std::string& name="DimuAnaRUS");
	virtual ~DimuAnaRUS();
	int Init(PHCompositeNode *startNode);
	int InitRun(PHCompositeNode *startNode);
	int process_event(PHCompositeNode *startNode);
	int End(PHCompositeNode *startNode);
	void SetOutputFileName(const std::string name) { m_file_name = name; }
    void EnableSQHit(bool enable) { sqhit_flag = enable; }

	void SetTreeName(const std::string& name) { m_tree_name = name; }
	void SetFileName(const std::string& name) { m_file_name = name; }
	static double GetPtReWeight(double pT, double xF, double mass);
	void ResetHitBranches();
	void ResetTrueBranches();
	void ResetRecoDimuBranches();

	private:

	// Scalar variables
	int eventID;
	double weight;
	int spillID;
	int runID;
	int evtPassTrigEmu = 0; // 1 if any reco dimuon passes the trigger emulation
	std::vector<int> rec_dimuon_pass_trig; // per-dimuon trigger emulation result

	// True dimuon (first one in SQTruthDimuonVector)
	double true_dimuon_mass, true_dimuon_pT, true_dimuon_xF;
	double ptReWeight;

	std::vector<int> hitID;
	std::vector<int> processID;
	std::vector<int> hitTrackID;
	std::vector<int> detectorID;
	std::vector<int> elementID;
	std::vector<double> tdcTime;
	std::vector<double> driftDistance;
	//std::vector<bool> hitsInTime;
	// True track data

	std::vector<int> gCharge, gTrackID;

	std::vector<double>
    gvx, gvy, gvz,
    gpx, gpy, gpz,
    gx_st1, gy_st1, gz_st1,
    gpx_st1, gpy_st1, gpz_st1,
    gx_st3, gy_st3, gz_st3,
    gpx_st3, gpy_st3, gpz_st3;

	std::vector<int>
    rec_dimuon_id, rec_dimuon_true_id, rec_dimuon_track_id_pos, rec_dimuon_track_id_neg;


	 std::vector<double>
    rec_dimuon_x, rec_dimuon_y, rec_dimuon_z,
    rec_dimuon_px_pos, rec_dimuon_py_pos, rec_dimuon_pz_pos,
    rec_dimuon_px_neg, rec_dimuon_py_neg, rec_dimuon_pz_neg,
    rec_dimuon_px_pos_tgt, rec_dimuon_py_pos_tgt, rec_dimuon_pz_pos_tgt,
    rec_dimuon_px_neg_tgt, rec_dimuon_py_neg_tgt, rec_dimuon_pz_neg_tgt,
    rec_dimuon_mass_tgt,
    rec_dimuon_x_pos_vtx, rec_dimuon_y_pos_vtx, rec_dimuon_z_pos_vtx,
    rec_dimuon_px_pos_vtx, rec_dimuon_py_pos_vtx, rec_dimuon_pz_pos_vtx,
    rec_dimuon_x_neg_vtx, rec_dimuon_y_neg_vtx, rec_dimuon_z_neg_vtx,
    rec_dimuon_px_neg_vtx, rec_dimuon_py_neg_vtx, rec_dimuon_pz_neg_vtx,
    rec_dimuon_x_pos_st1, rec_dimuon_y_pos_st1, rec_dimuon_z_pos_st1,
    rec_dimuon_x_neg_st1, rec_dimuon_y_neg_st1, rec_dimuon_z_neg_st1,
    rec_dimuon_x_pos_st3, rec_dimuon_y_pos_st3, rec_dimuon_z_pos_st3,
    rec_dimuon_x_neg_st3, rec_dimuon_y_neg_st3, rec_dimuon_z_neg_st3,
    rec_dimuon_px_pos_st1, rec_dimuon_py_pos_st1, rec_dimuon_pz_pos_st1,
    rec_dimuon_px_neg_st1, rec_dimuon_py_neg_st1, rec_dimuon_pz_neg_st1,
    rec_dimuon_px_pos_st3, rec_dimuon_py_pos_st3, rec_dimuon_pz_pos_st3,
    rec_dimuon_px_neg_st3, rec_dimuon_py_neg_st3, rec_dimuon_pz_neg_st3;


};

#endif // _DimuAnaRUS.h_
