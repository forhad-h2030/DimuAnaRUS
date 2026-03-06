#ifndef _DIMU_ANA_RUS_H_
#define _DIMU_ANA_RUS_H_
#include <fun4all/SubsysReco.h>
#include <UtilAna/TrigRoadset.h>
#include <unordered_map>


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
class SQHit;

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
	bool saveDimuonOnly;
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
		mc_trig_mode = !enable; // Automatically disable mc_trig_mode when data_trig_mode is enabled
	}   

	void SetMCTriggerEmu(bool enable) { 
		mc_trig_mode = enable; 
		data_trig_mode = !enable; // Automatically disable data_trig_mode when mc_trig_mode is enabled
	}   

	void SetSaveOnlyDimuon(bool enable) { saveDimuonOnly = enable; }
	void SetSourceFlag(int flag) { SourceFlag = flag; }
	void SetProcessId(int proc_id_) { proc_id = proc_id_; }


    std::pair<int, int> GetDetElemIDFromHitID(int hit_id) const;
    std::pair<int, int> GetDetElemIDFromHitID(int hit_id, const std::unordered_map<int, SQHit*>& hit_map) const;


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
	void ResetHitBranches();
	void ResetTrueBranches();
	void ResetRecoBranches();
	void ResetRecoDimuBranches();

	private:

	// Scalar variables
	int eventID;
	int spillID;
	int runID;
	int rfID;
	int turnID;
	int rfIntensity[33];
	int fpgaTrigger[5] = {0};
	int nimTrigger[5] = {0};

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
    rec_dimuon_px_pos_dump, rec_dimuon_py_pos_dump, rec_dimuon_pz_pos_dump,
    rec_dimuon_px_neg_dump, rec_dimuon_py_neg_dump, rec_dimuon_pz_neg_dump,
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


	std::vector<int>
    rec_track_id,
    rec_track_charge,
    rec_track_num_hits;
	std::vector<std::vector<int>> rec_hit_ids;
	std::vector<std::vector<double>> rec_track_hit_x;
	std::vector<std::vector<double>> rec_track_hit_y;

	std::vector<double>
    rec_track_vx, rec_track_vy, rec_track_vz,
    rec_track_px, rec_track_py, rec_track_pz,
    rec_track_x_st1, rec_track_y_st1, rec_track_z_st1,
    rec_track_px_st1, rec_track_py_st1, rec_track_pz_st1,
    rec_track_x_st3, rec_track_y_st3, rec_track_z_st3,
    rec_track_px_st3, rec_track_py_st3, rec_track_pz_st3,
    rec_track_chisq, rec_track_chisq_tgt, rec_track_chisq_dump, rec_track_chisq_upstream,
    rec_track_x_tgt, rec_track_y_tgt, rec_track_z_tgt,
    rec_track_px_tgt, rec_track_py_tgt, rec_track_pz_tgt,
    rec_track_x_dump, rec_track_y_dump, rec_track_z_dump,
    rec_track_px_dump, rec_track_py_dump, rec_track_pz_dump;
};

#endif // _DimuAnaRUS.h_
