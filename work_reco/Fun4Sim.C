R__LOAD_LIBRARY(libdimu_ana_rus)
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libPHPythia8)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libdptrigger)
R__LOAD_LIBRARY(libembedding)
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libcalibrator)

int Fun4Sim(const string Vect_in, const string DST_out, const int n_evt=20){
	const double FMAGSTR = -1.044;
	const double KMAGSTR = -1.025;

	recoConsts *rc = recoConsts::instance();
	rc->set_IntFlag("RUNNUMBER", 5433); /// The geometry is selected based on run number.
	rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
	rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
	rc->set_DoubleFlag("SIGX_BEAM", 0.3);
	rc->set_DoubleFlag("SIGY_BEAM", 0.3);
	rc->set_DoubleFlag("Z_UPSTREAM", -700.);

	Fun4AllServer* se = Fun4AllServer::instance();
	se->Verbosity(100);

	CalibHitElementPos* cal_ele_pos = new CalibHitElementPos();
	cal_ele_pos->CalibTriggerHit(false);
	se->registerSubsystem(cal_ele_pos);

	SQGeomAcc* geom_acc = new SQGeomAcc();
	geom_acc->SetMuonMode(SQGeomAcc::PAIR_TBBT); // PAIR, PAIR_TBBT, SINGLE, SINGLE_T, etc.
	geom_acc->SetPlaneMode(SQGeomAcc::HODO_CHAM); // HODO, CHAM or HODO_CHAM
	geom_acc->SetNumOfH1EdgeElementsExcluded(4); // Exclude 4 elements at H1 edges
	se->registerSubsystem(geom_acc);

	SQReco* reco = new SQReco();
	reco->Verbosity(1);
	reco->set_legacy_rec_container(false); 
	reco->set_geom_file_name((string)gSystem->Getenv("E1039_RESOURCE") + "/geometry/geom_run005433.root");
	reco->set_enable_KF(true);
	reco->setInputTy(SQReco::E1039);
	reco->setFitterTy(SQReco::KFREF);
	reco->set_evt_reducer_opt("none");
	reco->set_enable_eval_dst(true);
	for (int ii = 0; ii <= 3; ii++) reco->add_eval_list(ii);
	reco->set_enable_eval(true);
	se->registerSubsystem(reco);

	SQVertexing* vtx = new SQVertexing();
	vtx->Verbosity(1);
	se->registerSubsystem(vtx);


	Fun4AllRUSEventInputManager* in = new Fun4AllRUSEventInputManager("VectIn");
	in->Verbosity(99);
	in->set_tree_name("tree");
	in->SetMCMode(false);
	in->fileopen(Vect_in);
	//in->fileopen("/seaquest/users/mhossain/TrackMerging/out/combined_tracks_1.root");
	se->registerInputManager(in);
	///////////////////////////////////////////
	// Output
	// DST output manager
	//Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
	//se->registerOutputManager(out);
	//RUS output manager

	DimuAnaRUS* dimuAna = new DimuAnaRUS();
	dimuAna->SetTreeName("tree");
	dimuAna->SetOutputFileName(DST_out);
	//dimuAna->SetMCTrueMode(true);
	dimuAna->SetDataTriggerEmu(false); 
	dimuAna->SetRecoMode(true);
	se->registerSubsystem(dimuAna);
	//se->registerSubsystem(new DimuAnaRUS());

	se->run(n_evt);
	se->End();
	se->PrintTimer();
	rc->WriteToFile("recoConsts.tsv");
	std::cout << "All done" << std::endl;
	delete se;
	gSystem->Exit(0);
	return 0;
}
