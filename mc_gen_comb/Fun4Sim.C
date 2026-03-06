#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
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
R__LOAD_LIBRARY(fun4all_rus_file_manager)

using namespace std;

int Fun4Sim(const int nevent = 10)
{
	const double target_coil_pos_z = -300;
	const int nmu = 1;
	int embedding_opt = 0;
	const bool legacy_rec_container = true;

	const bool do_collimator = true;
	const bool do_target = true;
	const bool do_e1039_shielding = true;
	const bool do_fmag = true;
	const bool do_kmag = true;
	const bool do_absorber = true;
	const bool do_dphodo = true;
	const bool do_station1DC = false;       //station-1 drift chamber should be turned off by default

	const double target_l = 7.9; //cm
	const double target_z = (7.9-target_l)/2.; //cm
	const int use_g4steps = 1;

	const double FMAGSTR = -1.044;
	const double KMAGSTR = -1.025;

	//! Particle generator flag.  Only one of these must be true.
	const bool gen_particle_with_exp_pdf =false;
	const bool read_hepmc   = false;

	//! Use SQPrimaryVertexGen or not.
	const bool SQ_vtx_gen = true;

	recoConsts *rc = recoConsts::instance();
	rc->set_IntFlag("RUNNUMBER", 5433); /// The geometry is selected based on run number.
	rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
	rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
	rc->set_DoubleFlag("SIGX_BEAM", 0.3);
	rc->set_DoubleFlag("SIGY_BEAM", 0.3);
	rc->set_DoubleFlag("Z_UPSTREAM", -700.);

	if(SQ_vtx_gen) { // cf. SQPrimaryVertexGen
		rc->set_CharFlag("VTX_GEN_MATERIAL_MODE", "Target"); // All, Target, Dump, TargetDumpGap or Manual
		//rc->set_CharFlag("VTX_GEN_MATERIAL_MODE", "Target"); // All, Target, Dump, TargetDumpGap or Manual
		//rc->set_DoubleFlag("VTX_GEN_Z_START",  50.0); // For "Manual"
		//rc->set_DoubleFlag("VTX_GEN_Z_STOP" , 100.0); // For "Manual"
	}
	rc->Print();

	GeomSvc *geom_svc = GeomSvc::instance();
	//const double x0_shift = 0.0; //cm 
	//std::cout << "D2X::X0: " << geom_svc->getDetectorX0("D2X") << std::endl;
	//geom_svc->setDetectorX0("D2X", geom_svc->getDetectorX0("D2X")+x0_shift);
	//std::cout << "D2X::X0: " << geom_svc->getDetectorX0("D2X") << std::endl;

	///////////////////////////////////////////
	// Make the Server
	//////////////////////////////////////////
	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(0);


	if(read_hepmc) {
		HepMCNodeReader *hr = new HepMCNodeReader();
		hr->set_particle_filter_on(true);
		hr->insert_particle_filter_pid(13);
		hr->insert_particle_filter_pid(-13);
		se->registerSubsystem(hr);
	}

		PHG4E1039TrackPairGen *comb = new PHG4E1039TrackPairGen("MUP");
		comb->set_seed(125);
		comb->add_particles("mu+", nmu);  // mu+,e+,proton,pi+,Upsilon
		comb->set_seed(120);
		comb->add_particles("mu-", nmu);
		if (SQ_vtx_gen) comb->enableLegacyVtxGen();
		else{
			comb->set_vertex_distribution_function(PHG4E1039TrackPairGen::Uniform,
					PHG4E1039TrackPairGen::Uniform,
					PHG4E1039TrackPairGen::Uniform);
			comb->set_vertex_distribution_mean(0.0, 0.0, target_coil_pos_z);
			comb->set_vertex_distribution_width(0.0, 0.0, 0.0);
			comb->set_vertex_size_function(PHG4E1039TrackPairGen::Uniform);
			comb->set_vertex_size_parameters(0.0, 0.0);
		}
	if(gen_particle_with_exp_pdf) comb->SetExpPDFMode(true);
    else{
		comb->set_par1_pxpypz_range(-6.0,6.0, -4,4, 15,80);
		comb->set_par2_pxpypz_range(-6.0,6.0, -4,4, 15,80);
		comb->set_max_opening_angle(5);
        comb->set_bend_range(-0.15, 0.0,   // mu+: lo, hi
                              0.0,  0.15); // mu-: lo, hi
	}
        //comb->Verbosity(1);
		se->registerSubsystem(comb);

	// Fun4All G4 module
	PHG4Reco *g4Reco = new PHG4Reco();
	//PHG4Reco::G4Seed(123);
	g4Reco->set_field_map();
	// size of the world - every detector has to fit in here
	g4Reco->SetWorldSizeX(1000);
	g4Reco->SetWorldSizeY(1000);
	g4Reco->SetWorldSizeZ(5000);
	// shape of our world - it is a tube
	g4Reco->SetWorldShape("G4BOX");
	// this is what our world is filled with
	g4Reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
	// Geant4 Physics list to use
	g4Reco->SetPhysicsList("FTFP_BERT");

	// insensitive elements of the spectrometer
	SetupInsensitiveVolumes(g4Reco, do_e1039_shielding, do_fmag, do_kmag, do_absorber);

	// collimator, targer and shielding between target and FMag
	SetupBeamline(g4Reco, do_collimator, target_coil_pos_z - 302.36); // Is the position correct??

	if (do_target) {
		SetupTarget(g4Reco, target_coil_pos_z, target_l, target_z, use_g4steps);
	}

	// sensitive elements of the spectrometer
	SetupSensitiveDetectors(g4Reco, do_dphodo, do_station1DC);

	se->registerSubsystem(g4Reco);

	// save truth info to the Node Tree
	PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
	g4Reco->registerSubsystem(truth);

	// digitizer
	SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
	//digitizer->Verbosity(99);
	digitizer->set_enable_st1dc(do_station1DC);    // these two lines need to be in sync with the parameters used
	digitizer->set_enable_dphodo(do_dphodo);       // in the SetupSensitiveVolumes() function call above
	se->registerSubsystem(digitizer);

	se->registerSubsystem(new TruthNodeMaker());

	//Apply additonal cut after event generation
        //MuonTrackFilter* muon_filter = new MuonTrackFilter();
        //muon_filter->SetAngleThreshold(0.0, 50.0); //in degree
        //se->registerSubsystem(muon_filter);
	/// Save only events that are in the geometric acceptance.
	SQGeomAcc* geom_acc = new SQGeomAcc();
	geom_acc->SetMuonMode(SQGeomAcc::PAIR_TBBT); // PAIR, PAIR_TBBT, SINGLE, SINGLE_T, etc.
	geom_acc->SetPlaneMode(SQGeomAcc::HODO_CHAM); // HODO, CHAM or HODO_CHAM
	geom_acc->SetNumOfH1EdgeElementsExcluded(4); // Exclude 4 elements at H1 edges
	se->registerSubsystem(geom_acc);
	se->registerSubsystem(new TruthNodeMaker());

	if(embedding_opt == 1) {
		SRawEventEmbed *embed = new SRawEventEmbed("SRawEventEmbed");
		embed->set_in_name("digit_016070_R007.root");
		embed->set_in_tree_name("save");
		embed->set_trigger_bit((1<<0));
		//embed->set_in_name("random_run3a_1.root");
		//embed->set_in_tree_name("mb");
		//embed->set_trigger_bit((1<<7));
		embed->Verbosity(0);
		se->registerSubsystem(embed);
	}

	// Trigger Emulator
	DPTriggerAnalyzer* dptrigger = new DPTriggerAnalyzer();
	dptrigger->set_road_set_file_name("$E1039_RESOURCE/trigger/trigger_67.txt");
	se->registerSubsystem(dptrigger);

	// Event Filter
	//EvtFilter *evt_filter = new EvtFilter();
	//evt_filter->Verbosity(10);
	//evt_filter->set_trigger_req(1<<5);
	//se->registerSubsystem(evt_filter);
	// Tracking module
	// input - we need a dummy to drive the event loop

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
	if(read_hepmc) {
		Fun4AllHepMCInputManager *in = new Fun4AllHepMCInputManager("HEPMCIN");
		in->Verbosity(10);
		in->set_vertex_distribution_mean(0,0,target_coil_pos_z,0);
		se->registerInputManager(in);
		in->fileopen("hepmcout.txt");
	} else {
		Fun4AllInputManager *in = new Fun4AllDummyInputManager("DUMMY");
		se->registerInputManager(in);
	}

	///////////////////////////////////////////
	// Output
	///////////////////////////////////////////
	DimuAnaRUS* dimuAna = new DimuAnaRUS();
    dimuAna->SetTreeName("tree");
    dimuAna->SetMCTrueMode(true);
    dimuAna->SetOutputFileName("RUS.root");
    dimuAna->SetSaveOnlyDimuon(true);
    dimuAna->SetRecoMode(true);
    dimuAna->SetRecoDimuMode(true);
    dimuAna->EnableSQHit(true);
	dimuAna->SetProcessId(15);
    se->registerSubsystem(dimuAna);

    //Fun4AllRUSOutputManager* tree = new Fun4AllRUSOutputManager();
    //tree->SetTreeName("tree");
    //tree->SetFileName("RUS.root");
    //se->registerOutputManager(tree);

	const bool count_only_good_events = true;
	se->run(nevent, count_only_good_events);
	PHGeomUtility::ExportGeomtry(se->topNode(),"geom.root");

	// finish job - close and save output files
	se->End();
	se->PrintTimer();
	rc->WriteToFile("recoConsts.tsv");
	std::cout << "All done" << std::endl;
	// cleanup - delete the server and exit
	delete se;
	gSystem->Exit(0);
	return 0;
}
