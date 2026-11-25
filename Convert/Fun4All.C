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

int Fun4All(const string DST_in, const string DST_out, const int n_evt=0){
	Fun4AllServer* se = Fun4AllServer::instance();
	Fun4AllInputManager *in = new Fun4AllDstInputManager("DUMMY");
	se->registerInputManager(in);
	in->fileopen(DST_in);
	
	DimuAnaRUS* dimuAna = new DimuAnaRUS();
        dimuAna->SetTreeName("tree");
	dimuAna->SetMCTrueMode(false);
        dimuAna->SetOutputFileName("RUS.root");
        dimuAna->SetRecoMode(false);
        se->registerSubsystem(dimuAna);

	se->run(n_evt);
	se->End();
	se->PrintTimer();
	std::cout << "All done" << std::endl;
	delete se;
	gSystem->Exit(0);
	return 0;
}
