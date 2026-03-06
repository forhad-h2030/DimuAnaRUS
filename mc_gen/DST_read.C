R__LOAD_LIBRARY(interface_main)

int DST_read(const unsigned int n_evt_max=0)
{
  string fn_dst = "DST.root";

  TFile* file = new TFile(fn_dst.c_str());
  if (! file->IsOpen()) {
    cout << "Cannot open '" << fn_dst << "'.  Abort." << endl;
    exit(1);
  }

  ///
  /// Get the RUN node.
  ///
  TTree* tree_run = (TTree*)file->Get("T1");
  if (! tree_run) {
    cout << "Cannot get the RUN tree.  Abort." << endl;
    exit(1);
  }
  if (tree_run->GetEntries() != 1) {
    cout << "Find " << tree_run->GetEntries() << " entries (!= 1) in the RUN tree.  Abort." << endl;
    exit(1);
  }
  SQRun* sq_run = 0;
  tree_run->SetBranchAddress("RUN.SQRun", &sq_run);
  tree_run->GetEntry(0);
  int run_id = sq_run->get_run_id();
  cout << "Run ID = " << run_id << endl;

  ///
  /// Get the DST node.
  ///
  TTree* tree = (TTree*)file->Get("T");
  if (! tree) {
    cout << "Cannot get the DST tree.  Abort." << endl;
    exit(1);
  }
  unsigned int n_evt = tree->GetEntries();
  cout << "N of events = " << n_evt << endl;
  if (n_evt_max != 0 && n_evt > n_evt_max) n_evt = n_evt_max;

  SQEvent* sq_evt = 0;
  SQHitVector* sq_hit_vec = 0;
  tree->SetBranchAddress("DST.SQEvent"    , &sq_evt);
  tree->SetBranchAddress("DST.SQHitVector", &sq_hit_vec);

  for (int i_evt = 0; i_evt < n_evt; i_evt++) {
    tree->GetEntry(i_evt);
    if (i_evt == 0) cout << "Class " << sq_evt->ClassName() << " " << sq_hit_vec->ClassName() << endl;
    int evt_id    = sq_evt->get_event_id();
    int trig_bits = sq_evt->get_trigger();
    size_t n_hit = sq_hit_vec->size();
    cout << "E " << i_evt << " " << evt_id << " " << trig_bits << " " << n_hit << endl;
    for (size_t i_hit = 0; i_hit < n_hit; i_hit++) {
      SQHit* hit = sq_hit_vec->at(i_hit);
      int hit_id = hit->get_hit_id();
      int det_id = hit->get_detector_id();
      int ele_id = hit->get_element_id();
      double tdc = hit->get_tdc_time();
      cout << "  H " << hit_id << " " << det_id << " " << ele_id << " " << tdc << endl;
    }
    if (i_evt > 10) break;
  }

  exit(0);
}
