#ifndef PHG4E1039TrackPairGen_H__
#define PHG4E1039TrackPairGen_H__

#include <g4main/PHG4ParticleGeneratorBase.h>
#include <map>
#include <vector>

#include <TGeoManager.h>
#include <phgeom/PHGeomUtility.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TLorentzVector.h>


class PHG4InEvent;
class PHCompositeNode;
class SQEvent;
class SQMCEvent;
class SQPrimaryVertexGen;

/// Event generator to generate a set of particles under a condition of the vertex and the momentum.
class PHG4E1039TrackPairGen : public PHG4ParticleGeneratorBase {

public:

  //! supported function distributions
  enum FUNCTION {Uniform,Gaus};

  PHG4E1039TrackPairGen(const std::string &name="EVTGENERATOR");
  virtual ~PHG4E1039TrackPairGen(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  //! interface for adding particles by name
  void add_particles(const std::string &name, const unsigned int count);

  //! interface for adding particle by pid
  void add_particles(const int pid, const unsigned int count);

  //! set the starting time for the event
  void set_t0(const double t0);
  
  //! range of randomized eta values
  void set_eta_range(const double eta_min, const double eta_max);

  //! range of randomized phi values
  void set_phi_range(const double phi_min, const double phi_max);

  //! range of randomized pt values
  //! \param[in] pt_gaus_width   if non-zero, further apply a Gauss smearing to the pt_min - pt_max flat distribution
  void set_pt_range(const double pt_min, const double pt_max, const double pt_gaus_width = 0);

  //! range of randomized p values
  //! \param[in] p_gaus_width   if non-zero, further apply a Gauss smearing to the p_min - p_max flat distribution
  void set_p_range(const double p_min, const double p_max, const double p_gaus_width = 0);

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z);

  //! set the mean value of the vertex distribution
  void set_vertex_distribution_mean(const double x, const double y, const double z);

  //! set the width of the vertex distribution function about the mean
  void set_vertex_distribution_width(const double x, const double y, const double z);

  //! set an offset vector from the existing vertex
  void set_existing_vertex_offset_vector(const double x, const double y, const double z);
  
  //! set the distribution function of particles about the vertex
  void set_vertex_size_function(FUNCTION r);
  bool GenerateMomentumWithoutExpPDF(int muon_counter, int pdgcode, TLorentzVector &muon1, TLorentzVector &muon2, double &px, double &py, double &pz);
  bool GenerateMomentumWithExpPDF(int muon_counter, double &px, double &py, double &pz);


  //! set the dimensions of the distribution of particles about the vertex
  void set_vertex_size_parameters(const double mean, const double width);
  bool gen_with_exp_pdf;
  void SetExpPDFMode(bool enable) { gen_with_exp_pdf = enable; }

  //!
  void set_pxpypz_range(
  		const double x_min, const double x_max,
			const double y_min, const double y_max,
			const double z_min, const double z_max
			);
void set_par1_pxpypz_range(
                const double x_min, const double x_max,
                const double y_min, const double y_max,
                const double z_min, const double z_max
                );

void set_par2_pxpypz_range(
                const double x_min, const double x_max,
                const double y_min, const double y_max,
                const double z_min, const double z_max
                );

  void set_max_opening_angle(const double x_max);
  void set_max_muon_angle_with_z(const double x_max);
   //! Enable legacy vertex gen
   void enableLegacyVtxGen() { _legacy_vertexgenerator = true; } //Abi
  double _maxProbability;
  void set_bend_range(const double lo_pos, const double hi_pos,
                      const double lo_neg, const double hi_neg)
  {
    _bend_lo_pos = lo_pos;
    _bend_hi_pos = hi_pos;
    _bend_lo_neg = lo_neg;
    _bend_hi_neg = hi_neg;
    _apply_bend_cut = true;
  }

private:

  double smearvtx(const double position, const double width, FUNCTION dist) const;
  // these need to be stored separately until run time when the names
  // can be translated using the GEANT4 lookup
  std::vector<std::pair<int, unsigned int> > _particle_codes; // <pdgcode, count>
  std::vector<std::pair<std::string, unsigned int> > _particle_names; // <names, count>  
  FUNCTION _vertex_func_x;
  FUNCTION _vertex_func_y;
  FUNCTION _vertex_func_z;
  double _t0;
  double _vertex_x;
  double _vertex_y;
  double _vertex_z;
  double _vertex_width_x;
  double _vertex_width_y;
  double _vertex_width_z;
  double _vertex_offset_x;
  double _vertex_offset_y;
  double _vertex_offset_z;
  FUNCTION _vertex_size_func_r;
  double _vertex_size_mean;
  double _vertex_size_width;
  double _eta_min;
  double _eta_max;
  double _phi_min;
  double _phi_max;
  double _pt_min;
  double _pt_max;
  double _pt_gaus_width;
  double _p_min;
  double _p_max; 
  double _p_gaus_width;

  double _px_min, _px_max;
  double _py_min, _py_max;
  double _pz_min, _pz_max;

  // For mup
  double _px_par1_min, _px_par1_max;
  double _py_par1_min, _py_par1_max;
  double _pz_par1_min, _pz_par1_max;

  // For mum
  double _px_par2_min, _px_par2_max;
  double _py_par2_min, _py_par2_max;
  double _pz_par2_min, _pz_par2_max;

  double _theta_min, _theta_max, _theta_mu_max;
  int    _eventcount;
  PHG4InEvent* _ineve;
  SQEvent* _evt; //< An output node
  SQMCEvent* _mcevt; //< An output node
  TH2F *hProb;
  TFile* _normFlowFile;         // ROOT file for norm flow momentum map
  THnSparseD* _hMomentumMap;
  static constexpr double MUON_MASS = 0.1056;    // GeV/c^2
  static constexpr double BEAM_ENERGY = 120.0;   // GeV
  static constexpr double PROTON_MASS = 0.938;   // GeV/c^2


   // Bend-plane (px/pz) cuts
  double _bend_lo_pos;   // mu+ lower bound
  double _bend_hi_pos;   // mu+ upper bound
  double _bend_lo_neg;   // mu- lower bound
  double _bend_hi_neg;   // mu- upper bound
  bool   _apply_bend_cut;

  // Dimuon pz < BEAM_ENERGY cut
  bool   _apply_dimuon_pz_cut;


  bool _legacy_vertexgenerator;
  SQPrimaryVertexGen* _vertexGen;
};

#endif
