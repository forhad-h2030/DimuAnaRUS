#include "PHG4E1039TrackPairGen.h"
#include <cmath>
#include "TFile.h"
#include "TH2F.h"
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4InEvent.h>
#include "g4main/PHG4VtxPoint.h"
#include "g4main/PHG4TruthInfoContainer.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQMCEvent_v1.h>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <E906LegacyVtxGen/SQPrimaryVertexGen.h>

using namespace std;

PHG4E1039TrackPairGen::PHG4E1039TrackPairGen(const string &name): 
  PHG4ParticleGeneratorBase(name),
  _particle_codes(),
  _particle_names(),
  _vertex_func_x(Uniform),
  _vertex_func_y(Uniform),
  _vertex_func_z(Uniform),
  _t0(0.0),
  _vertex_x(0.0),
  _vertex_y(0.0),
  _vertex_z(0.0),
  _vertex_width_x(0.0),
  _vertex_width_y(0.0),
  _vertex_width_z(0.0),
  _vertex_offset_x(0.0),
  _vertex_offset_y(0.0),
  _vertex_offset_z(0.0),
  _vertex_size_func_r(Uniform),
  _vertex_size_mean(0.0),
  _vertex_size_width(0.0),
  _eta_min(-1.25),
  _eta_max(1.25),
  _phi_min(-M_PI),
  _phi_max(M_PI),
  _pt_min(0.0),
  _pt_max(5.0),
  _pt_gaus_width(0.0),
  _p_min(0.0),
  _p_max(NAN),
  _p_gaus_width(NAN),
  _px_min(NAN), _px_max(NAN),
  _py_min(NAN), _py_max(NAN),
  _pz_min(NAN), _pz_max(NAN),
  _theta_max(NAN),
  _px_par1_min(-6.0),
  _py_par1_min(-4.0),
  _pz_par1_min(10.0),
  _px_par1_max(6.0),
  _py_par1_max(4.0),
  _pz_par1_max(100.0),
  _px_par2_min(-6.0),
  _py_par2_min(-4.0),
  _pz_par2_min(10.0),
  _px_par2_max(6.0),
  _py_par2_max(4.0),
  _pz_par2_max(100.0),
  _eventcount(0),
  _ineve(NULL), 
  _evt(NULL),
  _mcevt(NULL),
  gen_with_exp_pdf(false),
  _legacy_vertexgenerator(false),
  _bend_lo_pos(-1.0),
  _bend_hi_pos( 1.0),
  _bend_lo_neg(-1.0),
  _bend_hi_neg( 1.0),
  _apply_bend_cut(false),
  _apply_dimuon_pz_cut(true)
{

	_vertexGen = new SQPrimaryVertexGen();
	return;
}

void PHG4E1039TrackPairGen::add_particles(const std::string &name, const unsigned int num) {
  _particle_names.push_back(std::make_pair(name,num));
  return;
}

void PHG4E1039TrackPairGen::add_particles(const int pid, const unsigned int num) {
  _particle_codes.push_back(std::make_pair(pid,num));
  return;
}

void PHG4E1039TrackPairGen::set_t0(const double t0) {
  _t0 = t0;
  return;
}

void PHG4E1039TrackPairGen::set_eta_range(const double min, const double max) {
  if (min > max)
    {
      cout << "not setting eta bc etamin " << min << " > etamax: " << max << endl;
      return;
    }
  _eta_min = min;
  _eta_max = max;
  return;
}

void PHG4E1039TrackPairGen::set_phi_range(const double min, const double max) {
  if (min > max)
    {
      cout << "not setting phi bc phimin " << min << " > phimax: " << max << endl;
      return;
    }
  _phi_min = min;
  _phi_max = max;
  return;
}

void PHG4E1039TrackPairGen::set_pt_range(const double min, const double max, const double pt_gaus_width) {
  if (min > max)
    {
      cout << "not setting pt bc ptmin " << min << " > ptmax: " << max << endl;
      return;
    }
  assert(pt_gaus_width >=0 );

  _pt_min = min;
  _pt_max = max;
  _pt_gaus_width = pt_gaus_width;
  _p_min = NAN;
  _p_max = NAN;
  _p_gaus_width = NAN;
  return;
}

void PHG4E1039TrackPairGen::set_p_range(const double min, const double max, const double p_gaus_width) {
  if (min > max)
    {
      cout << "not setting p bc ptmin " << min << " > ptmax: " << max << endl;
      return;
    }
  assert(p_gaus_width >=0 );
  _pt_min = NAN;
  _pt_max = NAN;
  _pt_gaus_width = NAN;
  _p_min = min;
  _p_max = max;
  _p_gaus_width = p_gaus_width;
  return;
}

void PHG4E1039TrackPairGen::set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z) {
  _vertex_func_x = x;
  _vertex_func_y = y;
  _vertex_func_z = z;
  return;
}

void PHG4E1039TrackPairGen::set_vertex_distribution_mean(const double x, const double y, const double z) {
  _vertex_x = x;
  _vertex_y = y;
  _vertex_z = z;
  return;
}

void PHG4E1039TrackPairGen::set_vertex_distribution_width(const double x, const double y, const double z) {
  _vertex_width_x = x;
  _vertex_width_y = y;
  _vertex_width_z = z;
  return;
}

void PHG4E1039TrackPairGen::set_existing_vertex_offset_vector(const double x, const double y, const double z) {
  _vertex_offset_x = x;
  _vertex_offset_y = y;
  _vertex_offset_z = z;
  return;
}

void PHG4E1039TrackPairGen::set_vertex_size_function(FUNCTION r) {
  _vertex_size_func_r = r;
  return;
}

void PHG4E1039TrackPairGen::set_vertex_size_parameters(const double mean, const double width) {
  _vertex_size_mean = mean;
  _vertex_size_width = width;
  return;
}

int PHG4E1039TrackPairGen::InitRun(PHCompositeNode *topNode) {

if(gen_with_exp_pdf){
	_normFlowFile = TFile::Open("NormMomentumMap.root");
	if (!_normFlowFile || _normFlowFile->IsZombie()) {
		cout << PHWHERE << "::Error - cannot open NormMomentumMap.root" << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}
	_hMomentumMap = dynamic_cast<THnSparseD*>(_normFlowFile->Get("MomentumMap"));
	if (!_hMomentumMap) {
		cout << PHWHERE << "::Error - MomentumMap histogram not found in file" << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}
	// Verify histogram dimensions
	if (_hMomentumMap->GetNdimensions() != 6) {
		cout << PHWHERE << "::Error - MomentumMap histogram does not have 6 dimensions" << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}
}

  if ((_vertex_func_x != Uniform)&&(_vertex_func_x != Gaus)) {
    cout << PHWHERE << "::Error - unknown vertex distribution function requested" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if ((_vertex_func_y != Uniform)&&(_vertex_func_y != Gaus)) {
    cout << PHWHERE << "::Error - unknown vertex distribution function requested" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if ((_vertex_func_z != Uniform)&&(_vertex_func_z != Gaus)) {
    cout << PHWHERE << "::Error - unknown vertex distribution function requested" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing.  ABORTRUN." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  _ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!_ineve) {
    _ineve = new PHG4InEvent();
    dstNode->addNode(new PHIODataNode<PHObject>(_ineve, "PHG4INEVENT", "PHObject"));
  }
  _evt = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_evt) {
    _evt = new SQEvent_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(_evt, "SQEvent", "PHObject"));
  }
  _mcevt = findNode::getClass<SQMCEvent>(topNode, "SQMCEvent");
  if (!_mcevt) {
    _mcevt = new SQMCEvent_v1();
    dstNode->addNode(new PHIODataNode<PHObject>(_mcevt, "SQMCEvent", "PHObject"));
  }

  if (verbosity > 0) {
    cout << "================ PHG4E1039TrackPairGen::InitRun() ======================" << endl;
    cout << " Random seed = " << get_seed() << endl;
    cout << " Particles:" << endl;
    for (unsigned int i=0; i<_particle_codes.size(); ++i) {
      cout << "    " << _particle_codes[i].first << ", count = " << _particle_codes[i].second << endl;
    }
    for (unsigned int i=0; i<_particle_names.size(); ++i) {
      cout << "    " << _particle_names[i].first << ", count = " << _particle_names[i].second << endl;
    }
    if (get_reuse_existing_vertex()) {
      cout << " Vertex Distribution: Set to reuse a previously generated sim vertex" << endl;
      cout << " Vertex offset vector (x,y,z) = (" << _vertex_offset_x << ","<< _vertex_offset_y << ","<< _vertex_offset_z << ")" << endl;
    } else {
      cout << " Vertex Distribution Function (x,y,z) = ("; 
      if (_vertex_func_x == Uniform) cout << "Uniform,";
      else if (_vertex_func_x == Gaus) cout << "Gaus,";
      if (_vertex_func_y == Uniform) cout << "Uniform,";
      else if (_vertex_func_y == Gaus) cout << "Gaus,";
      if (_vertex_func_z == Uniform) cout << "Uniform";
      else if (_vertex_func_z == Gaus) cout << "Gaus";
      cout << ")" << endl;
      cout << " Vertex mean (x,y,z) = (" << _vertex_x << ","<< _vertex_y << ","<< _vertex_z << ")" << endl;
      cout << " Vertex width (x,y,z) = (" << _vertex_width_x << ","<< _vertex_width_y << ","<< _vertex_width_z << ")" << endl;
    }
    cout << " Vertex size function (r) = (";
    if (_vertex_size_func_r == Uniform) cout << "Uniform";
    else if (_vertex_size_func_r == Gaus) cout << "Gaus";
    cout << ")" << endl;
    cout << " Vertex size (mean) = (" << _vertex_size_mean << ")" << endl;
    cout << " Vertex size (width) = (" << _vertex_size_width << ")" << endl;
    cout << " Eta range = " << _eta_min << " - " << _eta_max << endl;
    cout << " Phi range = " << _phi_min << " - " << _phi_max << endl;
    cout << " pT range = " << _pt_min << " - " << _pt_max << endl;
    cout << " t0 = " << _t0 << endl;

	if (_apply_bend_cut) {
      cout << " Bend cut (mu+): [" << _bend_lo_pos << ", " << _bend_hi_pos << "]" << endl;
      cout << " Bend cut (mu-): [" << _bend_lo_neg << ", " << _bend_hi_neg << "]" << endl;
    }
    cout << " Dimuon pz < BEAM_ENERGY cut: " << (_apply_dimuon_pz_cut ? "enabled" : "disabled") << endl;
    cout << "===========================================================================" << endl;
  }

  // the definition table should be filled now, so convert codes into names
  for (unsigned int i=0;i<_particle_codes.size();++i) {
    int pdgcode = _particle_codes[i].first;
    unsigned int count = _particle_codes[i].second;
    string pdgname = get_pdgname(pdgcode);
    _particle_names.push_back(make_pair(pdgname,count));
  }
  _vertexGen->InitRun(topNode);//abi

  recoConsts* rc = recoConsts::instance();
  ostringstream oss;
  for (unsigned int ii = 0; ii < _particle_names.size(); ++ii) {
    string      name = _particle_names[ii].first;
    unsigned int num = _particle_names[ii].second;
    if (ii > 0) oss << ":";
    oss << name << "*" << num;
  }
  rc->set_CharFlag("PHG4SEG_"+Name()+"_PARTICLES", oss.str());

  if(! std::isnan(_px_min)) {
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_PX_MIN", _px_min);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_PX_MAX", _px_max);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_PY_MIN", _py_min);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_PY_MAX", _py_max);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_PZ_MIN", _pz_min);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_PZ_MAX", _pz_max);
  }
  if (! _legacy_vertexgenerator) {
    rc->set_IntFlag   ("PHG4SEG_"+Name()+"_VTX_FUNC_X"  , _vertex_func_x);
    rc->set_IntFlag   ("PHG4SEG_"+Name()+"_VTX_FUNC_Y"  , _vertex_func_y);
    rc->set_IntFlag   ("PHG4SEG_"+Name()+"_VTX_FUNC_Z"  , _vertex_func_z);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_X"       , _vertex_x);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_Y"       , _vertex_y);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_Z"       , _vertex_z);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_WIDTH_X" , _vertex_width_x);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_WIDTH_Y" , _vertex_width_y);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_WIDTH_Z" , _vertex_width_z);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_OFFSET_X", _vertex_offset_x);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_OFFSET_Y", _vertex_offset_y);
    rc->set_DoubleFlag("PHG4SEG_"+Name()+"_VTX_OFFSET_Z", _vertex_offset_z);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4E1039TrackPairGen::process_event(PHCompositeNode *topNode) {
    if (verbosity > 0) {
        std::cout << "====================== PHG4E1039TrackPairGen::process_event() =====================" << std::endl;
        std::cout << "PHG4E1039TrackPairGen::process_event - reuse_existing_vertex = " << reuse_existing_vertex << std::endl;
    }

    if (!_ineve || !_evt || !_mcevt) {
        std::cout << PHWHERE << "::Error - PHG4InEvent, SQEvent, or SQMCEvent is null" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    if (_particle_names.size() != 2) {
        std::cout << PHWHERE << "::Error - Exactly 2 particles required, got " << _particle_names.size() << ": ";
        for (const auto& p : _particle_names) std::cout << p.first << " ";
        std::cout << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    int vtxindex = -1;
    int trackid = -1;
    TLorentzVector muon1;
    TLorentzVector muon2;
    int muon_counter = 0;

    // Vertex generation
    if (_legacy_vertexgenerator) {
        TVector3 vtx = _vertexGen->generateVertex();
        vtx_x = vtx.X();
        vtx_y = vtx.Y();
        vtx_z = vtx.Z();
    } else if (!ReuseExistingVertex(topNode)) {
        vtx_x = smearvtx(_vertex_x, _vertex_width_x, _vertex_func_x);
        vtx_y = smearvtx(_vertex_y, _vertex_width_y, _vertex_func_y);
        vtx_z = smearvtx(_vertex_z, _vertex_width_z, _vertex_func_z);
    }
    vtx_x += _vertex_offset_x;
    vtx_y += _vertex_offset_y;
    vtx_z += _vertex_offset_z;

    if (verbosity > 0) {
        std::cout << "PHG4E1039TrackPairGen::process_event - vertex center " << reuse_existing_vertex
                  << " (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ") cm" << std::endl;
    }

    if (vtx_z < -690) {
        std::cout << PHWHERE << "::Error - Vertex z = " << vtx_z << " < -690 cm" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    for (unsigned int i = 0; i < _particle_names.size(); ++i) {
        muon_counter++;
        std::string pdgname = _particle_names[i].first;
        int pdgcode = get_pdgcode(pdgname);
        if (pdgcode == 0) {
            std::cout << PHWHERE << "::Error - Invalid PDG code for particle " << pdgname << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
        unsigned int nparticles = _particle_names[i].second;

        for (unsigned int j = 0; j < nparticles; ++j) {
            if ((_vertex_size_width > 0.0) || (_vertex_size_mean != 0.0)) {
                if (!RandomGenerator) {
                    std::cout << PHWHERE << "::Error - RandomGenerator is null" << std::endl;
                    return Fun4AllReturnCodes::ABORTRUN;
                }
                if (_vertex_size_func_r != Uniform && _vertex_size_func_r != Gaus) {
                    std::cout << PHWHERE << "::Error - Invalid vertex size function" << std::endl;
                    return Fun4AllReturnCodes::ABORTRUN;
                }
                double r = smearvtx(_vertex_size_mean, _vertex_size_width, _vertex_size_func_r);
                double x = 0.0, y = 0.0, z = 0.0;
                gsl_ran_dir_3d(RandomGenerator, &x, &y, &z);
                x *= r;
                y *= r;
                z *= r;
                vtxindex = _ineve->AddVtx(vtx_x + x, vtx_y + y, vtx_z + z, _t0);
            } else if ((i == 0) && (j == 0)) {
                vtxindex = _ineve->AddVtx(vtx_x, vtx_y, vtx_z, _t0);
            }
            ++trackid;

            double px, py, pz;
            // Momentum generation
            if (gen_with_exp_pdf) {
                if (!GenerateMomentumWithExpPDF(muon_counter, px, py, pz)) {
                    return Fun4AllReturnCodes::ABORTRUN;
                }
            } else {
                if (!GenerateMomentumWithoutExpPDF(muon_counter, pdgcode, muon1, muon2, px, py, pz)) {
                    return Fun4AllReturnCodes::ABORTEVENT;
                }
            }

            if (muon_counter == 1) {
                muon1.SetXYZM(px, py, pz, MUON_MASS);
            } else if (muon_counter == 2) {
                muon2.SetXYZM(px, py, pz, MUON_MASS);
            }

            double m = get_mass(pdgcode);
            double e = sqrt(px * px + py * py + pz * pz + m * m);
            PHG4Particle *particle = new PHG4Particlev2();
            particle->set_track_id(trackid);
            particle->set_vtx_id(vtxindex);
            particle->set_parent_id(0);
            particle->set_name(pdgname);
            particle->set_pid(pdgcode);
            particle->set_px(px);
            particle->set_py(py);
            particle->set_pz(pz);
            particle->set_e(e);

            _ineve->AddParticle(vtxindex, particle);
            if (embedflag != 0) _ineve->AddEmbeddedParticle(particle, embedflag);
        }
    }

    _evt->set_run_id(0);
    _evt->set_spill_id(0);
    _evt->set_event_id(++_eventcount);
    _mcevt->set_cross_section(0.0);
    _mcevt->set_weight(1.0);

    if (verbosity > 0) {
        _ineve->identify();
        std::cout << "======================================================================================" << std::endl;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

bool PHG4E1039TrackPairGen::GenerateMomentumWithoutExpPDF(int muon_counter, int pdgcode,
        TLorentzVector &muon1, TLorentzVector &muon2,
        double &px, double &py, double &pz)
{
    if (!RandomGenerator) {
        std::cout << PHWHERE << "::Error - RandomGenerator is null" << std::endl;
        return false;
    }

    if (std::isnan(_theta_max)) {
        std::cout << PHWHERE << "::Error - _theta_max is not set (NAN)" << std::endl;
        return false;
    }

    int iteration = 0;
    const int max_iterations = 10000;
    double angle, xF;

    if (muon_counter == 1) {
        // ----------------------------------------------------------------
        // Generate mu+ (pdgcode == -13)
        // Cuts applied per iteration:
        //   1. pT in [_pt_min, _pt_max]
        //   2. Bend-plane ratio bend = px/pz in [_bend_lo_pos, _bend_hi_pos]  (if enabled)
        // ----------------------------------------------------------------
        if (verbosity > 0)
            std::cout << "Mu+ Vertex: (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ")" << std::endl;

        do {
            if (verbosity > 0)
                std::cout << "par1 pdgcode: " << pdgcode
                          << ", _pz_par1_min: " << _pz_par1_min
                          << ", _pz_par1_max: " << _pz_par1_max << std::endl;

            px = _px_par1_min + (_px_par1_max - _px_par1_min) * gsl_rng_uniform_pos(RandomGenerator);
            py = _py_par1_min + (_py_par1_max - _py_par1_min) * gsl_rng_uniform_pos(RandomGenerator);
            pz = _pz_par1_min + (_pz_par1_max - _pz_par1_min) * gsl_rng_uniform_pos(RandomGenerator);
            muon1.SetXYZM(px, py, pz, MUON_MASS);

            // pT cut
            bool pass_pt = (muon1.Pt() >= _pt_min && muon1.Pt() <= _pt_max);

            // Bend-plane cut for mu+ (pid == -13)
            bool pass_bend = true;
            if (_apply_bend_cut) {
                double bend = px / pz;
                pass_bend = (bend >= _bend_lo_pos && bend <= _bend_hi_pos);
            }

            if (verbosity > 0)
                std::cout << "par1 (px, py, pz): (" << px << ", " << py << ", " << pz
                          << "), Pt: " << muon1.Pt()
                          << ", bend: " << px/pz
                          << ", pass_pt: " << pass_pt
                          << ", pass_bend: " << pass_bend << std::endl;

            if (pass_pt && pass_bend) break;
            iteration++;
        } while (iteration < max_iterations);

        if (iteration >= max_iterations) {
            std::cout << PHWHERE << "::Warning - Max iterations reached for muon 1"
                      << ", Pt=" << muon1.Pt()
                      << ", bend=" << px/pz << std::endl;
            return false;
        }

    } else if (muon_counter == 2) {
        // ----------------------------------------------------------------
        // Generate mu- (pdgcode == 13)
        // Cuts applied per iteration:
        //   1. pT in [_pt_min, _pt_max]
        //   2. Opening angle between mu+ and mu- < _theta_max
        //   3. xF < 1.0
        //   4. Bend-plane ratio bend = px/pz in [_bend_lo_neg, _bend_hi_neg]  (if enabled)
        //   5. Dimuon pz (mu+ pz + mu- pz) < BEAM_ENERGY                      (if enabled)
        // ----------------------------------------------------------------
        if (verbosity > 0)
            std::cout << "Mu- Vertex: (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ")" << std::endl;

        // Pre-compute CMS boost (does not depend on mu- kinematics)
        TLorentzVector p_beam(0., 0., sqrt(BEAM_ENERGY * BEAM_ENERGY - PROTON_MASS * PROTON_MASS), BEAM_ENERGY);
        TLorentzVector p_target(0., 0., 0., PROTON_MASS);
        TLorentzVector p_cms  = p_beam + p_target;
        double s              = p_cms.M2();
        TVector3 bv_cms       = p_cms.BoostVector();

        do {
            if (verbosity > 0) {
                std::cout << "loop 2 pdgcode: " << pdgcode << std::endl;
                std::cout << "_pz_par2_min: " << _pz_par2_min
                          << ", _pz_par2_max: " << _pz_par2_max << std::endl;
                std::cout << "Vertex: (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ")" << std::endl;
            }

            px = _px_par2_min + (_px_par2_max - _px_par2_min) * gsl_rng_uniform_pos(RandomGenerator);
            py = _py_par2_min + (_py_par2_max - _py_par2_min) * gsl_rng_uniform_pos(RandomGenerator);
            pz = _pz_par2_min + (_pz_par2_max - _pz_par2_min) * gsl_rng_uniform_pos(RandomGenerator);
            muon2.SetXYZM(px, py, pz, MUON_MASS);

            // pT cut
            bool pass_pt = (muon2.Pt() >= _pt_min && muon2.Pt() <= _pt_max);

            // Opening-angle cut
            angle = muon1.Vect().Angle(muon2.Vect()) * (180.0 / M_PI);
            bool pass_angle = (angle < _theta_max);

            // xF cut
            TLorentzVector p_sum = muon1 + muon2;
            double mass = p_sum.M();
            TLorentzVector p_sum_cms = p_sum;
            p_sum_cms.Boost(-bv_cms);
            xF = 2. * p_sum_cms.Pz() / TMath::Sqrt(s) / (1. - mass * mass / s);
            bool pass_xF = (xF < 1.0);

            // Bend-plane cut for mu- (pid == 13)
            bool pass_bend = true;
            if (_apply_bend_cut) {
                double bend = px / pz;
                pass_bend = (bend >= _bend_lo_neg && bend <= _bend_hi_neg);
            }

            // Dimuon pz < BEAM_ENERGY cut
            // Uses the scalar sum of the two muon pz values as a proxy for
            // the dimuon longitudinal momentum in the lab frame.
            bool pass_dimuon_pz = true;
            if (_apply_dimuon_pz_cut) {
                double dimuon_pz = muon1.Pz() + pz;   // mu+ pz + mu- pz
                pass_dimuon_pz = (dimuon_pz < BEAM_ENERGY);
            }

            if (verbosity > 0) {
                std::cout << "par2 (px, py, pz): (" << px << ", " << py << ", " << pz << ")"
                          << ", Pt: "        << muon2.Pt()
                          << ", angle: "     << angle
                          << ", xF: "        << xF
                          << ", bend: "      << px/pz
                          << ", dimuon_pz: " << muon1.Pz() + pz
                          << ", pass_pt: "        << pass_pt
                          << ", pass_angle: "     << pass_angle
                          << ", pass_xF: "        << pass_xF
                          << ", pass_bend: "      << pass_bend
                          << ", pass_dimuon_pz: " << pass_dimuon_pz << std::endl;
            }

            if (pass_pt && pass_angle && pass_xF && pass_bend && pass_dimuon_pz) break;
            iteration++;
        } while (iteration < max_iterations);

        if (iteration >= max_iterations) {
            std::cout << PHWHERE << "::Warning - Max iterations reached for muon 2"
                      << ", Pt="        << muon2.Pt()
                      << ", angle="     << angle
                      << ", xF="        << xF
                      << ", bend="      << px/pz
                      << ", dimuon_pz=" << muon1.Pz() + pz << std::endl;
            return false;
        }
    }
    return true;
}

bool PHG4E1039TrackPairGen::GenerateMomentumWithExpPDF(int muon_counter, double &px, double &py, double &pz) {
    if (!_hMomentumMap) {
        std::cout << PHWHERE << "::Error - _hMomentumMap is null" << std::endl;
        return false;
    }

    Double_t point[6];
    _hMomentumMap->GetRandom(point);

    if (muon_counter == 1) {
        px = point[0];  // mu_pos_px
        py = point[1];  // mu_pos_py
        pz = point[2];  // mu_pos_pz
    } else if (muon_counter == 2) {
        px = point[3];  // mu_neg_px
        py = point[4];  // mu_neg_py
        pz = point[5];  // mu_neg_pz
    }

    if (verbosity > 0) {
        std::cout << (muon_counter == 1 ? "mu+" : "mu-") << " (px, py, pz): ("
                  << px << ", " << py << ", " << pz << ")" << std::endl;
    }

    return true;
}

int PHG4E1039TrackPairGen::End(PHCompositeNode *topNode)
{
  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("PHG4SEG_EVENT_COUNT", _eventcount); // Without "Name()" here to make this variable unique

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4E1039TrackPairGen::set_pxpypz_range(const double x_min,
		const double x_max, const double y_min, const double y_max,
		const double z_min, const double z_max) {
	_px_min = x_min; _px_max = x_max;
	_py_min = y_min; _py_max = y_max;
	_pz_min = z_min; _pz_max = z_max;
}

void PHG4E1039TrackPairGen::set_par1_pxpypz_range(const double x_min,
                const double x_max, const double y_min, const double y_max,
                const double z_min, const double z_max) {
        _px_par1_min = x_min; _px_par1_max = x_max;
        _py_par1_min = y_min; _py_par1_max = y_max;
        _pz_par1_min = z_min; _pz_par1_max = z_max;
}

void PHG4E1039TrackPairGen::set_par2_pxpypz_range(const double x_min,
                const double x_max, const double y_min, const double y_max,
                const double z_min, const double z_max) {
        _px_par2_min = x_min; _px_par2_max = x_max;
        _py_par2_min = y_min; _py_par2_max = y_max;
        _pz_par2_min = z_min; _pz_par2_max = z_max;
}


void PHG4E1039TrackPairGen::set_max_opening_angle(const double x_max){
         _theta_max = x_max;
}

double
PHG4E1039TrackPairGen::smearvtx(const double position, const double width, FUNCTION dist) const
{
  double res = position;
  if (dist == Uniform)
    {
      res = (position-width) + 2*gsl_rng_uniform_pos(RandomGenerator)*width;
    }
  else if (dist == Gaus)
    {
      res = position + gsl_ran_gaussian(RandomGenerator,width);
    }
  return res;
}
