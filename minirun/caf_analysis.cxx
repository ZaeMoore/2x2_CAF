/*
Analysis of MiniRun CAF files

Zae Moore

Can be run on DUNE GPVMs or NERSC
On NERSC:
Run the alias cafshifter, then source setup.sh, then source compile.sh
Then run using command ./minirun_analysis
Compiled script cannot be run without being in the cafshifter and sourcing setup.sh first

Input: None

Output: A ROOT file with information on all events that pass the
reco and/or truth level cuts
*/

#include <iostream>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TVector.h>
#include <TVector3.h>
#include "TEfficiency.h"
#include "TMath.h"
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <chrono>
#include <set>
#include <algorithm>
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "FluxCovarianceReweight.h"
#define Dimension 3
//#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_05_00/include/duneanaobj/StandardRecord/StandardRecord.h"
//#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include/duneanaobj/StandardRecord/Proxy/SRProxy.h"

/*
Dot product function between 2 vectors
*/
float dot_product(std::vector<float> vector_a, std::vector<float> vector_b) 
{
  float product = 0;
  for (int i = 0; i < Dimension; i++)
    product = product + vector_a[i] * vector_b[i];
  return product;
}

/*
Is the point vertex contained within the TPC fiducial volume? 
Fiducial volume: each module's active region shrunk by `margin` cm on every side.
Regions are {{xmin, ymin, zmin}, {xmax, ymax, zmax}}.
*/
static const double kActiveRegions[4][2][3] = {
    {{  3.069, -62.076,   2.4620001}, { 63.931, 62.076,  64.538   }},
    {{  3.069, -62.076, -64.538    }, { 63.931, 62.076,  -2.4620001}},
    {{-63.931, -62.076,   2.462    }, { -3.069, 62.076,  64.537994 }},
    {{-63.931, -62.076, -64.538    }, { -3.069, 62.076,  -2.4620001}},
};

inline bool is_contained(std::vector<double> point, double margin = 5.0)
{
    double x = point[0];
    double y = point[1];
    double z = point[2];

    for (const auto& r : kActiveRegions)
    {
        if (r[0][0] + margin <= x && x <= r[1][0] - margin &&
            r[0][1] + margin <= y && y <= r[1][1] - margin &&
            r[0][2] + margin <= z && z <= r[1][2] - margin)
            return true;
    }
    return false;
}

// Straight-line length between a track's start and end points (cm).
// Returns -999 if either point is NaN, so it never passes a length threshold.
template <typename P>
double dist(const P& start, const P& end)
{
    const double dx = double(end.x) - double(start.x);
    const double dy = double(end.y) - double(start.y);
    const double dz = double(end.z) - double(start.z);
    const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    return std::isnan(len) ? -999.0 : len;
}

// Distance between two points given as {x, y, z}; huge if either has a NaN
inline double dist3(const std::vector<double>& a, const std::vector<double>& b)
{
    const double dx = b[0] - a[0];
    const double dy = b[1] - a[1];
    const double dz = b[2] - a[2];
    const double d = std::sqrt(dx*dx + dy*dy + dz*dz);
    return std::isnan(d) ? 1e9 : d;
}

/*
Truth level cuts
Only primary particles are passed to this function, but all particles are saved in the output tree
*/
bool truth_cuts(std::vector<double> vtx,
                std::vector<std::vector<double>> track_start,
                std::vector<double> track_length,
                std::vector<int> pdg)
{
    int nproton = 0;
    int nmuon = 0;
    int npion = 0;

    for(unsigned long i = 0; i < pdg.size(); ++i)
    {
        double vtx_dist = abs(dist3(vtx, track_start[i]));
        if(pdg[i] == 2212)
        {
            if (is_contained(track_start[i]) == true && vtx_dist < 5.0)
                nproton++;
        }
        else if(pdg[i] == 13 || pdg[i] == -13)
        {
            if (is_contained(track_start[i]) == true && vtx_dist < 5.0)
                nmuon++;
        }
        else if((pdg[i] == 111 || pdg[i] == 211 || pdg[i] == -211))
        {
            if (is_contained(track_start[i]) == true && vtx_dist < 5.0)
                npion++;
        }
    }

    if((nproton >= 2) && (nmuon == 1) && (npion == 0) && is_contained(vtx) == true)
    {
        return true;
    }
    return false;
}
/*
Reco level cuts 
Applied in line 561
*/
bool reco_cuts(std::vector<double> vtx,
                std::vector<int> pdg,
                std::vector<double> p,
                std::vector<bool> prim,
                std::vector<std::vector<double>> track_start,
                std::vector<std::vector<double>> track_end,
                std::vector<double> track_length,
                double mu_minerva_track_end_z)
{
    int nproton = 0;
    int nmuon = 0;
    int npion = 0;
    // track_end = {{x, y, z}, {x, y, z}, ...}

    for(unsigned long i = 0; i < pdg.size(); ++i)
    {
        const double vtx_dist = std::min(dist3(vtx, track_start[i]), dist3(vtx, track_end[i]));
        if(pdg[i] == 2212)
        {
            if(is_contained(track_end[i]) == true && vtx_dist < 5.0)
                nproton++;
        }
        else if((pdg[i] == 13 || pdg[i] == -13))
        {
            if(mu_minerva_track_end_z > 62.0 && vtx_dist < 5.0 && is_contained(track_start[i]) == true)
                nmuon++;
        }
        else if((pdg[i] == 111 || pdg[i] == 211 || pdg[i] == -211))
        {
            if(is_contained(track_start[i]) == true && vtx_dist < 5.0)
                npion++;
        }
    }

    if((nproton >= 2) && (nmuon == 1) && (npion == 0) && is_contained(vtx) == true)
    {
        return true;
    }
    return false;

}

double get_efficiency(int num_signal, int num_truth)
{
    if(num_truth == 0)
        return 0.0;

    return static_cast<double>(num_signal) / static_cast<double>(num_truth);
}

double get_purity(int num_signal, int num_reco)
{
    if(num_reco == 0)
        return 0.0;

    return static_cast<double>(num_signal) / static_cast<double>(num_reco);
}

const double kFluxWMin  = 0.01;
const double kFluxWMax  = 20.0;
const int    kNFluxBins = 628;   // matches FluxCovarianceReweight

std::vector<double> get_flux_weights(FluxCovarianceReweight& rw, int nu_pdg, double nu_E,
                                     int& bin_out, int nUniv)
{
    std::vector<double> w(nUniv, 1.0);
    bin_out = rw.FindBinIndex(nu_pdg, nu_E);   // RHC bins only (hardcoded in the class)
    if (bin_out < 0 || bin_out >= kNFluxBins)
        return w;

    std::vector<double> raw(nUniv);
    rw.FillAllThrowsForBin(bin_out, raw.data());
    for (int u = 0; u < nUniv; ++u)
        if (raw[u] > kFluxWMin && raw[u] < kFluxWMax)
            w[u] = raw[u];
    return w;
}

/*
Define interaction data
This will be called at the beginning of the caf plotter function to define these variables for the output tree
Then will be re-called at the beginning of every interaction loop to reset the variables for each interaction
Output root file will have a tree with one entry = one interaction
*/
struct InteractionData
{
    // Particle level reco variables
    std::vector< bool >    reco_prim;
    std::vector< double >  reco_energy;
    std::vector< double >  reco_p_x; 
    std::vector< double >  reco_p_y; 
    std::vector< double >  reco_p_z;
    std::vector< double >  reco_p_mag;
    std::vector< double >  reco_length;
    std::vector< double >  reco_angle;
    std::vector< double >  reco_angle_rot;
    std::vector< double >  reco_angle_incl;
    std::vector< double >  reco_angle_x;
    std::vector< double >  reco_angle_y;
    std::vector< double >  reco_angle_z;
    std::vector< double >  reco_track_start_x;
    std::vector< double >  reco_track_start_y;
    std::vector< double >  reco_track_start_z;
    std::vector< double >  reco_track_end_x;
    std::vector< double >  reco_track_end_y;
    std::vector< double >  reco_track_end_z;
    std::vector< int >     reco_pdg;
    std::vector< double >  reco_enu_calo;

    // Interaction level reco variables
    double reco_vtx_x;
    double reco_vtx_y;
    double reco_vtx_z;

    // Particle level truth variables
    std::vector< bool >    true_prim;
    std::vector< double >  true_energy;
    std::vector< double >  true_p_x; 
    std::vector< double >  true_p_y; 
    std::vector< double >  true_p_z;
    std::vector< double >  true_p_mag;
    std::vector< double >  true_length;
    std::vector< double >  true_angle;
    std::vector< double >  true_angle_rot;
    std::vector< double >  true_angle_incl;
    std::vector< double >  true_angle_x;
    std::vector< double >  true_angle_y;
    std::vector< double >  true_angle_z;
    std::vector< double >  true_track_start_x;
    std::vector< double >  true_track_start_y;
    std::vector< double >  true_track_start_z;
    std::vector< double >  true_track_end_x;
    std::vector< double >  true_track_end_y;
    std::vector< double >  true_track_end_z;
    std::vector< int >     true_pdg;

    // Interaction level truth variables
    double true_vtx_x;
    double true_vtx_y;
    double true_vtx_z;
    int mode;
    double nu_momentum_x;
    double nu_momentum_y;
    double nu_momentum_z;

    // Minerva track variables
    std::vector< double >  minerva_track_E;
    std::vector< double >  minerva_track_dir_x;
    std::vector< double >  minerva_track_dir_y;
    std::vector< double >  minerva_track_dir_z;
    std::vector< double >  minerva_track_enddir_x;
    std::vector< double >  minerva_track_enddir_y;
    std::vector< double >  minerva_track_enddir_z;
    std::vector< double >  minerva_track_start_x;
    std::vector< double >  minerva_track_start_y;
    std::vector< double >  minerva_track_start_z;
    std::vector< double >  minerva_track_end_x;
    std::vector< double >  minerva_track_end_y;
    std::vector< double >  minerva_track_end_z;
    std::vector< double >  minerva_track_len_cm;
    
    // Interaction level variables
    int events_count;
    double true_ixn_index;
    double reco_ixn_index;
    int spill_index;
    int file_index;
    int genie_index;
    int event;
    int run;
    int subrun;

    // Systematics (genie and flux)
    std::vector<double> genie_weights;   // 100 weights per saved interaction
    std::vector<double> NormCCMEC; // Normalization weight for CC MEC events
    std::vector<double> MaCCRES;
    std::vector<double> MvCCRES;
    std::vector<double> FSIReweight_hN;
    std::vector<double> FSIReweight_INCL;
    std::vector<double> FSIReweight_G4BC;
    std::vector<float> flux_weights;   // nUniv weights per saved interaction
    int flux_bin = -1;                 // -1 = outside the binning
};

/*
Main function to loop through CAF files
*/
int caf_plotter(bool is_flat = true)
{
    gErrorIgnoreLevel = kError; 
    BkgWhy bkg_why;
    
    // Create output file
    std::string file_name = "multip_analysis_m6.5_1.5";

    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");

    // Create the data object
    InteractionData data;

    const int nUniv = 1000;
    std::vector<double> truth_pass_u(nUniv, 0.0), reco_pass_u(nUniv, 0.0), signal_u(nUniv, 0.0);

    // Single entry variables for the output tree
    double purity;
    double efficiency;
    int num_reco_pass = 0; // Number of reco interactions that pass cuts
    int num_background = 0; // Pass reco but not truth
    int num_signal = 0; // Pass reco and truth
    int num_truth_pass = 0; // Number of truth interactions that pass cuts, regardless of if they have a matched reco interaction
    int num_signal_unique = 0;

    // Flux systematics outputs: declare before the branches that point at them
    double Phi_nom = 0;
    std::vector<double> Phi_uValues;
    
    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fCafTree=new TTree("CafTree", "Caf reco and truth variables");
    fCafTree->Branch("reco_prim", &data.reco_prim);
    fCafTree->Branch("reco_energy", &data.reco_energy);
    fCafTree->Branch("reco_p_x", &data.reco_p_x);
    fCafTree->Branch("reco_p_y", &data.reco_p_y);
    fCafTree->Branch("reco_p_z", &data.reco_p_z);
    fCafTree->Branch("reco_p_mag", &data.reco_p_mag);
    fCafTree->Branch("reco_length", &data.reco_length);
    fCafTree->Branch("reco_angle", &data.reco_angle);
    fCafTree->Branch("reco_angle_rot", &data.reco_angle_rot);
    fCafTree->Branch("reco_angle_incl", &data.reco_angle_incl);
    fCafTree->Branch("reco_angle_x", &data.reco_angle_x);
    fCafTree->Branch("reco_angle_y", &data.reco_angle_y);
    fCafTree->Branch("reco_angle_z", &data.reco_angle_z);
    fCafTree->Branch("reco_vtx_x", &data.reco_vtx_x);
    fCafTree->Branch("reco_vtx_y", &data.reco_vtx_y);
    fCafTree->Branch("reco_vtx_z", &data.reco_vtx_z);
    fCafTree->Branch("reco_track_start_x", &data.reco_track_start_x);
    fCafTree->Branch("reco_track_start_y", &data.reco_track_start_y);
    fCafTree->Branch("reco_track_start_z", &data.reco_track_start_z);
    fCafTree->Branch("reco_track_end_x", &data.reco_track_end_x);
    fCafTree->Branch("reco_track_end_y", &data.reco_track_end_y);
    fCafTree->Branch("reco_track_end_z", &data.reco_track_end_z);
    fCafTree->Branch("reco_pdg", &data.reco_pdg);
    fCafTree->Branch("reco_ixn_index", &data.reco_ixn_index);
    fCafTree->Branch("reco_enu_calo", &data.reco_enu_calo);

    fCafTree->Branch("true_prim", &data.true_prim);
    fCafTree->Branch("true_energy", &data.true_energy);
    fCafTree->Branch("true_p_x", &data.true_p_x);
    fCafTree->Branch("true_p_y", &data.true_p_y);
    fCafTree->Branch("true_p_z", &data.true_p_z);
    fCafTree->Branch("true_p_mag", &data.true_p_mag);
    fCafTree->Branch("true_length", &data.true_length);
    fCafTree->Branch("true_angle", &data.true_angle);
    fCafTree->Branch("true_angle_rot", &data.true_angle_rot);
    fCafTree->Branch("true_angle_incl", &data.true_angle_incl);
    fCafTree->Branch("true_angle_x", &data.true_angle_x);
    fCafTree->Branch("true_angle_y", &data.true_angle_y);
    fCafTree->Branch("true_angle_z", &data.true_angle_z);
    fCafTree->Branch("true_vtx_x", &data.true_vtx_x);
    fCafTree->Branch("true_vtx_y", &data.true_vtx_y);
    fCafTree->Branch("true_vtx_z", &data.true_vtx_z);
    fCafTree->Branch("true_track_start_x", &data.true_track_start_x);
    fCafTree->Branch("true_track_start_y", &data.true_track_start_y);
    fCafTree->Branch("true_track_start_z", &data.true_track_start_z);
    fCafTree->Branch("true_track_end_x", &data.true_track_end_x);
    fCafTree->Branch("true_track_end_y", &data.true_track_end_y);
    fCafTree->Branch("true_track_end_z", &data.true_track_end_z);
    fCafTree->Branch("true_pdg", &data.true_pdg);
    fCafTree->Branch("true_ixn_index", &data.true_ixn_index);
    fCafTree->Branch("mode", &data.mode);
    fCafTree->Branch("nu_momentum_x", &data.nu_momentum_x);
    fCafTree->Branch("nu_momentum_y", &data.nu_momentum_y);
    fCafTree->Branch("nu_momentum_z", &data.nu_momentum_z);

    fCafTree->Branch("minerva_track_E", &data.minerva_track_E);
    fCafTree->Branch("minerva_track_dir_x", &data.minerva_track_dir_x);
    fCafTree->Branch("minerva_track_dir_y", &data.minerva_track_dir_y);
    fCafTree->Branch("minerva_track_dir_z", &data.minerva_track_dir_z);
    fCafTree->Branch("minerva_track_enddir_x", &data.minerva_track_enddir_x);
    fCafTree->Branch("minerva_track_enddir_y", &data.minerva_track_enddir_y);
    fCafTree->Branch("minerva_track_enddir_z", &data.minerva_track_enddir_z);
    fCafTree->Branch("minerva_track_start_x", &data.minerva_track_start_x);
    fCafTree->Branch("minerva_track_start_y", &data.minerva_track_start_y);
    fCafTree->Branch("minerva_track_start_z", &data.minerva_track_start_z);
    fCafTree->Branch("minerva_track_end_x", &data.minerva_track_end_x);
    fCafTree->Branch("minerva_track_end_y", &data.minerva_track_end_y);
    fCafTree->Branch("minerva_track_end_z", &data.minerva_track_end_z);
    fCafTree->Branch("minerva_track_len_cm", &data.minerva_track_len_cm);

    fCafTree->Branch("genie_weights", &data.genie_weights);
    fCafTree->Branch("NormCCMEC", &data.NormCCMEC);
    fCafTree->Branch("MaCCRES", &data.MaCCRES);
    fCafTree->Branch("MvCCRES", &data.MvCCRES);
    fCafTree->Branch("FSIReweight_hN", &data.FSIReweight_hN);
    fCafTree->Branch("FSIReweight_INCL", &data.FSIReweight_INCL);
    fCafTree->Branch("FSIReweight_G4BC", &data.FSIReweight_G4BC);
    fCafTree->Branch("flux_bin", &data.flux_bin);
    fCafTree->Branch("flux_weights", &data.flux_weights);

    fCafTree->Branch("events_count", &data.events_count);
    fCafTree->Branch("spill_index", &data.spill_index);
    fCafTree->Branch("file_index", &data.file_index);
    fCafTree->Branch("genie_index", &data.genie_index);
    fCafTree->Branch("event", &data.event);
    fCafTree->Branch("run", &data.run);
    fCafTree->Branch("subrun", &data.subrun);

    // Single variable outputs
    TTree *fEffTree=new TTree("EffTree", "Caf efficiency variables");
    fEffTree->Branch("purity", &purity);
    fEffTree->Branch("efficiency", &efficiency);
    fEffTree->Branch("num_reco_pass", &num_reco_pass);
    fEffTree->Branch("num_background", &num_background);
    fEffTree->Branch("num_signal", &num_signal);
    fEffTree->Branch("num_truth_pass", &num_truth_pass);
    fEffTree->Branch("num_signal_unique", &num_signal_unique);

    // Flux systematics information
    TTree *fFluxTree=new TTree("FluxTree", "Caf flux systematics variables");
    fFluxTree->Branch("Phi_nom", &Phi_nom);
    fFluxTree->Branch("Phi_uValues", &Phi_uValues);

    // Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    // z-direction (roughly beam dir)
    const auto z_plus_dir = TVector3(0, 0, 1.0);
    const auto y_plus_dir = TVector3(0, 1.0, 0.0);
    const auto x_plus_dir = TVector3(1.0, 0, 0.0);

    // negative y-direction 
    const auto y_minus_dir = TVector3(0, -1.0, 0.0);

    double minTrkLength = 3;

    int num_events = 0;

    //const int nUniv = 1000;

    /*
    Flux systematics information
    */

    // Needs numiFluxSystFull.root and flux_covariance_binning_NuMI_GeV.txt in the cwd
    FluxCovarianceReweight rw;

    if (!rw.GenerateThrows(nUniv, 12345) || rw.GetIntegratedFluxNominal() <= 0)
    {
        std::cerr << "Flux setup failed (check input files in cwd)\n";
        return 1;
    }

    Phi_nom = rw.GetIntegratedFluxNominal();   // no "double": assign to the outer variable
    const auto& phis = rw.GetIntegratedFluxThrows();
    if ((int)phis.size() != nUniv)
    {
        std::cerr << "Bad Phi_u size: " << phis.size() << "\n";
        return 1;
    }
    Phi_uValues.assign(phis.begin(), phis.end());
    fFluxTree->Fill();                         // one entry

    /*
    ----------------------------------------------------------------------------------------
    Loop through files in the CAF file list
    Save relevant information for events that pass truth or reco level cuts
    ----------------------------------------------------------------------------------------
    */
    const auto t_start{std::chrono::steady_clock::now()};
    for(unsigned long file_num = 0; file_num < 1000; ++file_num)
    {
        std::cout << "Processing file_num " << file_num << " / 1000" << std::endl;

        std::string file_path = "/global/cfs/cdirs/dune/www/data/2x2/simulation/productions";

        std::string caf_path = file_path + "/MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF.flat/0000000/" +
            Form("MiniRun6.5_1E19_RHC.caf.%07lu.CAF.flat.root", file_num);

        // Open file and attach SRProxy Object
        TFile* caf_file = TFile::Open(caf_path.c_str(), "READ");
        
        if(!caf_file || caf_file->IsZombie())
        {
            std::cerr << "Error opening CAF file: " << caf_path << std::endl;
            continue;
        }
        
        TTree* caf_tree = (TTree*)caf_file->Get("cafTree");

        if(!caf_tree)
        {
            std::cerr << "Error: cafTree not found in file: " << caf_path << std::endl;
            caf_file->Close();
            continue;
        }

        std::string tree_name = is_flat ? "rec" : "";
        auto sr = new caf::SRProxy(caf_tree, tree_name);

        // Open the companion GENIE systematics file
        std::string genie_path = file_path + "/systematics/nusystematics/MiniRun6.5.nusyst/" +
            Form("MiniRun6.5_1E19_RHC.nuweights.%07lu.nusyst.root", file_num);

        TFile* genie_rw_file = TFile::Open(genie_path.c_str(), "READ");
        if (!genie_rw_file || genie_rw_file->IsZombie())
        {
            std::cerr << "Error opening GENIE file: " << genie_path << std::endl;
            delete sr;
            caf_file->Close();
            continue;
        }
        TTree* genie_rw_tree = (TTree*)genie_rw_file->Get("SystWeights");
        if (!genie_rw_tree)
        {
            std::cerr << "SystWeights not found in: " << genie_path << std::endl;
            delete sr;
            genie_rw_file->Close();
            caf_file->Close();
            continue;
        }
        Double_t totWeight[100];
        Double_t NormCCMEC[100];
        Double_t MaCCRES[100];
        Double_t MvCCRES[100];
        Double_t FSIReweight_hN[100];
        Double_t FSIReweight_INCL[100];
        Double_t FSIReweight_G4BC[100];
        genie_rw_tree->SetBranchAddress("totWeight", totWeight);
        genie_rw_tree->SetBranchAddress("NormCCMEC", NormCCMEC);
        genie_rw_tree->SetBranchAddress("MaCCRES", MaCCRES);
        genie_rw_tree->SetBranchAddress("MvCCRES", MvCCRES);
        genie_rw_tree->SetBranchAddress("FSIReweight_hN", FSIReweight_hN);
        genie_rw_tree->SetBranchAddress("FSIReweight_INCL", FSIReweight_INCL);
        genie_rw_tree->SetBranchAddress("FSIReweight_G4BC", FSIReweight_G4BC);

        /*
        ----------------------------------------------------------------------------------------
        Loop over each spill
        ----------------------------------------------------------------------------------------
        */
        const unsigned long nspills = caf_tree->GetEntries();
        const unsigned long incr = std::max(1UL, nspills / 10);
        for (unsigned long i = 0; i < nspills; ++i)
        {
            caf_tree->GetEntry(i);

            // Keep track of spill # with print statement
            if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

            int spill_num = i;

            const auto num_ixn = sr->common.ixn.ndlp;

            // Don't repeat truth interactions after they've been matched
            std::set<long> truth_counted, signal_counted;

            /*
            ----------------------------------------------------------------------------------------
            Loop over each reco interaction
            Check if the interaction passes cuts
            Primary particles are the focus of the cuts, but all particles are saved in output tree
            ----------------------------------------------------------------------------------------
            */
            for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
            {
                // This gives me a fresh set of vectors for each interaction
                // I don't have to worry about clearing them at the end of the loop
                data = InteractionData{};
                
                bool reco_passes = false;
                bool truth_passes = false;

                int partMult = 0;
                double longestTrk = -9999;
                int trackMult = 0;
                int trackMultExit = 0;

                // Reco interaction
                const auto& reco_ixn = sr->common.ixn.dlp[ixn];
                const auto& vtx = reco_ixn.vtx;

                // Get the truth interaction(s) corresponding to this reco interaction
                const auto& vec_truth_ixn = reco_ixn.truth;
                const auto& vec_overlap_ixn = reco_ixn.truthOverlap;

                if(vec_overlap_ixn.empty())
                    continue;

                // Find the truth interaction with the largest overlap
                double current_max = 0;
                unsigned int max_overlap = 0;
                for(unsigned int i = 0; i < vec_overlap_ixn.size(); i++)
                {
                    auto val = vec_overlap_ixn.at(i);
                    if(val > current_max)
                    {
                        current_max = val;
                        max_overlap = i;
                    }
                }

                // Matched truth interaction
                const auto truth_idx = vec_truth_ixn.at(max_overlap);
                const auto& truth_ixn = sr->mc.nu[truth_idx];

                std::vector<double> vtx_vec = {vtx.x, vtx.y, vtx.z};
                std::vector<double> truth_vtx_vec = {truth_ixn.vtx.x, truth_ixn.vtx.y, truth_ixn.vtx.z};

                // If target is not argon, skip interaction (vtx cut is in reco/truth cuts)
                if(truth_ixn.targetPDG != 1000180400)
                    continue;

                // Initialize vectors to hold reco particle information
                std::vector<int> reco_pdg;
                std::vector<double> reco_p;
                std::vector<bool> reco_prim;
                std::vector<std::vector<double>> reco_track_start;
                std::vector<std::vector<double>> reco_track_end;
                std::vector<double> reco_length;
                double mu_minerva_track_end_z = -9999;
                std::vector<int> true_pdg;
                std::vector<bool> true_prim;
                std::vector<std::vector<double>> truth_track_start;
                std::vector<double> true_length;

                // Get truth particle information for ALL truth particles, not just matching ones
                for(unsigned long ipart = 0; ipart < truth_ixn.prim.size(); ++ipart)
                    {
                        const auto& part = truth_ixn.prim[ipart];
                        true_pdg.push_back(part.pdg);
                        truth_track_start.push_back({part.start_pos.x, part.start_pos.y, part.start_pos.z});
                        true_length.push_back(dist(part.start_pos, part.end_pos));
                    }

                // Initialize information for the minerva track matching 
                // Reset for every interaction, fill for each particle
                const auto nparts = sr->common.ixn.dlp[ixn].part.dlp.size();
                std::vector<int> minerva_idx(nparts, -9999);
                std::vector<int> minerva_ixn(nparts, -9999);
                std::vector<bool> minerva_match(nparts, false);

                /*
                ----------------------------------------------------------------------------------------
                Loop over each particle in the reco interaction
                Check if the event passes reco cuts
                ----------------------------------------------------------------------------------------
                */
                for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
                {
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                    int pdg = part.pdg;
                    reco_pdg.push_back(pdg);
                    reco_p.push_back(TVector3(part.p.x, part.p.y, part.p.z).Mag());
                    reco_track_start.push_back({part.start.x, part.start.y, part.start.z});
                    reco_track_end.push_back({part.end.x, part.end.y, part.end.z});
                    reco_prim.push_back(part.primary);
                    reco_length.push_back(dist(part.start, part.end));
                    // Initialized per particle
                    double best_ang = -1;
                    int ixnM = -1;
                    int idxM = -1; 
        
                    // Get Minerva match
                    bool minerva_track = false;
                    // Loop over tracks
                    if ((abs(pdg) == 2212 || abs(pdg) == 13 || abs(pdg) == 211 || 
                        abs(pdg) == 111 || abs(pdg) == 321))
                    {
                        const auto& start_pos = part.start;
                        const auto& end_pos = part.end;
                        double diffVertexdZ = abs(start_pos.z - sr->common.ixn.dlp[ixn].vtx.z);
                        double diffVertexdX = abs(start_pos.x - sr->common.ixn.dlp[ixn].vtx.x);
                        double diffVertexdY = abs(start_pos.y - sr->common.ixn.dlp[ixn].vtx.y);
                        double diffVertex = TMath::Sqrt(diffVertexdZ * diffVertexdZ + 
                                                        diffVertexdX * diffVertexdX + 
                                                        diffVertexdY*diffVertexdY);
                        if (diffVertex > 5.0)
                            continue;
                        // Make sure it is near the vertex
                        double dX = (end_pos.x - start_pos.x);
                        double dY = (end_pos.y - start_pos.y);
                        double dZ = (end_pos.z - start_pos.z);
                        double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
                        // double dirX = dX / length;
                        // double dirY = dY / length;
                        // double dirZ = dZ / length;

                        if (std::isnan(start_pos.z))
                            length = -999;
                        if (length > longestTrk)
                            longestTrk = length;
                        // Make sure it is above the track threshold
                        if (part.primary == true && length > minTrkLength)
                        {
                            partMult++;
                            trackMult++;
                            // See if it punches out and match it to MINERvA
                            if ((start_pos.z) > 62 || (end_pos.z) > 62)
                                trackMultExit++;
                            if ((abs(start_pos.z) > 62 || abs(end_pos.z) > 62))
                            {
                                for(long unsigned int k=0; k<sr->nd.trkmatch.extrap.size(); k++)
                                {
                                    if (sr->nd.trkmatch.extrap[k].larid.ixn != static_cast<int>(ixn) || 
                                        sr->nd.trkmatch.extrap[k].larid.reco != 1) continue;
                                    int index=sr->nd.trkmatch.extrap[k].larid.idx;
                                    if (sr->nd.lar.dlp[ixn].tracks[index].start.z==sr->common.ixn.dlp[ixn].part.dlp[ipart].start.z && 
                                        sr->nd.lar.dlp[ixn].tracks[index].end.z==sr->common.ixn.dlp[ixn].part.dlp[ipart].end.z)
                                    {
                                        double ang = std::abs(sr->nd.trkmatch.extrap[k].angdispl);
                                        if (ang < 0.99 || ang <= best_ang) continue;
                                        best_ang = ang;
                                        ixnM = sr->nd.trkmatch.extrap[k].minervaid.ixn;
                                        idxM = sr->nd.trkmatch.extrap[k].minervaid.idx;
                                        minerva_track = true;
                                    }
                                }
                            }
                        }
                    } // End of Minerva matching

                    if(minerva_track)
                    {
                        if(pdg == 13 || pdg == -13)
                            mu_minerva_track_end_z = sr->nd.minerva.ixn[ixnM].tracks[idxM].end.z;

                        minerva_idx[ipart]   = idxM;
                        minerva_ixn[ipart]   = ixnM;
                        minerva_match[ipart] = true;
                    }
                } // End of reco particle loop

                // ----------------------------------------------------------------------------------------
                // ----------------------------------------------------------------------------------------
                // Count the number of events that pass 
                // 1. Reco cuts
                // 2. Truth cuts
                // 3. Signal: Pass both reco and truth cuts
                // 4. Background: Pass reco but not truth cuts
                // ----------------------------------------------------------------------------------------
                // ----------------------------------------------------------------------------------------

                reco_passes = reco_cuts(vtx_vec, reco_pdg, reco_p, reco_prim, reco_track_start, reco_track_end, reco_length, mu_minerva_track_end_z);
                num_reco_pass += reco_passes ? 1 : 0;

                truth_passes = truth_cuts(truth_vtx_vec, truth_track_start, true_length, true_pdg);
                // num_truth_pass += truth_passes ? 1 : 0;

                const long tidx = static_cast<long>(truth_idx);
                if (truth_passes && truth_counted.insert(tidx).second) num_truth_pass++;

                if (reco_passes && truth_passes)
                {
                    num_signal++;                                                 // purity numerator
                    if (signal_counted.insert(tidx).second) num_signal_unique++;  // efficiency numerator
                }
                else if (reco_passes && !truth_passes)
                {
                    num_background++;
                    std::vector<bool> nc;
                    for (size_t k = 0; k < reco_pdg.size(); ++k)
                        nc.push_back(dist3(vtx_vec, reco_track_start[k]) < 5.0 && is_contained(reco_track_start[k]));
                    tally_bkg(bkg_why, int(truth_ixn.npi0), int(truth_ixn.npip) + int(truth_ixn.npim), reco_pdg, nc);
                }
                    

                // If interaction does not pass reco OR truth cuts, go to next interaction      
                if(!reco_passes && !truth_passes)
                    continue;

                // Flux weights for this interaction, for each universe
                int flux_bin = -1;
                const auto flux_w = get_flux_weights(rw, static_cast<int>(truth_ixn.pdg),
                                                    static_cast<double>(truth_ixn.E), flux_bin, nUniv);
                for (int u = 0; u < nUniv; ++u)
                {
                    if (reco_passes)                 reco_pass_u[u] += flux_w[u];
                    if (reco_passes && truth_passes) signal_u[u]    += flux_w[u];
                }
                data.flux_weights.assign(flux_w.begin(), flux_w.end());
                data.flux_bin = flux_bin;
                                
                // Interaction passes, now count it as an event
                num_events++;

                /*
                ----------------------------------------------------------------------------------------
                Loop over each particle in the reco interaction 
                Now to save information.
                Save primaries and secondaries
                ----------------------------------------------------------------------------------------
                */
                for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
                {
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];
                    
                    // Get truth particle matches for reco particle
                    caf::Proxy<caf::SRTrueParticle>* truth_match = nullptr;
                    const auto& vec_truth_id = part.truth;
                    const auto& vec_overlap = part.truthOverlap;
                    bool truth_primary = false;

                    // If the truth overlap vector is empty, assume no truth match and skip
                    if(vec_overlap.empty())
                        continue;

                    // Find the truth particle with the largest overlap
                    double current_max = 0;
                    unsigned int max_overlap = 0;
                    for(unsigned int i=0; i < vec_overlap.size(); i++)
                    {
                        auto val = vec_overlap.at(i);
                        if(val > current_max)
                        {
                            current_max = val;
                            max_overlap = i;
                        }
                    }

                    const auto& truth_id = vec_truth_id.at(max_overlap);

                    // Get pointer to the corresponding truth particle
                    if(truth_id.type == 1)
                    {
                        truth_match = &(sr->mc.nu[truth_id.ixn].prim[truth_id.part]);
                        truth_primary = true;
                    }
                    else if(truth_id.type == 3)
                        truth_match = &(sr->mc.nu[truth_id.ixn].sec[truth_id.part]);
                    else
                    {
                        std::cout << "Invalid truth id type!" << std::endl;
                        continue;
                    }        

                    // Genie systematic weights for this truth interaction
                    const int genie_idx = truth_ixn.genieIdx;
                    if (genie_idx >= 0 && genie_idx < genie_rw_tree->GetEntries())
                    {
                        genie_rw_tree->GetEntry(genie_idx);
                        data.genie_weights.assign(totWeight, totWeight + 100);
                        data.NormCCMEC.assign(NormCCMEC, NormCCMEC + 100);
                        data.MaCCRES.assign(MaCCRES, MaCCRES + 100);
                        data.MvCCRES.assign(MvCCRES, MvCCRES + 100);
                        data.FSIReweight_hN.assign(FSIReweight_hN, FSIReweight_hN + 100);
                        data.FSIReweight_INCL.assign(FSIReweight_INCL, FSIReweight_INCL + 100);
                        data.FSIReweight_G4BC.assign(FSIReweight_G4BC, FSIReweight_G4BC + 100);
                    }
                    data.genie_index = genie_idx;
                    
                    

                    // Get/calculate various reco/truth quantities
                    auto pvec = TVector3(part.p.x, part.p.y, part.p.z);
                    auto dir = TVector3(part.end.x, part.end.y, part.end.z) - 
                                TVector3(part.start.x, part.start.y, part.start.z);
                    //auto cos_angle = TMath::Cos(dir.Angle(beam_dir)); //Calculate cos of angle wrt neutrino beam direction
                    TVector3 dir_rot = dir;
                    dir_rot.RotateY(-TMath::Pi()/2);
                    double length = dir.Mag();

                    auto true_pvec = TVector3(truth_match->p.px, truth_match->p.py, truth_match->p.pz);
                    auto true_dir = TVector3(truth_match->end_pos.x, truth_match->end_pos.y, truth_match->end_pos.z) - 
                                    TVector3(truth_match->start_pos.x, truth_match->start_pos.y, truth_match->start_pos.z);
                    TVector3 true_dir_rot = true_dir;
                    true_dir_rot.RotateY(-TMath::Pi()/2);
                    double true_length_val = true_dir.Mag();

                    // Currently unused
                    //auto T_diff = truth_match->p.E - part.E;
                    //auto p_diff = true_pvec.Mag() - pvec.Mag();
                    //auto length_diff = true_length_val - length;
                    //auto cos_angle_diff = true_cos_angle - cos_angle;

                    dir.RotateY(-TMath::Pi()/2);
                    true_dir.RotateY(-TMath::Pi()/2);

                    // Population information in vectors for tracks that have passed all cuts
                    // Reco
                    data.reco_prim.push_back(part.primary);
                    data.reco_energy.push_back(part.E);
                    data.reco_p_x.push_back(part.p.x);
                    data.reco_p_y.push_back(part.p.y);
                    data.reco_p_z.push_back(part.p.z);
                    data.reco_p_mag.push_back(pvec.Mag());
                    data.reco_length.push_back(length);
                    data.reco_angle.push_back(dir.Angle(beam_dir));
                    data.reco_angle_x.push_back(dir.Angle(x_plus_dir));
                    data.reco_angle_y.push_back(dir.Angle(y_plus_dir));
                    data.reco_angle_z.push_back(dir.Angle(z_plus_dir));
                    data.reco_angle_rot.push_back(dir_rot.Theta());
                    data.reco_angle_incl.push_back(dir_rot.Phi());
                    data.reco_vtx_x = (sr->common.ixn.dlp[ixn].vtx.x);
                    data.reco_vtx_y = (sr->common.ixn.dlp[ixn].vtx.y);
                    data.reco_vtx_z = (sr->common.ixn.dlp[ixn].vtx.z);
                    data.reco_track_start_x.push_back(part.start.x);
                    data.reco_track_start_y.push_back(part.start.y);
                    data.reco_track_start_z.push_back(part.start.z);
                    data.reco_track_end_x.push_back(part.end.x);
                    data.reco_track_end_y.push_back(part.end.y);
                    data.reco_track_end_z.push_back(part.end.z);
                    data.reco_pdg.push_back(part.pdg);
                    data.reco_enu_calo.push_back(sr->common.ixn.dlp[ixn].Enu.calo);
                    // Truth
                    data.true_prim.push_back(truth_primary);
                    data.true_energy.push_back(truth_match->p.E);
                    data.true_p_x.push_back(truth_match->p.px); 
                    data.true_p_y.push_back(truth_match->p.py); 
                    data.true_p_z.push_back(truth_match->p.pz);
                    data.true_p_mag.push_back(true_pvec.Mag());
                    data.true_length.push_back(true_length_val);
                    data.true_angle.push_back(true_dir.Angle(beam_dir));
                    data.true_angle_x.push_back(true_dir.Angle(x_plus_dir));
                    data.true_angle_y.push_back(true_dir.Angle(y_plus_dir));
                    data.true_angle_z.push_back(true_dir.Angle(z_plus_dir));
                    data.true_angle_rot.push_back(true_dir_rot.Theta());
                    data.true_angle_incl.push_back(true_dir_rot.Phi());
                    data.true_vtx_x = truth_ixn.vtx.x;
                    data.true_vtx_y = truth_ixn.vtx.y;
                    data.true_vtx_z = truth_ixn.vtx.z;
                    data.true_track_start_x.push_back(truth_match->start_pos.x);
                    data.true_track_start_y.push_back(truth_match->start_pos.y);
                    data.true_track_start_z.push_back(truth_match->start_pos.z);
                    data.true_track_end_x.push_back(truth_match->end_pos.x);
                    data.true_track_end_y.push_back(truth_match->end_pos.y);
                    data.true_track_end_z.push_back(truth_match->end_pos.z);
                    data.true_pdg.push_back(truth_match->pdg);
                    data.mode = truth_ixn.mode;
                    data.nu_momentum_x = truth_ixn.momentum.x;
                    data.nu_momentum_y = truth_ixn.momentum.y;
                    data.nu_momentum_z = truth_ixn.momentum.z;

                    // Minerva
                    if (minerva_match[ipart] == true) // in a reco particle loop
                    {
                        // Saved in the previous reco particle loop
                        int idxM = minerva_idx[ipart];
                        int ixnM = minerva_ixn[ipart];

                        data.minerva_track_E.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].E);
                        data.minerva_track_dir_x.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].dir.x);
                        data.minerva_track_dir_y.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].dir.y);
                        data.minerva_track_dir_z.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].dir.z);
                        data.minerva_track_enddir_x.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].enddir.x);
                        data.minerva_track_enddir_y.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].enddir.y);
                        data.minerva_track_enddir_z.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].enddir.z);
                        data.minerva_track_start_x.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].start.x);
                        data.minerva_track_start_y.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].start.y);
                        data.minerva_track_start_z.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].start.z);
                        data.minerva_track_end_x.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].end.x);
                        data.minerva_track_end_y.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].end.y);
                        data.minerva_track_end_z.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].end.z);
                        data.minerva_track_len_cm.push_back(sr->nd.minerva.ixn[ixnM].tracks[idxM].len_cm);
                    }
                    else
                    {
                        data.minerva_track_E.push_back(-999);
                        data.minerva_track_dir_x.push_back(-999);
                        data.minerva_track_dir_y.push_back(-999);
                        data.minerva_track_dir_z.push_back(-999);
                        data.minerva_track_enddir_x.push_back(-999);
                        data.minerva_track_enddir_y.push_back(-999);
                        data.minerva_track_enddir_z.push_back(-999);
                        data.minerva_track_start_x.push_back(-999);
                        data.minerva_track_start_y.push_back(-999);
                        data.minerva_track_start_z.push_back(-999);
                        data.minerva_track_end_x.push_back(-999);
                        data.minerva_track_end_y.push_back(-999);
                        data.minerva_track_end_z.push_back(-999);
                        data.minerva_track_len_cm.push_back(-999);
                    }
                    // Other
                    data.events_count = num_events;
                    data.true_ixn_index = truth_idx;
                    data.reco_ixn_index = ixn;
                    data.spill_index = spill_num;
                    data.file_index = file_num;
                    data.event = sr->meta.nd_lar.event;
                    data.run = sr->meta.nd_lar.run;
                    data.subrun = sr->meta.nd_lar.subrun;

                } // End of particle loop

                // Fill at the end of every interaction. One entry = One interaction
                fCafTree->Fill();
            } // End of reco interaction loop

            // Truth interaction loop
            

        } // End of spill loop
        f->cd();
        fCafTree->AutoSave("Overwrite");
        delete sr;               // first, while the tree still exists
        genie_rw_file->Close();
        caf_file->Close();
    } // End of file loop

    
    
    /*
    ----------------------------------------------------------------------------------------
    Save flux throws, efficiency, purity to output file
    Pass those values to the output ROOT file for later use in the real data analysis
    TTree *fFluxTree=new TTree("FluxTree", "Caf flux systematics variables");
    fFluxTree->Branch("Phi_nom", &Phi_nom);
    fFluxTree->Branch("Phi_uValues", &Phi_uValues);
    ----------------------------------------------------------------------------------------
    */
    efficiency = get_efficiency(num_signal_unique, num_truth_pass);
    purity = get_purity(num_signal, num_reco_pass);
    fEffTree->Fill();

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    // POPULATE: Write to output ROOT file
    fCafTree->Write("", TObject::kOverwrite);
    fEffTree->Write("", TObject::kOverwrite);
    fFluxTree->Write("", TObject::kOverwrite);
        
    std::cout << "Wrote TTree." << std::endl;

    // CLOSE: Output ROOT file
    f->Close();

    std::cout << "\nBackground (reco pass, truth fail), n=" << bkg_why.n << ":\n"
          << "  true pi0 in event:                    " << bkg_why.pi0 << "\n"
          << "  charged pions only:                   " << bkg_why.chpi_only << "\n"
          << "  no true pion (other reason):          " << bkg_why.no_pion_truth << "\n"
          << "  reco pion that escaped the veto:      " << bkg_why.reco_pion_outside_veto << "\n"
          << "  reco photon/electron present:         " << bkg_why.reco_shower << "\n"
          << "  neither seen by reco (hard to cut):   " << bkg_why.no_reco_hint << "\n";

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;


    return 0;
}

int main()
{
    caf_plotter(true);

    return 0;
}