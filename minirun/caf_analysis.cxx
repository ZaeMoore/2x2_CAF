/*
Analysis of MiniRun CAF files

Zae Moore

To be run on DUNE GPVMs
This needs to be adapted to run on NERSC in order to use the systematics files Richie has there
Will have to remake the .list file

Input: A .list file contained the paths to every CAF file

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
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
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
Fiducial volume defined as 8 cm from each TPC wall
*/
bool contained(double x, double y, double z)
{
    double tpc_dist = 8.0; // distance from the tpc walls for containment cuts
    double xbound = 63.931;
    double ybound = 62.076;
    double zbound = 64.3163;

    bool cont = abs(x) < xbound - tpc_dist &&
                abs(x) > tpc_dist &&
                abs(y) < ybound - tpc_dist &&
                abs(z) > tpc_dist &&
                abs(z) < zbound - tpc_dist;
    return cont;
}

/*
Truth level cuts
Applied in line 420
*/
bool truth_cuts(int nproton, int nmuon, int npion)
{
    if((nproton >= 2) && (nmuon == 1) && (npion == 0))
    {
        return true;
    }
    return false;
}

/*
Reco level cuts 
Applied in line 395
*/
bool reco_cuts(int nproton, int nmuon, int npion)
{
    if((nproton >= 2) && (nmuon == 1) && (npion == 0))
    {
        return true;
    }
    return false;
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
    int true_nproton;
    int true_nmuon;
    int true_npion;
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

    // Systematics
    std::vector<double> genie_weights;   // 100 weights per saved interaction
};

/*
Main function to loop through CAF files
*/
int caf_plotter(bool is_flat = true)
{
    // Create output file
    std::string file_name = "multip_analysis_m6.5";

    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");

    // Create the data object
    InteractionData data;
    
    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fCafTree=new TTree("CafTree", "Caf reco and truth variables");
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
    fCafTree->Branch("true_nproton", &data.true_nproton);
    fCafTree->Branch("true_nmuon", &data.true_nmuon);
    fCafTree->Branch("true_npion", &data.true_npion);
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

    fCafTree->Branch("events_count", &data.events_count);
    fCafTree->Branch("spill_index", &data.spill_index);
    fCafTree->Branch("file_index", &data.file_index);
    fCafTree->Branch("genie_index", &data.genie_index);
    fCafTree->Branch("event", &data.event);
    fCafTree->Branch("run", &data.run);
    fCafTree->Branch("subrun", &data.subrun);

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

    // Loop through files
    const auto t_start{std::chrono::steady_clock::now()};
    for(unsigned long file_num = 0; file_num < 1000; ++file_num)
    {
        std::string file_path = "/global/cfs/cdirs/dune/www/data/2x2/simulation/productions";

        // Open file and attach SRProxy Object
        TFile* caf_file = TFile::Open(Form(file_path + "MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF.flat/0000000/MiniRun6.5_1E19_RHC.caf.%07d.CAF.flat.root", file_num), "READ");
        
        if(!caf_file || caf_file->IsZombie())
        {
            std::cerr << "Error opening CAF file: " << Form(file_path + "MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF.flat/0000000/MiniRun6.5_1E19_RHC.caf.%07d.CAF.flat.root", file_num) << std::endl;
            continue;
        }
        
        TTree* caf_tree = (TTree*)caf_file->Get("cafTree");

        if(!caf_tree)
        {
            std::cerr << "Error: cafTree not found in file: " << Form(file_path + "MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF.flat/0000000/MiniRun6.5_1E19_RHC.caf.%07d.CAF.flat.root", file_num) << std::endl;
            continue;
        }

        std::string tree_name = is_flat ? "rec" : "";
        auto sr = new caf::SRProxy(caf_tree, tree_name);

        // Open the companion GENIE systematics file
        TFile* genie_rw_file = TFile::Open(Form(file_path + "/systematics/nusystematics/MiniRun6.5.nusyst/MiniRun6.5_1E19_RHC.nuweights.%07d.nusyst.root", file_num), "READ");
        TTree* genie_rw_tree = (TTree*)genie_rw_file->Get("SystWeights");
        Double_t totWeight[100];
        genie_rw_tree->SetBranchAddress("totWeight", totWeight);

        // Loop over each spill
        const unsigned long nspills = caf_tree->GetEntries();
        const unsigned int incr = nspills / 10;
        std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;
        for (unsigned long i = 0; i < nspills; ++i)
        {
            caf_tree->GetEntry(i);

            // Keep track of spill # with print statement
            if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

            int spill_num = i;

            const auto num_ixn = sr->common.ixn.ndlp;

            // Loop over each reco interaction
            for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
            {
                // Reset data for each interaction
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

                // Genie systematic weights for this truth interaction
                const int genie_idx = truth_ixn.genieIdx;
                genie_rw_tree->GetEntry(genie_idx);
                data.genie_weights.assign(totWeight, totWeight + 100);
                data.genie_index = genie_idx;

                // If vertex is not contained or target is not argon, skip interaction
                if(contained(vtx.x, vtx.y, vtx.z) == false || truth_ixn.targetPDG != 1000180400)
                    continue;

                // Count number of relevant (reco) particles
                auto reco_nproton = 0;
                auto reco_npion = 0;
                auto reco_nmuon = 0;

                // Loop through particles in reco interaction
                for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
                {
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                    int pdg = part.pdg;

                    // Count protons, muons, and pions
                    if(pdg == 2212)
                        reco_nproton++;
                    if(pdg == 13 || pdg == -13)
                        reco_nmuon++;
                    if(pdg == 111 || pdg == 211 || pdg == -211)
                        reco_npion++;
                }

                reco_passes = reco_cuts(reco_nproton, reco_nmuon, reco_npion);

                // Count number of relevant (truth) particles
                auto truth_nproton = 0;
                auto truth_npion = 0;
                auto truth_nmuon = 0;
                for(unsigned long ipart = 0; ipart < truth_ixn.prim.size(); ++ipart)
                {
                    const auto& part = truth_ixn.prim[ipart];

                    if(part.pdg == 2212) 
                        truth_nproton++;

                    if(part.pdg == 13 || part.pdg == -13) // Muon (neutrino) and anti muon (anti neutrino)
                        truth_nmuon++;

                    if(part.pdg == 111 || part.pdg == 211 || part.pdg == -211)
                        truth_npion++;
                }

                truth_passes = truth_cuts(truth_nproton, truth_nmuon, truth_npion);

                // If interaction is not CC2p1mu0pi (reco or truth), go to next interaction      
                if(!reco_passes && !truth_passes)
                    continue;

                // Interaction passes, now count it as an event
                num_events++;

                // Loop over particles in reco interaction
                // Now to save information
                for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
                {
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];
                    int pdg = part.pdg;

                    // Save info for Minerva interaction track matching
                    int ixnM;
                    int idxM; 

                    // Get truth particle matches for reco particle
                    caf::Proxy<caf::SRTrueParticle>* truth_match = nullptr;
                    const auto& vec_truth_id = part.truth;
                    const auto& vec_overlap = part.truthOverlap;

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
                        truth_match = &(sr->mc.nu[truth_id.ixn].prim[truth_id.part]);
                    else if(truth_id.type == 3)
                        truth_match = &(sr->mc.nu[truth_id.ixn].sec[truth_id.part]);
                    else
                    {
                        std::cout << "Invalid truth id type!" << std::endl;
                        continue;
                    }

                    // Get Minerva match
                    bool minerva_track = false;
                    // Loop over primary tracks
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
                        double dirX = dX / length;
                        double dirY = dY / length;
                        double dirZ = dZ / length;

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
                                for(int k=0; k<sr->nd.trkmatch.extrap.size(); k++)
                                {
                                    if (sr->nd.trkmatch.extrap[k].larid.ixn!=ixn || sr->nd.trkmatch.extrap[k].larid.reco!=1) continue;
                                    int index=sr->nd.trkmatch.extrap[k].larid.idx;
                                    if (sr->nd.lar.dlp[ixn].tracks[index].start.z==sr->common.ixn.dlp[ixn].part.dlp[ipart].start.z && 
                                        sr->nd.lar.dlp[ixn].tracks[index].end.z==sr->common.ixn.dlp[ixn].part.dlp[ipart].end.z)
                                    {
                                        double dotProductTemp=abs(sr->nd.trkmatch.extrap[k].angdispl);
                                        if (dotProductTemp<0.99) continue;
                                        ixnM=sr->nd.trkmatch.extrap[k].minervaid.ixn;
                                        idxM=sr->nd.trkmatch.extrap[k].minervaid.idx; 
                                        // rec.nd.trkmatch.extrap.minervaid.ixn
                                        // rec.nd.trkmatch.extrap.minervaid.idx
                                        minerva_track = true;
                                    }
                                }
                            }
                        }
                    } // End of Minerva matching


                    // Get/calculate various reco/truth quantities
                    auto pvec = TVector3(part.p.x, part.p.y, part.p.z);
                    auto dir = TVector3(part.end.x, part.end.y, part.end.z) - 
                                TVector3(part.start.x, part.start.y, part.start.z);
                    //auto cos_angle = TMath::Cos(dir.Angle(beam_dir)); //Calculate cos of angle wrt neutrino beam direction
                    dir.RotateY(-TMath::Pi()/2);
                    //auto cos_rot_anode_angle = TMath::Cos(dir.Theta()); //Calculate cos of track rotational angle (projection on anode)
                    //auto cos_incl_anode_angle = TMath::Cos(dir.Phi()); //Calculate cos of track inclination angle (off of anode)
                    dir.RotateY(TMath::Pi()/2);
                    auto length = dir.Mag();

                    auto true_pvec = TVector3(truth_match->p.px, truth_match->p.py, truth_match->p.pz);
                    auto true_dir = TVector3(truth_match->end_pos.x, truth_match->end_pos.y, truth_match->end_pos.z) - 
                                    TVector3(truth_match->start_pos.x, truth_match->start_pos.y, truth_match->start_pos.z);
                    //auto true_cos_angle = TMath::Cos(true_dir.Angle(beam_dir));
                    true_dir.RotateY(TMath::Pi()/2);
                    auto true_length_val = true_dir.Mag();

                    // Currently unused
                    //auto T_diff = truth_match->p.E - part.E;
                    //auto p_diff = true_pvec.Mag() - pvec.Mag();
                    //auto length_diff = true_length_val - length;
                    //auto cos_angle_diff = true_cos_angle - cos_angle;

                    dir.RotateY(-TMath::Pi()/2);
                    true_dir.RotateY(-TMath::Pi()/2);

                    // Population information in vectors for tracks that have passed all cuts
                    // Reco
                    data.reco_energy.push_back(part.E);
                    data.reco_p_x.push_back(part.p.x);
                    data.reco_p_y.push_back(part.p.y);
                    data.reco_p_z.push_back(part.p.z);
                    data.reco_p_mag.push_back(pvec.Mag());
                    data.reco_length.push_back(length);
                    dir.RotateY(TMath::Pi()/2);
                    data.reco_angle.push_back(dir.Angle(beam_dir));
                    data.reco_angle_x.push_back(dir.Angle(x_plus_dir));
                    data.reco_angle_y.push_back(dir.Angle(y_plus_dir));
                    data.reco_angle_z.push_back(dir.Angle(z_plus_dir));
                    dir.RotateY(-TMath::Pi()/2);
                    data.reco_angle_rot.push_back(dir.Theta());
                    data.reco_angle_incl.push_back(dir.Phi());
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
                    data.true_energy.push_back(truth_match->p.E);
                    data.true_p_x.push_back(truth_match->p.px); 
                    data.true_p_y.push_back(truth_match->p.py); 
                    data.true_p_z.push_back(truth_match->p.pz);
                    data.true_p_mag.push_back(true_pvec.Mag());
                    data.true_length.push_back(true_length_val);
                    true_dir.RotateY(TMath::Pi()/2);
                    data.true_angle.push_back(true_dir.Angle(beam_dir));
                    data.true_angle_x.push_back(true_dir.Angle(x_plus_dir));
                    data.true_angle_y.push_back(true_dir.Angle(y_plus_dir));
                    data.true_angle_z.push_back(true_dir.Angle(z_plus_dir));
                    true_dir.RotateY(-TMath::Pi()/2);
                    data.true_angle_rot.push_back(true_dir.Theta());
                    data.true_angle_incl.push_back(true_dir.Phi());
                    data.true_vtx_x = (sr->mc.nu[truth_id.ixn].vtx.x);
                    data.true_vtx_y = (sr->mc.nu[truth_id.ixn].vtx.y);
                    data.true_vtx_z = (sr->mc.nu[truth_id.ixn].vtx.z);
                    data.true_track_start_x.push_back(truth_match->start_pos.x);
                    data.true_track_start_y.push_back(truth_match->start_pos.y);
                    data.true_track_start_z.push_back(truth_match->start_pos.z);
                    data.true_track_end_x.push_back(truth_match->end_pos.x);
                    data.true_track_end_y.push_back(truth_match->end_pos.y);
                    data.true_track_end_z.push_back(truth_match->end_pos.z);
                    data.true_pdg = (sr->mc.nu[truth_id.ixn].pdg);
                    data.true_nproton = (truth_nproton); //rec.mc.nu.nproton
                    data.true_nmuon = (truth_nmuon); //rec.mc.nu.nmuon
                    data.true_npion = (truth_npion); //rec.mc.nu.npion
                    data.mode = (sr->mc.nu[truth_id.ixn].mode); //rec.mc.nu.mode
                    data.nu_momentum_x = (sr->mc.nu[truth_id.ixn].momentum.x);
                    data.nu_momentum_y = (sr->mc.nu[truth_id.ixn].momentum.y);
                    data.nu_momentum_z = (sr->mc.nu[truth_id.ixn].momentum.z);

                    // Minerva
                    if (minerva_track == true)
                    {
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

            } // End of interaction loop
        
        } // End of spill loop
        caf_file->Close();
        genie_rw_file->Close();
        delete sr;
    } // End of file loop

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    // POPULATE: Write to output ROOT file
    fCafTree->Write();
        
    std::cout << "Wrote TTree." << std::endl;

    // CLOSE: Output ROOT file
    f->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}

int main()
{
    caf_plotter(true);

    return 0;
}