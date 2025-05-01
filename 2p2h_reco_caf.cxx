#include <iostream>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TVector.h>
#include <TVector3.h>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_01_00/include/duneanaobj/StandardRecord/StandardRecord.h"
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_01_00/include/duneanaobj/StandardRecord/Proxy/SRProxy.h"
#define Dimension 3

float dot_product(std::vector<float> vector_a, std::vector<float> vector_b) {
  float product = 0;
  for (int i = 0; i < Dimension; i++)
    product = product + vector_a[i] * vector_b[i];
  return product;
}

bool contained(double x, double y, double z){
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

int caf_plotter(std::string file_list, bool is_flat = true)
{

    std::vector<std::string> root_list;
    std::ifstream fin(file_list, std::ios::in);

    //Check if file exists
    if(!fin.is_open())
    {
        std::cerr << "Failed to open " << file_list << std::endl;
        std::cerr << "Exiting" << std::endl;
        return 111;
    }
    else
    {
        std::cout << "Reading " << file_list << " for input ROOT files." << std::endl;
        std::string name;
        
        //Add name of each file to root_list
        while(std::getline(fin, name))
        {
            root_list.push_back(name);
        }
    }

    //Check if list of files is empty
    if(root_list.empty())
    {
        std::cerr << "No input ROOT files. Exiting." << std::endl;
        return 121;
    }

    std::cout << "Finished adding files..." << std::endl;

    // DEFINE: Vectors to hold information to keep in output TTree file
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
    std::vector< int >     true_nproton;
    std::vector< int >     true_nmuon;
    std::vector< int >     interaction_id;

    std::vector< double >  overlap;
    std::vector< double >  true_ixn_index;
    std::vector< double >  reco_ixn_index;
    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    std::vector< int >     event;
    std::vector< int >     run;
    std::vector< int >     subrun;
    std::vector<std::string> caf_file_name;


    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fRecoTree=new TTree("RecoTree", "Caf reco and truth variables");
    fRecoTree->Branch("reco_energy", &reco_energy);
    fRecoTree->Branch("reco_p_x", &reco_p_x);
    fRecoTree->Branch("reco_p_y", &reco_p_y);
    fRecoTree->Branch("reco_p_z", &reco_p_z);
    fRecoTree->Branch("reco_p_mag", &reco_p_mag);
    fRecoTree->Branch("reco_length", &reco_length);
    fRecoTree->Branch("reco_angle", &reco_angle);
    fRecoTree->Branch("reco_angle_rot", &reco_angle_rot);
    fRecoTree->Branch("reco_angle_incl", &reco_angle_incl);
    fRecoTree->Branch("reco_angle_x", &reco_angle_x);
    fRecoTree->Branch("reco_angle_y", &reco_angle_y);
    fRecoTree->Branch("reco_angle_z", &reco_angle_z);
    fRecoTree->Branch("reco_track_start_x", &reco_track_start_x);
    fRecoTree->Branch("reco_track_start_y", &reco_track_start_y);
    fRecoTree->Branch("reco_track_start_z", &reco_track_start_z);
    fRecoTree->Branch("reco_track_end_x", &reco_track_end_x);
    fRecoTree->Branch("reco_track_end_y", &reco_track_end_y);
    fRecoTree->Branch("reco_track_end_z", &reco_track_end_z);
    fRecoTree->Branch("reco_pdg", &reco_pdg);
    fRecoTree->Branch("reco_ixn_index", &reco_ixn_index);

    fRecoTree->Branch("true_energy", &true_energy);
    fRecoTree->Branch("true_p_x", &true_p_x);
    fRecoTree->Branch("true_p_y", &true_p_y);
    fRecoTree->Branch("true_p_z", &true_p_z);
    fRecoTree->Branch("true_p_mag", &true_p_mag);
    fRecoTree->Branch("true_length", &true_length);
    fRecoTree->Branch("true_angle", &true_angle);
    fRecoTree->Branch("true_angle_rot", &true_angle_rot);
    fRecoTree->Branch("true_angle_incl", &true_angle_incl);
    fRecoTree->Branch("true_angle_x", &true_angle_x);
    fRecoTree->Branch("true_angle_y", &true_angle_y);
    fRecoTree->Branch("true_angle_z", &true_angle_z);
    fRecoTree->Branch("true_track_start_x", &true_track_start_x);
    fRecoTree->Branch("true_track_start_y", &true_track_start_y);
    fRecoTree->Branch("true_track_start_z", &true_track_start_z);
    fRecoTree->Branch("true_track_end_x", &true_track_end_x);
    fRecoTree->Branch("true_track_end_y", &true_track_end_y);
    fRecoTree->Branch("true_track_end_z", &true_track_end_z);
    fRecoTree->Branch("true_pdg", &true_pdg);
    fRecoTree->Branch("true_nproton", &true_nproton);
    fRecoTree->Branch("true_nmuon", &true_nmuon);
    fRecoTree->Branch("interaction_id", &interaction_id);
    fRecoTree->Branch("true_ixn_index", &true_ixn_index);

    fRecoTree->Branch("overlap", &overlap);
    fRecoTree->Branch("spill_index", &spill_index);
    fRecoTree->Branch("file_index", &file_index);
    fRecoTree->Branch("event", &event);
    fRecoTree->Branch("run", &run);
    fRecoTree->Branch("subrun", &subrun);
    fRecoTree->Branch("caf_file_name", &caf_file_name);


    //Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    //z-direction (roughly beam dir)
    const auto z_plus_dir = TVector3(0, 0, 1.0);
    const auto y_plus_dir = TVector3(0, 1.0, 0.0);
    const auto x_plus_dir = TVector3(1.0, 0, 0.0);

    //negative y-direction 
    const auto y_minus_dir = TVector3(0, -1.0, 0.0);

    //Loop through files in list
    const auto t_start{std::chrono::steady_clock::now()};
    auto file_num = 0;
    for(const auto& f : root_list)
    {
        std::string current_file = f;
        std::cout << "Processing " << f << std::endl;
        file_num++;

        //Open file and attach SRProxy Object
        TFile* caf_file = TFile::Open(f.c_str(), "READ");
        TTree* caf_tree = (TTree*)caf_file->Get("cafTree");
        std::string tree_name = is_flat ? "rec" : "";
        auto sr = new caf::SRProxy(caf_tree, tree_name);

        //Loop over each spill
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
                const auto& vtx = sr->common.ixn.dlp[ixn].vtx;

                // Get the truth interaction(s) corresponding to this reco interaction
                const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
                const auto& vec_overlap_ixn = sr->common.ixn.dlp[ixn].truthOverlap;

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

                // Require vertex to be within 
                bool is_contained = true;
                is_contained = contained(truth_ixn.vtx.x, truth_ixn.vtx.y, truth_ixn.vtx.z);

                // truth_ixn.ix > 1E9
                if(is_contained == false || truth_ixn.targetPDG != 1000180400)
                    continue;


                // Count number of relevant particles
                // Don't trust nproton branch 
                auto truth_nproton = 0;
                auto truth_npion = 0;
                auto truth_nmuon = 0;

                // Loop over true particles and count number of muons and protons
                for(unsigned long ipart = 0; ipart < truth_ixn.prim.size(); ++ipart)
                {
                    const auto& part = truth_ixn.prim[ipart];

                    auto pmag = (TVector3(part.p.px, part.p.py, part.p.pz)).Mag();

                    if(part.pdg == 2212) // Cut low momentum protons under 10 MeV
                        truth_nproton++;

                    if(part.pdg == 13 || part.pdg == -13)
                        truth_nmuon++;

                    if(part.pdg == 111 || part.pdg == 211 || part.pdg == -211)
                        truth_npion++;
                }

                // If interaction is not CC2p1mu0pi, go to next interaction             
                if(truth_nproton != 2 || truth_nmuon != 1 || truth_npion > 0)
                    continue;

                // Loop over particles in reco interaction
                for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
                {
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

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

                    // Get/calculate various reco/truth quantities
                    auto pvec = TVector3(part.p.x, part.p.y, part.p.z);
                    auto dir = TVector3(part.end.x, part.end.y, part.end.z) - TVector3(part.start.x, part.start.y, part.start.z);
                    auto cos_angle = TMath::Cos(dir.Angle(beam_dir)); //Calculate cos of angle wrt neutrino beam direction
                    dir.RotateY(-TMath::Pi()/2);
                    auto cos_rot_anode_angle = TMath::Cos(dir.Theta()); //Calculate cos of track rotational angle (projection on anode)
                    auto cos_incl_anode_angle = TMath::Cos(dir.Phi()); //Calculate cos of track inclination angle (off of anode)
                    dir.RotateY(TMath::Pi()/2);
                    auto length = dir.Mag();

                    auto true_pvec = TVector3(truth_match->p.px, truth_match->p.py, truth_match->p.pz);
                    auto true_dir = TVector3(truth_match->end_pos.x, truth_match->end_pos.y, truth_match->end_pos.z) - TVector3(truth_match->start_pos.x, truth_match->start_pos.y, truth_match->start_pos.z);
                    auto true_cos_angle = TMath::Cos(true_dir.Angle(beam_dir));
                    true_dir.RotateY(TMath::Pi()/2);
                    auto true_length_val = true_dir.Mag();

                    auto T_diff = truth_match->p.E - part.E;
                    auto p_diff = true_pvec.Mag() - pvec.Mag();
                    auto length_diff = true_length_val - length;
                    auto cos_angle_diff = true_cos_angle - cos_angle;

                    dir.RotateY(-TMath::Pi()/2);
                    true_dir.RotateY(-TMath::Pi()/2);

                    //Population information in vectors for tracks that have pass all cuts
                    reco_energy.push_back(part.E);
                    reco_p_x.push_back(part.p.x);
                    reco_p_y.push_back(part.p.y);
                    reco_p_z.push_back(part.p.z);
                    reco_p_mag.push_back(pvec.Mag());
                    reco_length.push_back(length);
                    dir.RotateY(TMath::Pi()/2);
                    reco_angle.push_back(dir.Angle(beam_dir));
                    reco_angle_x.push_back(dir.Angle(x_plus_dir));
                    reco_angle_y.push_back(dir.Angle(y_plus_dir));
                    reco_angle_z.push_back(dir.Angle(z_plus_dir));
                    dir.RotateY(-TMath::Pi()/2);
                    reco_angle_rot.push_back(dir.Theta());
                    reco_angle_incl.push_back(dir.Phi());
                    reco_track_start_x.push_back(part.start.x);
                    reco_track_start_y.push_back(part.start.y);
                    reco_track_start_z.push_back(part.start.z);
                    reco_track_end_x.push_back(part.end.x);
                    reco_track_end_y.push_back(part.end.y);
                    reco_track_end_z.push_back(part.end.z);
                    reco_pdg.push_back(part.pdg);
                    true_energy.push_back(truth_match->p.E);
                    true_p_x.push_back(truth_match->p.px); 
                    true_p_y.push_back(truth_match->p.py); 
                    true_p_z.push_back(truth_match->p.pz);
                    true_p_mag.push_back(true_pvec.Mag());
                    true_length.push_back(true_length_val);
                    true_dir.RotateY(TMath::Pi()/2);
                    true_angle.push_back(true_dir.Angle(beam_dir));
                    true_angle_x.push_back(true_dir.Angle(x_plus_dir));
                    true_angle_y.push_back(true_dir.Angle(y_plus_dir));
                    true_angle_z.push_back(true_dir.Angle(z_plus_dir));
                    true_dir.RotateY(-TMath::Pi()/2);
                    true_angle_rot.push_back(true_dir.Theta());
                    true_angle_incl.push_back(true_dir.Phi());
                    true_track_start_x.push_back(truth_match->start_pos.x);
                    true_track_start_y.push_back(truth_match->start_pos.y);
                    true_track_start_z.push_back(truth_match->start_pos.z);
                    true_track_end_x.push_back(truth_match->end_pos.x);
                    true_track_end_y.push_back(truth_match->end_pos.y);
                    true_track_end_z.push_back(truth_match->end_pos.z);
                    true_pdg.push_back(truth_match->pdg);
                    true_nproton.push_back(sr->mc.nu[truth_id.ixn].nproton); //rec.mc.nu.nproton
                    overlap.push_back(current_max);
                    true_ixn_index.push_back(truth_idx);
                    reco_ixn_index.push_back(ixn);
                    spill_index.push_back(spill_num);
                    file_index.push_back(file_num);
                    event.push_back(sr->meta.nd_lar.event);
                    run.push_back(sr->meta.nd_lar.run);
                    subrun.push_back(sr->meta.nd_lar.subrun);
                    caf_file_name.push_back(current_file.erase(0, current_file.find_last_of("/")+1).c_str());

                } // End of particle loop

            } // End of interaction loop
        
        } // End of spill loop
        caf_file->Close();
    } // End of file loop

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    // Output TTree file name
    std::string file_name = "2p2h_efficiency_output";

    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");

    // POPULATE: Fill TTree and write to output ROOT file
    fRecoTree->Fill();
    fRecoTree->Write();
        
    std::cout << "Filled and wrote TTree." << std::endl;

    // CLOSE: Output ROOT file
    f->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}

int main(int argc, char** argv)
{
    std::string input_file_list = argv[1];

    caf_plotter(input_file_list, true);

    return 0;
}