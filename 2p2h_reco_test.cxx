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

    //Create TChain and add files to it
    //TChain* caf_chain = new TChain("cafTree");
    //for(const auto& file : root_list)
    //{
    //    std::cout << "Adding " << file << " to TChain." << std::endl;
    //    caf_chain->Add(file.c_str());
    //}

    std::cout << "Finished adding files..." << std::endl;

    // DEFINE: Vectors to hold information to keep in output TTree file 
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
    std::vector< int >     true_ixn_charged_track_mult;
    std::vector< int >     true_nproton;
    std::vector< int >     true_nmuon;
    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    std::vector< int >     interaction_id;
    std::vector< int >     event;
    std::vector< int >     run;
    std::vector< int >     subrun;
    std::vector<std::string> caf_file_name;


    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fCafTree=new TTree("CafTree", "Truth Variables");
    fCafTree->Branch("true_energy", &true_energy);
    fCafTree->Branch("true_p_x", &true_p_x);
    fCafTree->Branch("true_p_y", &true_p_y);
    fCafTree->Branch("true_p_z", &true_p_z);
    fCafTree->Branch("true_p_mag", &true_p_mag);
    fCafTree->Branch("true_length", &true_length);
    fCafTree->Branch("true_angle", &true_angle);
    fCafTree->Branch("true_angle_rot", &true_angle_rot);
    fCafTree->Branch("true_angle_incl", &true_angle_incl);
    fCafTree->Branch("true_angle_x", &true_angle_x);
    fCafTree->Branch("true_angle_y", &true_angle_y);
    fCafTree->Branch("true_angle_z", &true_angle_z);
    fCafTree->Branch("true_track_start_x", &true_track_start_x);
    fCafTree->Branch("true_track_start_y", &true_track_start_y);
    fCafTree->Branch("true_track_start_z", &true_track_start_z);
    fCafTree->Branch("true_track_end_x", &true_track_end_x);
    fCafTree->Branch("true_track_end_y", &true_track_end_y);
    fCafTree->Branch("true_track_end_z", &true_track_end_z);
    fCafTree->Branch("true_pdg", &true_pdg);
    fCafTree->Branch("true_ixn_charged_track_mult", &true_ixn_charged_track_mult);
    fCafTree->Branch("true_nproton", &true_nproton);
    fCafTree->Branch("true_nmuon", &true_nmuon);
    fCafTree->Branch("spill_index", &spill_index);
    fCafTree->Branch("file_index", &file_index);
    fCafTree->Branch("interaction_id", &interaction_id);
    fCafTree->Branch("event", &event);
    fCafTree->Branch("run", &run);
    fCafTree->Branch("subrun", &subrun);
    fCafTree->Branch("caf_file_name", &caf_file_name);


    //Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    //z-direction (roughly beam dir)
    const auto z_plus_dir = TVector3(0, 0, 1.0);
    const auto y_plus_dir = TVector3(0, 1.0, 0.0);
    const auto x_plus_dir = TVector3(1.0, 0, 0.0);

    //negative y-direction 
    const auto y_minus_dir = TVector3(0, -1.0, 0.0);

    //Center of the 2x2 LAr
    //const float tpc_x = 0.0;
    //const float tpc_y = 0.0; // -268.0;
    //const float tpc_z = 0.0; //(1333.5 + 1266.5) / 2.0;

    const auto t_start{std::chrono::steady_clock::now()};
    auto file_num = 0;
    for(const auto& f : root_list)
    {
        std::cout << "Processing " << f << std::endl;
        file_num++;

        //Open file and attach SRProxy object
        //Different tree name for structured and flat CAFs
        TFile* caf_file = TFile::Open(f.c_str(), "READ");
        TTree* caf_tree = (TTree*)caf_file->Get("cafTree");
        std::string tree_name = is_flat ? "rec" : "";
        auto sr = new caf::SRProxy(caf_tree, tree_name);
        
        //First loop over each spill, then each reco interaction, then each reco particle
        const unsigned long nspills = caf_tree->GetEntries();
        const unsigned int incr = nspills / 10;
        std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;

        // Loop over each spill
        for(unsigned long i = 0; i < nspills; ++i)
        {
            caf_tree->GetEntry(i);

            if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

            int spill_num = i;

            const auto num_ixn = sr->mc.nu.size();

            // Loop over each interaction
            for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
            {
                const auto& truth_ixn = sr->mc.nu[ixn];

                // Put cuts on true interaction quantities
                // Reject interactions not on argon
                if(truth_ixn.targetPDG != 1000180400)
                    continue;

                bool is_contained = true;
                unsigned int npi = truth_ixn.npi0 + truth_ixn.npim + truth_ixn.npip;

                //if(truth_ixn.nproton >= 2 and npi == 0 and truth_ixn.mode == 10) 
                // Reject interactions with less than 2 protons or more than 0 pions
                if(truth_ixn.nproton < 2 or npi != 0)
                    continue;
                
                int nproton = truth_ixn.nproton;
                int nmuon = 0;

                // Loop over each particle in interaction and cut events with low energy particles
                for(unsigned long ipart = 0; ipart < sr->mc.nu[ixn].prim.size(); ++ipart)
                {
                    auto pmag = (TVector3(sr->mc.nu[ixn].prim[ipart].p.px, sr->mc.nu[ixn].prim[ipart].p.py, sr->mc.nu[ixn].prim[ipart].p.pz)).Mag();

                    if(sr->mc.nu[ixn].prim[ipart].pdg == 2212 and pmag <= 0.05) // Cut low momentum protons < 50 MeV
                        nproton--;
                        
                    if(sr->mc.nu[ixn].prim[ipart].pdg == 13)
                        nmuon++;
                }   

                is_contained = contained(sr->mc.nu[ixn].vtx.x, sr->mc.nu[ixn].vtx.y, sr->mc.nu[ixn].vtx.z);

                // If # protons above 5 MeV are <2, the interaction isn't contained, or there are 0 muons, skip this event
                if(nproton < 2 or is_contained == false or nmuon == 0)
                    continue;

                bool goodInteraction=false;   
                std::vector<int> trueInteractionIndex;
                // ----------------------Find matching reco interaction----------------------
                // Loop through reco interactions
                for(long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++){

                    bool goodInteraction=false;
                    double biggestMatch=-999;
                    int biggestMatchIndex=-999;
                    int maxEventPar=-999; 
                    int maxEventTyp=-9999; 
                    int maxEventIxn=-999;

                    // Loop through all the truth match values for each interaction
                    for (int ntruth=0; ntruth<sr->common.ixn.dlp[nixn].truth.size(); ntruth++){
                        if (biggestMatch<sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)){
                        biggestMatch=sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
                        biggestMatchIndex=sr->common.ixn.dlp[nixn].truth.at(ntruth);
                        }
                    }
                }

                // Loop over each particle in interaction
                for(unsigned long ipart = 0; ipart < sr->mc.nu[ixn].prim.size(); ++ipart)
                {
                    // Store current true particle for easier access
                    const auto& true_part = sr->mc.nu[ixn].prim[ipart];

                    // Make cuts for individual particles if needed
                    // Cut phantom particles
                    if(true_part.pdg == 0 or std::isnan(sr->mc.nu[ixn].prim[ipart].start_pos.x))
                        continue;

                    // Cut protons below 5 MeV
                    auto pmag = (TVector3(sr->mc.nu[ixn].prim[ipart].p.px, sr->mc.nu[ixn].prim[ipart].p.py, sr->mc.nu[ixn].prim[ipart].p.pz)).Mag();
                    if(sr->mc.nu[ixn].prim[ipart].pdg == 2212 and pmag <= 0.005)
                        continue;

                    // Finally get or calculate various truth quantities
                    auto true_pvec = TVector3(true_part.p.px, true_part.p.py, true_part.p.pz);
                    auto true_dir = TVector3(true_part.end_pos.x, true_part.end_pos.y, true_part.end_pos.z)
                                    - TVector3(true_part.start_pos.x, true_part.start_pos.y, true_part.start_pos.z);
                    //auto true_cos_angle = TMath::Cos(true_dir.Angle(beam_dir));
                    true_dir.RotateY(-TMath::Pi()/2);
                    //auto true_cos_rot_anode_angle = TMath::Cos(true_dir.Theta()); //calculate cosine of track rotational angle (projection on anode)
                    //auto true_cos_incl_anode_angle = TMath::Cos(true_dir.Phi()); //calculate cosine of track inclination angle (off of anode)
                    true_dir.RotateY(TMath::Pi()/2);
                    auto true_length_val = true_dir.Mag();

                    true_dir.RotateY(-TMath::Pi()/2);

                    // POPULATE: Record information in vectors (defined above) for tracks that have passed all cuts
                    true_energy.push_back(true_part.p.E);
                    true_p_x.push_back(true_part.p.px); 
                    true_p_y.push_back(true_part.p.py); 
                    true_p_z.push_back(true_part.p.pz);
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
                    true_track_start_x.push_back(true_part.start_pos.x);
                    true_track_start_y.push_back(true_part.start_pos.y);
                    true_track_start_z.push_back(true_part.start_pos.z);
                    true_track_end_x.push_back(true_part.end_pos.x);
                    true_track_end_y.push_back(true_part.end_pos.y);
                    true_track_end_z.push_back(true_part.end_pos.z);
                    true_pdg.push_back(true_part.pdg);
                    // true_ixn_charged_track_mult.push_back(true_ixn_trks); Do the stuff to make this line work
                    true_nproton.push_back(nproton);
                    true_nmuon.push_back(nmuon);
                    spill_index.push_back(spill_num);
                    file_index.push_back(file_num);
                    interaction_id.push_back(true_part.interaction_id);
                    event.push_back(sr->meta.nd_lar.event);
                    run.push_back(sr->meta.nd_lar.run);
                    subrun.push_back(sr->meta.nd_lar.subrun);
                    caf_file_name.push_back(f.c_str());

                }

            }

        }
        caf_file->Close();
    }

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    // Output TTree file name
    std::string file_name = "2p2h_truth_minirun6";

    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");

    // POPULATE: Fill TTree and write to output ROOT file
    fCafTree->Fill();
    fCafTree->Write();
    
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