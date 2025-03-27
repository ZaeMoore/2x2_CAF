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
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_05_00/include/duneanaobj/StandardRecord/StandardRecord.h"
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_05_00/include/duneanaobj/StandardRecord/Proxy/SRProxy.h"
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
    std::vector< double >  reco_ixn_index;

    std::vector< double >  reco_pandora_energy;
    std::vector< double >  reco_pandora_p_x;
    std::vector< double >  reco_pandora_p_y;
    std::vector< double >  reco_pandora_p_z;
    std::vector< double >  reco_pandora_p_mag;
    std::vector< double >  reco_pandora_length;
    std::vector< double >  reco_pandora_angle;

    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    std::vector< int >     event;
    std::vector< int >     run;
    std::vector< int >     subrun;
    std::vector<std::string> caf_file_name;


    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fCafTree=new TTree("CafTree", "CAF Variables");
    fCafTree->Branch("reco_energy", &reco_energy);
    fCafTree->Branch("reco_p_x", &reco_p_x);
    fCafTree->Branch("reco_p_y", &reco_p_y);
    fCafTree->Branch("reco_p_z", &reco_p_z);
    fCafTree->Branch("reco_p_mag", &reco_p_mag);
    fCafTree->Branch("reco_length", &reco_length);
    fCafTree->Branch("reco_angle", &reco_angle);
    fCafTree->Branch("reco_angle_rot", &reco_angle_rot);
    fCafTree->Branch("reco_angle_incl", &reco_angle_incl);
    fCafTree->Branch("reco_angle_x", &reco_angle_x);
    fCafTree->Branch("reco_angle_y", &reco_angle_y);
    fCafTree->Branch("reco_angle_z", &reco_angle_z);
    fCafTree->Branch("reco_track_start_x", &reco_track_start_x);
    fCafTree->Branch("reco_track_start_y", &reco_track_start_y);
    fCafTree->Branch("reco_track_start_z", &reco_track_start_z);
    fCafTree->Branch("reco_track_end_x", &reco_track_end_x);
    fCafTree->Branch("reco_track_end_y", &reco_track_end_y);
    fCafTree->Branch("reco_track_end_z", &reco_track_end_z);
    fCafTree->Branch("reco_pdg", &reco_pdg);
    fCafTree->Branch("reco_ixn_index", &reco_ixn_index);

    fCafTree->Branch("spill_index", &spill_index);
    fCafTree->Branch("file_index", &file_index);
    fCafTree->Branch("event", &event);
    fCafTree->Branch("run", &run);
    fCafTree->Branch("subrun", &subrun);
    fCafTree->Branch("caf_file_name", &caf_file_name);


    // Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    // z-direction (roughly beam dir)
    const auto z_plus_dir = TVector3(0, 0, 1.0);
    const auto y_plus_dir = TVector3(0, 1.0, 0);
    const auto x_plus_dir = TVector3(1.0, 0, 0);

    // negative y-direction
    const auto y_minus_dir = TVector3(0, -1.0, 0);

    // Loop through files in list
    const auto t_start{std::chrono::steady_clock::now()};
    auto file_num = 0;
    for(const auto& f : root_list)
    { // Beginning of file loop
        std::string current_file = f;
        std::cout << "Processing " << f << std::endl;
        file_num++;

        // Open file and attach SRProxy Object
        TFile* caf_file = TFile::Open(f.c_str(), "READ");
        TTree* caf_tree = (TTree*)caf_file->Get("cafTree");
        std::string tree_name = is_flat ? "rec" : "";
        auto sr = new caf::SRProxy(caf_tree, tree_name);

        // Loop over each spill
        const unsigned long nspills = caf_tree->GetEntries();
        const unsigned int incr = nspills / 10;
        std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;
        for (unsigned long i = 0; i < nspills; ++i)
        { // Beginning of spill loop
            caf_tree->GetEntry(i);

            // Keep track of spill # with print statement
            if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

            int spill_num = i;

            const auto num_ixn = sr->common.ixn.ndlp; // spine
            //const auto num_ixn = sr->common.ixn.npandora; // pandora

            // Loop over each reco interaction
            for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
            { // Beginning of interaction loop
                const auto& reco_ixn = sr->common.ixn.dlp[ixn];

                // Require vertex to be within volume -------------------------------------------
                bool is_contained = contained(reco_ixn.vtx.x, reco_ixn.vtx.y, reco_ixn.vtx.z);
                if(is_contained == false)
                    continue;

                auto nproton = 0;
                auto npion = 0;
                auto nmuon = 0;

                // Require 2 protons and a muon to start from the same vertex, and no more than 1 pion from that vertex
                bool vtx_2pmu = false;
                std::vector<double> mu_start_x;
                std::vector<double> mu_start_y;
                std::vector<double> mu_start_z;
                std::vector<double> p_start_x;
                std::vector<double> p_start_y;
                std::vector<double> p_start_z;
                std::vector<double> pi_start_x;
                std::vector<double> pi_start_y;
                std::vector<double> pi_start_z;

                // Loop over particles in interaction, count the number of protons, pions, and muons, and cut non-viable interactions
                // Require that at least 2 protons and 1 muon originate from the same vertex
                for(unsigned long ipart = 0; ipart < reco_ixn.part.dlp.size(); ++ipart)
                { // Beginning of particle loop 1
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                    if(part.pdg == 2212) // proton
                    {
                        nproton++;
                        p_start_x.push_back(part.start.x);
                        p_start_y.push_back(part.start.y);
                        p_start_z.push_back(part.start.z);
                    }
                        

                    if(part.pdg == 13) // muon
                    {
                        nmuon++;
                        mu_start_x.push_back(part.start.x);
                        mu_start_y.push_back(part.start.y);
                        mu_start_z.push_back(part.start.z);
                    }
       
                    if(part.pdg == 111 || part.pdg == 211 || part.pdg == -211) // all pions
                    {
                        npion++;
                        pi_start_x.push_back(part.start.x);
                        pi_start_y.push_back(part.start.y);
                        pi_start_z.push_back(part.start.z);
                    }

                } // End of particle loop 1

                // Loop through muon start positions
                for(unsigned long i = 0; i < mu_start_x.size(); ++i)
                { // Beginning of muon start pos loop
                    auto vtx_match = 0; // Add 1 for each proton that starts in the same pos as muon
                    auto mu_start = TVector3(mu_start_x[i], mu_start_y[i], mu_start_z[i]);

                    for(unsigned long j = 0; j < p_start_x.size(); ++j)
                    { // Beginning of proton start pos loop
                        auto p_start = TVector3(p_start_x[j], p_start_y[j], p_start_z[j]);

                        // Compare distance, must be less than something
                        TVector3 diff = mu_start - p_start;
                        auto magnitude = diff.Mag();

                        // if distance less than something, ++vtx_match
                        if(magnitude < 5)
                            vtx_match++;
                    } // End of proton start pos loop
                    
                    if(vtx_match >= 2)
                        vtx_2pmu = true;
                } // End of muon start pos loop


                if(nproton < 2 || nmuon < 1 || npion > 0 || vtx_2pmu == false)
                    continue;

                // Loop over particles in interaction and save information
                for(unsigned long ipart = 0; ipart < reco_ixn.part.dlp.size(); ++ipart)
                { // Beginning of particle loop 2
                    const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                    // Get/calculate various reco quantities and save to output root file
                    auto pvec = TVector3(part.p.x, part.p.y, part.p.z);
                    auto dir = TVector3(part.end.x, part.end.y, part.end.z) - TVector3(part.start.x, part.start.y, part.start.z);
                    //auto cos_angle = TMath::Cos(dir.Angle(beam_dir)); // Cos of angle wrt neutrino beam direction
                    dir.RotateY(-TMath::Pi()/2);
                    //auto cos_rot_anode_angle = TMath::Cos(dir.Theta()); // Cos of track rotational angle (projection on anode)
                    //auto cos_incl_anode_angle = TMath::Cos(dir.Phi()); // Cos of track inclination angle (off of anode)
                    dir.RotateY(TMath::Pi()/2);
                    auto length = dir.Mag();

                    // Populate information in vectors for tracks that have passed all cuts
                    reco_energy.push_back(part.E);
                    reco_p_x.push_back(part.p.x);
                    reco_p_y.push_back(part.p.y);
                    reco_p_z.push_back(part.p.z);
                    reco_p_mag.push_back(pvec.Mag());
                    reco_length.push_back(length);
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
                    reco_ixn_index.push_back(ixn);
                    spill_index.push_back(spill_num);
                    file_index.push_back(file_num);
                    event.push_back(sr->meta.nd_lar.event);
                    subrun.push_back(sr->meta.nd_lar.subrun);
                    caf_file_name.push_back(current_file.erase(0, current_file.find_last_of("/")+1).c_str());

                } // End of particle loop 2

            } // End of interaction loop

        } // End of spill loop
        caf_file->Close();
    } // End of file loop

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    // Output TTree file name
    std::string file_name = "2x2_2p2h_output_1.3";

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