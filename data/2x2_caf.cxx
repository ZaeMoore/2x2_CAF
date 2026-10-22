/*
2x2 data analysis script

Analysis of 2x2 data (sandbox) files
Selects only multi-proton zero-pion final state events and saves them to a separate root file

To run:
Compile with compile.sh (edit the script to point to caf_plotter.cxx file)
Run with ./2x2_caf arg1

Inputs:
arg1 = path to text file containing list of input root files
soon: arg2 = path to text file containing list of input systematics root files

Output:
Root file with TTree containing only multi-proton zero-pion final state events
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
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#define Dimension 3

float dot_product(std::vector<float> vector_a, std::vector<float> vector_b) {
  float product = 0;
  for (int i = 0; i < Dimension; i++)
    product = product + vector_a[i] * vector_b[i];
  return product;
}

bool contained(double x, double y, double z){
    double tpc_dist = 5.0; // distance from the tpc walls for containment cuts
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


struct PartInfo {
    int pdg;
    TVector3 start;
    TVector3 end;
};

/*
Apply full cuts
- Counts 2+ protons, exactly 1 muon, and 0 pions
- Vertex is contained in TPC
- Requires all protons to be contained in the TPC
- Requires at least 2 protons to start within 3 cm of the muon start
- Requires muon to punch out the back of MINERvA
- Fewer than max_other other particles (not p, mu, or pi)
*/
bool cuts(const std::vector<PartInfo>& parts,
          const TVector3& vtx,
          double minerva_track_end_z,
          double minerva_back_z = 62.0,
          int max_other = 5)
{
    // Vertex containment
    if(!contained(vtx.X(), vtx.Y(), vtx.Z()))
        return false;

    int nprotons = 0, nmuons = 0, npions = 0, nother = 0;
    std::vector<TVector3> p_starts;
    TVector3 mu_start;

    for(const auto& p : parts)
    {
        const int apdg = std::abs(p.pdg);

        if(apdg == 2212)
        {
            // Proton must be fully contained in the TPC
            if(!contained(p.end.X(),   p.end.Y(),   p.end.Z()))
                return false;
            nprotons++;
            p_starts.push_back(p.start);
        }
        else if(apdg == 13)
        {
            nmuons++;
            mu_start = p.start;
        }
        else if(apdg == 211 || p.pdg == 111)
            npions++;
        else
            nother++;
    }

    // Topology: CC1mu, 2+ protons, 0 pions
    if(nprotons < 2 || nmuons != 1 || npions > 0)
        return false;

    // Fewer than max_other other particles
    if(nother >= max_other)
        return false;

    // At least 2 protons share a vertex with the muon (within 3 cm)
    int vtx_match = 0;
    for(const auto& ps : p_starts)
        if((mu_start - ps).Mag() < 3.0)
            vtx_match++;
    if(vtx_match < 2)
        return false;

    // Muon exits the back of MINERvA
    if(minerva_track_end_z < minerva_back_z)
        return false;

    return true;
}

/*
Define interaction data
This will be called at the beginning of the caf plotter function
to define these variables for the output tree
Then will be re-called at the beginning of every interaction loop
*/
struct InteractionData
{
    std::vector< double >  track_dlp_e;
    std::vector< double >  track_dlp_px; 
    std::vector< double >  track_dlp_py; 
    std::vector< double >  track_dlp_pz;
    std::vector< double >  track_dlp_pmag;
    std::vector< double >  track_dlp_length;
    std::vector< double >  track_dlp_angle;
    std::vector< double >  track_dlp_angle_rot;
    std::vector< double >  track_dlp_angle_incl;
    std::vector< double >  track_dlp_angle_x;
    std::vector< double >  track_dlp_angle_y;
    std::vector< double >  track_dlp_angle_z;
    std::vector< double >  track_dlp_start_x;
    std::vector< double >  track_dlp_start_y;
    std::vector< double >  track_dlp_start_z;
    std::vector< double >  track_dlp_end_x;
    std::vector< double >  track_dlp_end_y;
    std::vector< double >  track_dlp_end_z;
    std::vector< int >     track_dlp_pdg;
    std::vector< double >  track_dlp_ixn;

    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    std::vector< int >     event;
    std::vector< int >     run;
    std::vector< int >     subrun;
    std::vector<std::string> caf_file_name;
};

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

    // Create the data object
    InteractionData data;

    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fCafTree=new TTree("CafTree", "CAF Variables");
    fCafTree->Branch("track_dlp_e", &data.track_dlp_e);
    fCafTree->Branch("track_dlp_px", &data.track_dlp_px);
    fCafTree->Branch("track_dlp_py", &data.track_dlp_py);
    fCafTree->Branch("track_dlp_pz", &data.track_dlp_pz);
    fCafTree->Branch("track_dlp_pmag", &data.track_dlp_pmag);
    fCafTree->Branch("track_dlp_length", &data.track_dlp_length);
    fCafTree->Branch("track_dlp_angle", &data.track_dlp_angle);
    fCafTree->Branch("track_dlp_angle_rot", &data.track_dlp_angle_rot);
    fCafTree->Branch("track_dlp_angle_incl", &data.track_dlp_angle_incl);
    fCafTree->Branch("track_dlp_angle_x", &data.track_dlp_angle_x);
    fCafTree->Branch("track_dlp_angle_y", &data.track_dlp_angle_y);
    fCafTree->Branch("track_dlp_angle_z", &data.track_dlp_angle_z);
    fCafTree->Branch("track_dlp_start_x", &data.track_dlp_start_x);
    fCafTree->Branch("track_dlp_start_y", &data.track_dlp_start_y);
    fCafTree->Branch("track_dlp_start_z", &data.track_dlp_start_z);
    fCafTree->Branch("track_dlp_end_x", &data.track_dlp_end_x);
    fCafTree->Branch("track_dlp_end_y", &data.track_dlp_end_y);
    fCafTree->Branch("track_dlp_end_z", &data.track_dlp_end_z);
    fCafTree->Branch("track_dlp_pdg", &data.track_dlp_pdg);
    fCafTree->Branch("track_dlp_ixn", &data.track_dlp_ixn);

    fCafTree->Branch("spill_index", &data.spill_index);
    fCafTree->Branch("file_index", &data.file_index);
    fCafTree->Branch("event", &data.event);
    fCafTree->Branch("run", &data.run);
    fCafTree->Branch("subrun", &data.subrun);
    fCafTree->Branch("caf_file_name", &data.caf_file_name);


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

                std::vector<PartInfo> parts;
                for(unsigned long ipart = 0; ipart < reco_ixn.part.dlp.size(); ++ipart)
                {
                    const auto& part = reco_ixn.part.dlp[ipart];
                    parts.push_back({part.pdg,
                                    TVector3(part.start.x, part.start.y, part.start.z),
                                    TVector3(part.end.x,   part.end.y,   part.end.z)});
                }

                TVector3 vtx(reco_ixn.vtx.x, reco_ixn.vtx.y, reco_ixn.vtx.z);

                // TODO: fill from the matched MINERvA track for this interaction
                double minerva_track_end_z = get_minerva_track_end_z(sr, ixn);

                // If the interaction doesn't pass cuts, skip it
                if(!cuts(parts, vtx, minerva_track_end_z, minerva_back_z, max_other))
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
                    data.track_dlp_e.push_back(part.E);
                    data.track_dlp_px.push_back(part.p.x);
                    data.track_dlp_py.push_back(part.p.y);
                    data.track_dlp_pz.push_back(part.p.z);
                    data.track_dlp_pmag.push_back(pvec.Mag());
                    data.track_dlp_length.push_back(length);
                    data.track_dlp_angle.push_back(dir.Angle(beam_dir));
                    data.track_dlp_angle_x.push_back(dir.Angle(x_plus_dir));
                    data.track_dlp_angle_y.push_back(dir.Angle(y_plus_dir));
                    data.track_dlp_angle_z.push_back(dir.Angle(z_plus_dir));
                    dir.RotateY(-TMath::Pi()/2);
                    data.track_dlp_angle_rot.push_back(dir.Theta());
                    data.track_dlp_angle_incl.push_back(dir.Phi());
                    data.track_dlp_start_x.push_back(part.start.x);
                    data.track_dlp_start_y.push_back(part.start.y);
                    data.track_dlp_start_z.push_back(part.start.z);
                    data.track_dlp_end_x.push_back(part.end.x);
                    data.track_dlp_end_y.push_back(part.end.y);
                    data.track_dlp_end_z.push_back(part.end.z);
                    data.track_dlp_pdg.push_back(part.pdg);
                    data.track_dlp_ixn.push_back(ixn);
                    data.spill_index.push_back(spill_num);
                    data.file_index.push_back(file_num);
                    data.event.push_back(sr->meta.nd_lar.event);
                    data.subrun.push_back(sr->meta.nd_lar.subrun);
                    data.caf_file_name.push_back(current_file.erase(0, current_file.find_last_of("/")+1).c_str());

                } // End of particle loop

                // Fill at the end of every interaction. 
                // One entry = One interaction
                fCafTree->Fill();
            } // End of interaction loop

        } // End of spill loop
        delete sr;
        fCafTree->AutoSave("Overwrite");
        caf_file->Close();
        
    } // End of file loop

    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    // Output TTree file name
    std::string file_name = "2x2_multip_output";

    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");

    // POPULATE: Write to output ROOT file
    fCafTree->Write("", TObject::kOverwrite);
        
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