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


    
}