#include "FluxCovarianceReweight.h"

#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TDecompChol.h"
#include "TMatrixDSymEigen.h"
#include "TH1D.h"

#include <fstream>
#include <iostream>
#include <cmath>

// =====================================================
// Load binning
// =====================================================
bool FluxCovarianceReweight::LoadBinning(const char* binfile)
{
    std::ifstream fin("flux_covariance_binning_NuMI_GeV.txt");
    if (!fin.is_open()) {
        std::cerr << "ERROR: cannot open binning file\n";
        return false;
    }

    int isRHC, pdg;
    double elo, ehi;

    while (fin >> isRHC >> pdg >> elo >> ehi) {
        FluxBin b;
        b.pdg = pdg;
        b.elo = elo;
        b.ehi = ehi;
        b.isRHC = isRHC;
        fBins.push_back(b);
    }

    std::cout << "Loaded " << fBins.size() << " flux bins\n";
    return true;
}

// =====================================================
// Load covariance matrix
// =====================================================
bool FluxCovarianceReweight::LoadCovariance()
{

    TFile* f = TFile::Open("numiFluxSyst.root", "READ");
    if (!f || f->IsZombie()) return false;
    TH2D* h=nullptr;
    f->GetObject("hcov_hadProd_abs", h);
    TH1D* hflux_fhc_nue=nullptr;
    TH1D* hflux_fhc_nuebar=nullptr;
    TH1D* hflux_fhc_numu=nullptr;
    TH1D* hflux_fhc_numubar=nullptr;    
    
    TH1D* hflux_rhc_nue=nullptr;
    TH1D* hflux_rhc_nuebar=nullptr;
    TH1D* hflux_rhc_numu=nullptr;
    TH1D* hflux_rhc_numubar=nullptr;
    f->GetObject("hflux_fhc_nue_ppfx_corrected", hflux_fhc_nue);
    f->GetObject("hflux_fhc_nuebar_ppfx_corrected", hflux_fhc_nuebar);
    f->GetObject("hflux_fhc_numu_ppfx_corrected", hflux_fhc_numu);
    f->GetObject("hflux_fhc_numubar_ppfx_corrected", hflux_fhc_numubar);

    f->GetObject("hflux_rhc_nue_ppfx_corrected", hflux_rhc_nue);
    f->GetObject("hflux_rhc_nuebar_ppfx_corrected", hflux_rhc_nuebar);
    f->GetObject("hflux_rhc_numu_ppfx_corrected", hflux_rhc_numu);
    f->GetObject("hflux_rhc_numubar_ppfx_corrected", hflux_rhc_numubar);

    LoadBinning("flux_covariance_binning_NuMI_GeV.txt");

    LoadFluxHistograms(
      hflux_fhc_nue,
      hflux_fhc_nuebar,
      hflux_fhc_numu,
      hflux_fhc_numubar,
      hflux_rhc_nue,
      hflux_rhc_nuebar,
      hflux_rhc_numu,
      hflux_rhc_numubar
    );

    if (!h) {
        std::cerr << "ERROR: could not find hcov_total_abs\n";
        return false;
    }

    const int n = h->GetNbinsX();
    std::cout<<fBins.size()<<","<<n<<std::endl;
    if ((int)fBins.size() != n) {
        std::cerr << "ERROR: binning size does not match covariance\n";
        return false;
    }

    if (fNominalFlux.GetNrows() != n) {
        std::cerr << "ERROR: nominal flux not loaded before covariance\n";
        return false;
    }

    fCov.ResizeTo(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            const double cov_abs = h->GetBinContent(i+1, j+1);

            const double phi_i = fNominalFlux[i];
            const double phi_j = fNominalFlux[j];

            if (phi_i > 0.0 && phi_j > 0.0) {
                fCov(i,j) = cov_abs / (phi_i * phi_j);
            } else {
                fCov(i,j) = 0.0;
            }
        }
    }

    return true;
}
bool FluxCovarianceReweight::LoadCombinedCovariance()
{
    // Load hadron production covariance (FHC + RHC)
    TFile* f = TFile::Open("numiFluxSystFull.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open numiFluxSyst.root\n";
        return false;
    }
    
    TH2D* h_hadprod = nullptr;
    f->GetObject("hcov_hadProd_abs", h_hadprod);
    if (!h_hadprod) {
        std::cerr << "ERROR: could not find hcov_total_abs\n";
        return false;
    }
    

    
    TH2D* h_beamfocus = nullptr;
    f->GetObject("hcov_beamFocus_abs", h_beamfocus);
    if (!h_beamfocus) {
        std::cerr << "ERROR: could not find cov_beam_focusing_rhc\n";
        return false;
    }
    
    // Load flux histograms
    TH1D* hflux_fhc_nue=nullptr;
    TH1D* hflux_fhc_nuebar=nullptr;
    TH1D* hflux_fhc_numu=nullptr;
    TH1D* hflux_fhc_numubar=nullptr;    
    
    TH1D* hflux_rhc_nue=nullptr;
    TH1D* hflux_rhc_nuebar=nullptr;
    TH1D* hflux_rhc_numu=nullptr;
    TH1D* hflux_rhc_numubar=nullptr;
    f->GetObject("hflux_fhc_nue_ppfx_corrected", hflux_fhc_nue);
    f->GetObject("hflux_fhc_nuebar_ppfx_corrected", hflux_fhc_nuebar);
    f->GetObject("hflux_fhc_numu_ppfx_corrected", hflux_fhc_numu);
    f->GetObject("hflux_fhc_numubar_ppfx_corrected", hflux_fhc_numubar);

    f->GetObject("hflux_rhc_nue_ppfx_corrected", hflux_rhc_nue);
    f->GetObject("hflux_rhc_nuebar_ppfx_corrected", hflux_rhc_nuebar);
    f->GetObject("hflux_rhc_numu_ppfx_corrected", hflux_rhc_numu);
    f->GetObject("hflux_rhc_numubar_ppfx_corrected", hflux_rhc_numubar);

    
    // Load binning and flux
    LoadBinning("flux_covariance_binning_NuMI_GeV.txt");
    LoadFluxHistograms(
        hflux_fhc_nue,
        hflux_fhc_nuebar,
        hflux_fhc_numu,
        hflux_fhc_numubar,
        hflux_rhc_nue,
        hflux_rhc_nuebar,
        hflux_rhc_numu,
        hflux_rhc_numubar
    );
    
    const int n = h_hadprod->GetNbinsX();
    std::cout << "Total bins: " << fBins.size() << ", Covariance size: " << n << std::endl;
    
    if ((int)fBins.size() != n) {
        std::cerr << "ERROR: binning size does not match covariance\n";
        return false;
    }
    
    if (fNominalFlux.GetNrows() != n) {
        std::cerr << "ERROR: nominal flux not loaded before covariance\n";
        return false;
    }
    
    // Check beam focusing covariance size (should be 314 bins for RHC)
    const int n_rhc = h_beamfocus->GetNbinsX();
    const int n_fhc = n - n_rhc;  // First half is FHC
    
    if (n_rhc * 2 != n) {
        std::cerr << "WARNING: Expected beam focusing bins (" << n_rhc 
                  << ") to be half of total (" << n << ")\n";
    }
    
    std::cout << "FHC bins: " << n_fhc << ", RHC bins: " << n_rhc << std::endl;
    
    // Resize and fill combined covariance
    fCov.ResizeTo(n, n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Start with hadron production covariance
            double cov_abs = h_hadprod->GetBinContent(i+1, j+1);
            
            // Add beam focusing contribution for RHC bins only
            // RHC bins are in the second half: [n_fhc, n)
            if (i >= n_fhc && j >= n_fhc) {
                int i_rhc = i - n_fhc;  // 0-indexed RHC bin
                int j_rhc = j - n_fhc;
                cov_abs += h_beamfocus->GetBinContent(i_rhc+1, j_rhc+1);
            }
            
            // Convert to fractional covariance
            const double phi_i = fNominalFlux[i];
            const double phi_j = fNominalFlux[j];
            
            if (phi_i > 0.0 && phi_j > 0.0) {
                fCov(i, j) = cov_abs / (phi_i * phi_j);
            } else {
                fCov(i, j) = 0.0;
            }
        }
    }
    
    f->Close();

    
    std::cout << "Successfully loaded and combined covariances\n";
    return true;
}
// =====================================================
// Generate correlated throws
// =====================================================
bool FluxCovarianceReweight::GenerateThrows(int nThrows, int seed)
{
    LoadCombinedCovariance(); 
    TRandom3 rng(seed);
    //AddDiagonalEpsilon();
    TDecompChol chol(fCov);
    if (!chol.Decompose()) {
        std::cerr << "ERROR: Cholesky failed\n";
        return false;
    }

    TMatrixD* U = (TMatrixD*)chol.GetU().Clone();
    U->Transpose(chol.GetU()); // L = U^T

    int n = fCov.GetNrows();
    fThrows.clear();
    fThrows.reserve(nThrows);

    for (int t = 0; t < nThrows; ++t) {
        TVectorD z(n);
        for (int i = 0; i < n; ++i)
            z(i) = rng.Gaus(0,1);

        TVectorD delta = *U * z;

        // Convert to reweights: w = 1 + delta
        for (int i = 0; i < n; ++i)
            delta(i) = 1.0 + delta(i);

        fThrows.push_back(delta);
    }

    std::cout << "Generated " << nThrows << " throws\n";
    ComputeIntegratedFluxes();

    TMatrixD test=(TMatrixD) ReconstructCovarianceFromThrows();

    for (int i = 100; i < n; ++i){
    //std::cout<<i<<","<<std::sqrt(fCov(i,i))<<std::endl;

    }

    return true;
}

// =====================================================
// Find bin index
// =====================================================
int FluxCovarianceReweight::FindBinIndex(int pdg, double Enu) const
{
    for (size_t i = 0; i < fBins.size(); ++i)
        if (fBins[i].Contains(pdg, Enu,1))
            return i;

    return -1;
}

// =====================================================
// Public lookup
// =====================================================
double FluxCovarianceReweight::GetWeight(int ithrow,
                                          int pdg,
                                          double Enu) const
{
    if (ithrow < 0 || ithrow >= (int)fThrows.size())
        return 1.0;

    int idx = FindBinIndex(pdg, Enu);
    if (idx < 0)
        return 1.0;

    return fThrows[ithrow](idx);
}
bool FluxCovarianceReweight::GenerateThrows2(int nThrows, int seed)
{
    TRandom3 rng(seed);

    const int n = fCov.GetNrows();

    // Eigen decomposition
    TMatrixDSymEigen eig(fCov);
    TVectorD eigenValues = eig.GetEigenValues();
    TMatrixD eigenVectors = eig.GetEigenVectors();

    // Build sqrt(D) with protection
    TMatrixD sqrtD(n, n);
    sqrtD.Zero();

    int nClipped = 0;
    for (int i = 0; i < n; ++i) {
        if (eigenValues(i) > 0) {
            sqrtD(i,i) = std::sqrt(eigenValues(i));
        } else {
            sqrtD(i,i) = 0.0;  // clip negative modes
            nClipped++;
        }
    }

    if (nClipped > 0) {
        std::cout << "WARNING: clipped " << nClipped
                  << " non-positive eigenvalues\n";
    }

    // Transformation matrix
    TMatrixD A = eigenVectors * sqrtD;

    fThrows.clear();
    fThrows.reserve(nThrows);

    for (int t = 0; t < nThrows; ++t) {

        // Uncorrelated standard normals
        TVectorD z(n);
        for (int i = 0; i < n; ++i)
            z(i) = rng.Gaus(0,1);

        // Correlated variation
        TVectorD delta = A * z;

        // Convert to reweights
        for (int i = 0; i < n; ++i)
            delta(i) = 1.0 + delta(i);

        fThrows.push_back(delta);
    }
       ComputeIntegratedFluxes();

    //std::cout << "Generated " << nThrows << " throws using eigen decomposition\n";
    return true;
}
const TVectorD& FluxCovarianceReweight::GetAllWeights(int pdg,
                                                       double Enu) const
{
    static TVectorD unity;  // fallback

    int idx = FindBinIndex(pdg, Enu);
    if (idx < 0) {
        unity.ResizeTo(fThrows.size());
        unity = 1.0;
        return unity;
    }

    static TVectorD weights;
    weights.ResizeTo(fThrows.size());

    for (size_t i = 0; i < fThrows.size(); ++i)
        weights(i) = fThrows[i](idx);

    return weights;
}
void FluxCovarianceReweight::AddDiagonalEpsilon(double rel)
{
    double maxDiag = 0.0;
    for (int i = 0; i < fCov.GetNrows(); ++i)
        maxDiag = std::max(maxDiag, fCov(i,i));

    double eps = rel * maxDiag;

    for (int i = 0; i < fCov.GetNrows(); ++i)
        fCov(i,i) += eps;

    std::cout << "Added diagonal epsilon = " << eps << "\n";
}
int FluxCovarianceReweight::GetBinIndex(int pdg, double Enu) const
{
    return FindBinIndex(pdg, Enu);
}

double FluxCovarianceReweight::GetWeightFast(int ithrow,
                                               int binIndex) const
{
    if (ithrow < 0 || ithrow >= (int)fThrows.size())
        return 1.0;

    if (binIndex < 0 || binIndex >= fThrows[ithrow].GetNrows())
        return 1.0;

    return fThrows[ithrow](binIndex);
}
void FluxCovarianceReweight::FillAllThrowsForBin(int binIndex,
                                                  Double_t* out) const
{
    const int nThrows = fThrows.size();

    for (int ithrow = 0; ithrow < nThrows; ++ithrow)
        out[ithrow] = 1.0;

    if (binIndex < 0)
        return;

    for (int ithrow = 0; ithrow < nThrows; ++ithrow) {
        if (binIndex < fThrows[ithrow].GetNrows())
            out[ithrow] = fThrows[ithrow](binIndex);
    }
}
TMatrixD
FluxCovarianceReweight::ReconstructCovarianceFromThrows() const
{
    const int nUniv = fThrows.size();
    const int nBins = fThrows[0].GetNrows();

    TMatrixD covEmp(nBins, nBins);
    covEmp.Zero();

    // Compute mean weight per bin (should be ~1)
    std::vector<double> mean(nBins, 0.0);

    for (int u = 0; u < nUniv; ++u) {
        for (int i = 0; i < nBins; ++i) {
            mean[i] += fThrows[u](i);
        }
    }

    for (int i = 0; i < nBins; ++i)
        mean[i] /= nUniv;

    // Build covariance
    for (int u = 0; u < nUniv; ++u) {
        for (int i = 0; i < nBins; ++i) {
            const double di = fThrows[u](i) - mean[i];
            for (int j = 0; j < nBins; ++j) {
                const double dj = fThrows[u](j) - mean[j];
                covEmp(i,j) += di * dj;
            }
        }
    }

    covEmp *= (1.0 / (nUniv - 1));

    return covEmp;
}
std::vector<double>
FluxCovarianceReweight::GetRelativeUncertaintyFromCovariance() const
{
    const int nBins = fCov.GetNrows();
    std::vector<double> rel(nBins, 0.0);

    for (int i = 0; i < nBins; ++i) {
        if (fNominalFlux[i] > 0.0 && fCov(i,i) > 0.0) {
            rel[i] = std::sqrt(fCov(i,i)); /// fNominalFlux[i];
        } else {
            rel[i] = 0.0;
        }
    }

    for (int i = 507; i < nBins; ++i)     std::cout<<rel[i]<<std::endl;
    return rel;
}

std::vector<double> FluxCovarianceReweight::BuildFluxBinCenterVector(
    TH1D* h_fhc_nue,
    TH1D* h_fhc_nuebar,
    TH1D* h_fhc_numu,
    TH1D* h_fhc_numubar,
    TH1D* h_rhc_nue,
    TH1D* h_rhc_nuebar,
    TH1D* h_rhc_numu,
    TH1D* h_rhc_numubar
) const {
    std::vector<double> w;
    w.reserve(628); // total bins

    TH1D* hists[8] = {h_fhc_nue, h_fhc_nuebar, h_fhc_numu, h_fhc_numubar,
                      h_rhc_nue, h_rhc_nuebar, h_rhc_numu, h_rhc_numubar};

    for (int h = 0; h < 8; ++h) {
        for (int i = 1; i <= hists[h]->GetNbinsX(); ++i)
            w.push_back(hists[h]->GetBinCenter(i));
    }

    return w;
}

TVectorD FluxCovarianceReweight::BuildNominalFluxVector(
    TH1D* h_fhc_nue,
    TH1D* h_fhc_nuebar,
    TH1D* h_fhc_numu,
    TH1D* h_fhc_numubar,
    TH1D* h_rhc_nue,
    TH1D* h_rhc_nuebar,
    TH1D* h_rhc_numu,
    TH1D* h_rhc_numubar
) const {
    TVectorD fluxNom(628);
    int bin = 0;

    TH1D* hists[8] = {h_fhc_nue, h_fhc_nuebar, h_fhc_numu, h_fhc_numubar,
                      h_rhc_nue, h_rhc_nuebar, h_rhc_numu, h_rhc_numubar};

    for (int h = 0; h < 8; ++h) {
        for (int i = 1; i <= hists[h]->GetNbinsX(); ++i)
            fluxNom[bin++] = hists[h]->GetBinContent(i);
    }

    return fluxNom;
}
double FluxCovarianceReweight::GetFluxWeight(
    int binIndex,
    int ithrow
) const {
    return fThrows[ithrow](binIndex) / fNominalFlux[binIndex];
}
void FluxCovarianceReweight::LoadFluxHistograms(
    TH1D* h_fhc_nue,
    TH1D* h_fhc_nuebar,
    TH1D* h_fhc_numu,
    TH1D* h_fhc_numubar,
    TH1D* h_rhc_nue,
    TH1D* h_rhc_nuebar,
    TH1D* h_rhc_numu,
    TH1D* h_rhc_numubar
) {
    const int NBINS = 628;

    fNominalFlux.ResizeTo(NBINS);
    fFluxBinCenter.clear();
    fFluxBinCenter.reserve(NBINS);

    TH1D* hists[8] = {h_fhc_nue, h_fhc_nuebar, h_fhc_numu, h_fhc_numubar,
                      h_rhc_nue, h_rhc_nuebar, h_rhc_numu, h_rhc_numubar};

    int bin = 0;
    for (int h = 0; h < 8; ++h) {
        for (int i = 1; i <= hists[h]->GetNbinsX(); ++i, ++bin) {
            fNominalFlux[bin] = hists[h]->GetBinContent(i);
            fFluxBinCenter.push_back(hists[h]->GetBinCenter(i));
        }
    }
}



void FluxCovarianceReweight::ComputeIntegratedFluxes()
{
    const int nUniv = (int)fThrows.size();
    const int nBins = fThrows[1].GetNrows();
    //GetRelativeUncertaintyFromCovariance();
    
    // Nominal
    fIntegratedFluxNominal = 0.0;
    for (int b = 386; b < nBins; ++b){
            if (fFluxBinCenter[b]<2) continue;
        fIntegratedFluxNominal +=
            fNominalFlux[b];// * fFluxBinCenter[b];


    }
    // Universes
    fIntegratedFluxThrows.assign(nUniv, 0.0);

    for (int u = 0; u < nUniv; ++u) {
        for (int b = 386; b < nBins; ++b) {
            if (fFluxBinCenter[b]<2) continue;
          //  if (b>506 && b<527){  
           //     continue;}
            fIntegratedFluxThrows[u] +=
                fThrows[u](b)* fNominalFlux[b];// * fFluxBinCenter[b];
        }
    }
}