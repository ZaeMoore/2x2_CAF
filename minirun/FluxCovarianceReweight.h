#ifndef FLUXCOVARIANCEREWEIGHT_H
#define FLUXCOVARIANCEREWEIGHT_H

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TH1D.h"

#include <vector>

// -------------------------------------
struct FluxBin {
    int pdg;
    double elo;
    double ehi;
    bool isRHC;  // true = RHC, false = FHC

    bool Contains(int p, double e, bool rhc) const {
        return (p == pdg && e >= elo && e < ehi && rhc == isRHC);
    }
};

// -------------------------------------
class FluxCovarianceReweight {
public:
    FluxCovarianceReweight() {}

    // Binning & covariance
    bool LoadBinning(const char* binfile);
    bool LoadCovariance();
    void AddDiagonalEpsilon(double rel = 1e-8);

    // Throw generation
    bool GenerateThrows(int nThrows, int seed = 0);
    bool GenerateThrows2(int nThrows, int seed = 0);
    TMatrixD ReconstructCovarianceFromThrows() const;
    std::vector<double> GetRelativeUncertaintyFromCovariance() const;
    // Accessors
    double GetWeightFast(int ithrow, int binIndex) const;
    int    GetBinIndex(int pdg, double Enu) const;
    int    FindBinIndex(int pdg, double Enu) const;
    double GetWeight(int ithrow, int pdg, double Enu) const;
    const TVectorD& GetAllWeights(int pdg, double Enu) const;

    double GetFluxWeight(int fluxBin, int universe) const;
    // Nominal flux
    const TVectorD& GetNominalFlux() const { return fNominalFlux; }
    const std::vector<double>& GetFluxBinCenter() const { return fFluxBinCenter; }

    // Integrated flux
    double GetIntegratedFluxNominal() const { return fIntegratedFluxNominal; }
    const std::vector<double>& GetIntegratedFluxThrows() const {
        return fIntegratedFluxThrows;
    }
    void ComputeIntegratedFluxes();

    // Load 8 histograms for FHC & RHC
    void LoadFluxHistograms(
        TH1D* h_fhc_nue,
        TH1D* h_fhc_nuebar,
        TH1D* h_fhc_numu,
        TH1D* h_fhc_numubar,
        TH1D* h_rhc_nue,
        TH1D* h_rhc_nuebar,
        TH1D* h_rhc_numu,
        TH1D* h_rhc_numubar
    );

    // Build flux vectors
    std::vector<double> BuildFluxBinCenterVector(
        TH1D* h_fhc_nue,
        TH1D* h_fhc_nuebar,
        TH1D* h_fhc_numu,
        TH1D* h_fhc_numubar,
        TH1D* h_rhc_nue,
        TH1D* h_rhc_nuebar,
        TH1D* h_rhc_numu,
        TH1D* h_rhc_numubar
    ) const;

    TVectorD BuildNominalFluxVector(
        TH1D* h_fhc_nue,
        TH1D* h_fhc_nuebar,
        TH1D* h_fhc_numu,
        TH1D* h_fhc_numubar,
        TH1D* h_rhc_nue,
        TH1D* h_rhc_nuebar,
        TH1D* h_rhc_numu,
        TH1D* h_rhc_numubar
    ) const;

    void FillAllThrowsForBin(int binIndex, Double_t* out) const;
    bool LoadCombinedCovariance();
private:
    std::vector<FluxBin> fBins;       // Binning including FHC/RHC
    std::vector<TVectorD> fThrows;    // Generated throws
    TMatrixDSym fCov;                 // Covariance matrix

    TVectorD fNominalFlux;            // Nominal flux (628 bins)
    std::vector<double> fFluxBinCenter; // Bin centers (628 bins)

    double fIntegratedFluxNominal;    // Integrated nominal flux
    std::vector<double> fIntegratedFluxThrows; // Integrated flux per universe
};

#endif