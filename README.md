# 2x2_CAF
Analysis of 2x2 CAF files

Main script is minirun/caf_analysis.cxx \
You must run the compiled script from the topmost directory, which also must contain numiFluxSystFull.root
and flux_covariance_binning_NuMI_GeV.txt \
minirun directory must also contain FluxCovarianceReweight.cc and FluxCovarianceReweight.h

DUNE GPVMS:

Run setup.sh \
Run compile.sh \
./minirun_analysis 2> errors.log


NERSC:
Enter shifter container "shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash" \
Run setup.sh \
Run compile.sh \
./minirun_analysis 2> errors.log

