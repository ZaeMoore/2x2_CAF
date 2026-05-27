source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cmake v3_22_2
setup gcc v12_1_0
setup pycurl
setup ifdhc
setup geant4 v4_11_0_p01c -q e20:debug
setup dk2nugenie   v01_10_01k -q debug:e20
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau
setup jobsub_client
setup eigen v3_3_5
setup duneanaobj v03_06_01b -q e20:prof
setup srproxy v00.43 -q py3913
setup hdf5 v1_10_5a -q e20
setup fhiclcpp v4_15_03 -q debug:e20
#setup root v6_28_12 -q e26:p3915:prof 
setup root v6_26_06b -q e20:p3913:prof
#export DISPLAY=0:0


