export mywd=/global/u2/r/rdiurba/nusystematics

setup boost v1_80_0 -q e26:prof
setup tbb v2021_7_0 -q e26
setup sqlite v3_39_02_00

export fhiclcpp_ROOT=${mywd}/fhicl-cpp-standalone/build/fhicl-cpp-install/
export PATH=${fhiclcpp_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${fhiclcpp_ROOT}/lib/:${LD_LIBRARY_PATH}

export cetlib_ROOT=${mywd}/fhicl-cpp-standalone/build/cetlib-install/
export PATH=${cetlib_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${cetlib_ROOT}/lib/:${LD_LIBRARY_PATH}

export cetlib_except_ROOT=${mywd}/fhicl-cpp-standalone/build/cetlib-except-install/
export PATH=${cetlib_except_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${cetlib_except_ROOT}/lib/:${LD_LIBRARY_PATH}

export hep_concurrency_ROOT=${mywd}/fhicl-cpp-standalone/build/hep-concurrency-install
export PATH=${hep_concurrency_ROOT}/bin/:${PATH}
export LD_LIBRARY_PATH=${hep_concurrency_ROOT}/lib/:${LD_LIBRARY_PATH}


source ${NUSYST}/build/Linux/bin/setup.nusystematics.sh
export PATH=/global/u2/r/rdiurba/nusystematics/ndnusyst/ndnusyst-src/build/Linux/bin/:${PATH}