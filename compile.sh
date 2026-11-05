#compile the CAF plotter
#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` minirun/FluxCovarianceReweight.h -o FluxCovarianceReweight.o -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict


#g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` minirun/caf_analysis.cxx -o minirun_analysis -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict
#if [ $? -eq 0 ]; then
#    echo "Compilation successful"
#else
#    echo "Compilation failed"
#fi

# 1. Compile the flux class to an object file (-c: compile only, no linking, no main needed)
g++ -std=c++17 -O2 -g -pedantic -Wall -c `root-config --cflags` \
    minirun/FluxCovarianceReweight.cc -o FluxCovarianceReweight.o
if [ $? -ne 0 ]; then
    echo "Flux class compilation failed"
    return 1
fi

# 2. Compile the analysis and link the flux object file into it
g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` \
    minirun/caf_analysis.cxx FluxCovarianceReweight.o -o minirun_analysis \
    -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC \
    -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -lduneanaobj_StandardRecord -lduneanaobj_StandardRecord_dict
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi