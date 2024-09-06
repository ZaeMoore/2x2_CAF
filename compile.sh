#compile the CAF plotter
g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` 2p2h_caf_plotter.cxx -o plotter2p2h_nonmec -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_01_00/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi
