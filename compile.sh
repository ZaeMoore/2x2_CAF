#compile the CAF plotter
g++ -std=c++17 -O2 -g -pedantic -Wall `root-config --cflags --glibs` minirun_caf_2p2h.cxx -o plotter_minirun_2p2h -I/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_06_01b/include -I$SRPROXY_INC -L$DUNEANAOBJ_LIB -lduneanaobj_StandardRecordProxy -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord -L $DUNEANAOBJ_LIB -lduneanaobj_StandardRecord_dict
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi