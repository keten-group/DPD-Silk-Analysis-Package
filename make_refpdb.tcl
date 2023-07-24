package require topotools
topo readlammpsdata ref.data bond
set protein [atomselect top "type 1 2 4 5 6"]
$protein writepdb ref.pdb
exit