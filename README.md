# Crystal_Classification_Project_1stdraft

Reads xyz format file containing sets of co-ordinates defining multiple crystal lattices

Each atom has close neighbours (atoms within a certain distance) and a quantum q6 vector.

For each atom, calculates the number of close neighbours and their q6 vector

For each lattice, calculates the number of atoms with >= 7 close neighbours and a q6 vector > 0.5

Lattices with >= 95% consituent atoms with close neighbours can be classified as solid

Latices with <= 5% consituent atoms with close neighbours can be classified as liquid

Lattices may fail to be identified as solid or liquid if the data is incorrect 
