import ase

import quippy
import numpy as np
import pickle
import scipy.special
from ase.io import read

#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

rho = 1
L = 3/(rho)**(1/3)
l = 6
neighbour_dist = 4.25/(rho)**(1/3)


def soap_descriptors(filename):    

    my_atoms = read(filename)

    #my_atoms.set_cell([L, L, L]) # Assume cubic box - dont want this

    #my_atoms.pbc=True # Periodic boundary conditions

# Define the descriptor we want to calculate

#to only calculate SOAP descriptors for a given central atom: Z=32
#to make SOAP distinguish between two neighbours: n_species=2 species_Z={32 52}


    desc = quippy.descriptors.Descriptor("soap cutoff=4.25 cutoff_transition_width=0.0 atom_sigma=0.7 n_species=2 Z=32 species_Z={32 52} n_max=8 l_max=6")
    #
    
    #for descriptors only for one atom type, use line 58 instead, with extra z component specified
    #desc = quippy.descriptors.Descriptor("soap cutoff={} cutoff_transition_width=0.0 atom_sigma=0.7 Z=32 n_species=2 species_Z={32 52} n_max=8 l_max=6".format(neighbour_dist))
# Pass the my_atoms object to this descriptor

    descriptors = desc.calc(my_atoms)["data"]


    #with open('soap(n_max=8)_GeTe_alpha.txt', 'ab') as filehandle:
       
       #store the data as binary data stream
       #ab instead of wb, so it doesnt overwrite parameters
           
     #  pickle.dump(descriptors, filehandle)

    return descriptors

def read_xyz(filename):

    start2 = time.time()
    

    i = 0
    xyz = open(filename)  

    while i < 1:

        p = soap_descriptors(filename)

        
        print("Processed a configuration")
        

        
        i += 1 
        
    xyz.close()        
    
    end2 = time.time()

    print(end2-start2)
    
    return p

t= read_xyz("quenched_99.xyz")






