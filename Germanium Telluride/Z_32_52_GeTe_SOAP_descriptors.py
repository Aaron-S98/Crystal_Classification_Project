import ase

import quippy
import numpy as np
import pickle
import scipy.special
from ase.io import read


rho = 1
L = 3/(rho)**(1/3)
l = 6
neighbour_dist = 4.25/(rho)**(1/3)

#This function reads the ASE file and calculates the SOAP descriptor throuh the interface provided.
#This function has 2 extra components that specify the Z of the atoms involved and the Z value of the descriptor you wan to produce


def soap_descriptors(filename):    

    my_atoms = read(filename)

    #my_atoms.set_cell([L, L, L]) # Assume cubic box - dont want this

    #my_atoms.pbc=True # Periodic boundary conditions

# Define the descriptor we want to calculate

#to only calculate SOAP descriptors for a given central atom: Z=32
#to make SOAP distinguish between two neighbours: n_species=2 species_Z={32 52}


    desc = quippy.descriptors.Descriptor("soap cutoff=4.25 cutoff_transition_width=0.0 atom_sigma=0.7 n_species=2 Z=32 species_Z={32 52} n_max=8 l_max=6")
    
    
    
# Pass the my_atoms object to this descriptor

    descriptors = desc.calc(my_atoms)["data"]


    with open('soap(n_max=8)_GeTe_alpha.txt', 'ab') as filehandle:
       
       #stores the data as binary data stream
       #ab instead of wb, so it doesnt overwrite parameters
           
       pickle.dump(descriptors, filehandle)

    return descriptors

#This function repeats this process over the number of configurations. 
#Given that there is only 1 configuration, the while condition is set to execute only once.
#timings were done using the time.time() function

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






