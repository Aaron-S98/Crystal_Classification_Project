import ase
import time
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

    desc = quippy.descriptors.Descriptor("soap cutoff={} cutoff_transition_width=0.0 atom_sigma=0.7 n_max=9 l_max=6".format(neighbour_dist))

# Pass the my_atoms object to this descriptor

    descriptors = desc.calc(my_atoms)["data"]


  #  with open('soap(n_max=9)_GeTe_quenched.txt', 'ab') as filehandle:
       
       #store the data as binary data stream
       #ab instead of wb, so it doesnt overwrite parameters
           
  #     pickle.dump(descriptors, filehandle)
        
    return descriptors


def read_xyz(filename):

    
    time1 = time.time()
    
    i=0
    xyz = open(filename)  

    while i < 1:

        p = soap_descriptors(filename)

        
        print("Processed a configuration")
        
        i+=1
    
    
    xyz.close()        
    
    time2 = time.time()
    
    print(time2-time1)
    
    return p

t= read_xyz("quenched_9.xyz")






