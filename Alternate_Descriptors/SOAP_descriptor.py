import ase
import quippy
import numpy as np
import pickle
import scipy.special
import time

#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop


rho = 1.2
L = 7.93701/(rho)**(1/3)
l = 6
neighbour_dist = 1.5/(rho)**(1/3)


def soap_descriptors(coordinates):    

    coordinates_array = np.stack(coordinates)

    labels = ['Ar' for x in range(0,500)]



    my_atoms = ase.Atoms(labels,positions=coordinates_array)

    #Regarding the parameters of SOAP:

        #cutoff: as it says on the tin

        #cutoff_transition_width: range within the cutoff smoothly changes from 1 to 0 (leave at 0 for us)

       # atom_sigma: width of Gaussian representing each neighbour

        #n_max: number of radial basis functions (radial resolution)

        #l_max: number of angular basis functions (angular resolution), bandwidth of spherical harmonics
    
    my_atoms.set_cell([L, L, L]) # Assume cubic box

    my_atoms.pbc=True # Periodic boundary conditions

# Define the descriptor we want to calculate

    desc = quippy.descriptors.Descriptor("soap cutoff={} cutoff_transition_width=0.0 atom_sigma=0.25 n_max=9 l_max=6".format(neighbour_dist))

# Pass the my_atoms object to this descriptor

    descriptors = desc.calc(my_atoms)["data"]


   # with open('soap(n_max=9)_liquid(T=2.1)_test.txt', 'ab') as filehandle:
       
       #store the data as binary data stream
       #ab instead of wb, so it doesnt overwrite parameters
           
    #   pickle.dump(descriptors, filehandle)
    

   
    return descriptors


def read_xyz(filename):

    start2 = time.time()
        
    xyz = open(filename)

   
    while(True):

        
        atoms = []
        coordinates = []
    
        try:
            n_atoms = int(xyz.readline())
        except:
            print("Reached end of file")
            break
        
        title = xyz.readline()
    
        for line in xyz:
            
            atom,x,y,z = line.split()
            atoms.append(atom)
            coordinates.append(np.array([float(x), float(y), float(z)]))
            if len(coordinates) == n_atoms:
                break
            
   
        
        p = soap_descriptors(coordinates)
        
        print("Processed a configuration")
        

        
    end2 = time.time()

    print(end2-start2)
    
    xyz.close()        
    return p

t= read_xyz("liquid_T2.1rho1.0.xyz")

