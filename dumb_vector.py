import scipy as sc
import numpy as np
import pickle
import scipy.special


#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop


rho = 1
L = 7.93701/(rho)**(1/3)
l = 6
neighbour_dist = 1.8/(rho)**(1/3)
neighbour_dist_sqrd = neighbour_dist**2


def reduce_vector(i,j,L):
    
# DQ constantly converting lists to numpy arrays is slow. Instead I've stored
# the lists as numpy arrays right from the start.
    
    #vi = np.array(i) 
    #vj = np.array(j)            
    r = j - i
    d = r/L
                
# DQ: This is faster than a loop as it can be vectorized in numpy
    d=d-np.floor(d+0.5)  

#this corrects for the atoms clsoe to edge of box           
        
    r = d*L
     
    return r
   
#creates reduced vector                    
  
def dumb_descriptors(coordinates):

    descriptors = []

#gives list of all spherical harmonics for each particle 
# DQ using enumerate is more "pythonic" than that you were doing....
    
    for a, i in enumerate(coordinates):
        
        p = 12            
        total_vectors = []       
        total_vectors_magnitudes = []       
        neighbour_vectors = []

#these are vectors with positive and negative values, so want to make into magnitudes
      
#this puts the indicies of the first 12 values of the magnitude of the r vlaues.

        for b, j in enumerate(coordinates):
            
            if a != b:
                
                r = reduce_vector(i,j,L)

                y = np.dot(r,r) 
                #y is magnitude squared
                
                if y <= neighbour_dist_sqrd:
                    
                    total_vectors.append(r)
                    total_vectors_magnitudes.append(y)
       
#appends dot products to a list, so magnitudes can be compared        
        indicies_needed = np.argpartition(total_vectors_magnitudes, p)[0:12]

                
#Gives the index position of the first 12 values as an array 

        for i in indicies_needed:
                 
            neighbour_vectors.append(total_vectors[i])  
        
        descriptors.append(neighbour_vectors)

    with open('dumb_vector_quenched.txt', 'ab') as filehandle:
       
#store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
        pickle.dump(descriptors, filehandle)       

    return descriptors

def read_xyz(filename): 

    m=[]
    
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
        
        p = dumb_descriptors(coordinates)

        
        print("Processed a configuration")
        
        m.append(sc)
        
    xyz.close()        
    return m

t= read_xyz("solid_T0.2_rho1.00.xyz")
