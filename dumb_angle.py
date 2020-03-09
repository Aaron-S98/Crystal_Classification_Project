import scipy as sc
import numpy as np
import pickle
import scipy.special
import time

#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

L = 7.93701
l = 6
neighbour_dist = 1.8
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
        
        dumb_a_total = []
        
        mag_dumb_a_total = [] 
        
        dumb_a_12 = []
          
        angles = []

#these are vectors with positive and negative values, so want to make into magnitudes
      

        
#this puts the indicies of the first 12 values of the magnitude of the r vlaues.
        for b, j in enumerate(coordinates):
            
                if a != b:
                
                    r = reduce_vector(i,j,L)

                    y = np.dot(r,r) 
                #y is magnitude squared
                
                    if y <= neighbour_dist_sqrd:
                    
                        dumb_a_total.append(r)
                        mag_dumb_a_total.append(y)
       
#appends dot products to a list, so magnitudes can be compared        
                
        indicies_needed = np.argpartition(mag_dumb_a_total, p)
        indicies_needed_12 = indicies_needed[0:12]
                
#Gives the index position of the first 12 values as an array 

        for i in indicies_needed_12:
            
            
            d_a_12 = dumb_a_total[i]  
            
            dumb_a_12.append(d_a_12)
        
        s = dumb_a_12
    
        ref_particle = dumb_a_12[0]
        rest_of_particles = dumb_a_12[1:]
        
        for j in rest_of_particles:
        
            v = reduce_vector(ref_particle,j,L)
            
            the = np.arctan2(v[1],v[0])
            phi = np.arccos(v[2]/np.linalg.norm(v))
            angles.append(the)
            angles.append(phi)

        descriptors.append(angles)

 #   with open('dumb_angle_misaligned_quenched.txt', 'ab') as filehandle:
       
#store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
 #       pickle.dump(descriptors, filehandle)
            
    
    return descriptors


def read_xyz(filename): 
    
    time1 = time.time()
    
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
    
    time2 = time.time()
    
    print(time2-time1)
    
    xyz.close()        
    return m

t= read_xyz("q1.xyz")
