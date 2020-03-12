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
l2 = 4
neighbour_dist = 1.5/(rho)**(1/3)
neighbour_dist_sqrd = neighbour_dist**2


def reduce_vector(i,j,L):
    
# constantly converting lists to numpy arrays is slow. Instead I've stored
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
  
def obtain_parameters(coordinates):

    
    mrange = range(-l,l+1)
    nrange = range(-l2,l2+1)
    parameters = []


# gives list of all spherical harmonics for each particle 


    for a, i in enumerate(coordinates):
        s_harm1 = []
        s_harm2 = []
        
        conc = []
    
    
        for b, j in enumerate(coordinates):
            
            if a != b:
                
                r = reduce_vector(i,j,L)

                if np.dot(r,r) <= neighbour_dist_sqrd:
                    
                    the = np.arctan2(r[1],r[0])
                    phi = np.arccos(r[2]/np.linalg.norm(r))
                    
                    m_vec = np.array([sc.special.sph_harm(m,l,the,phi) for m in mrange ])
                    
                    
                    n_vec = np.array([sc.special.sph_harm(m,l2,the,phi) for m in nrange ])
                    
                    # gives q4 spherical harmonic aswell
                    
                    s_harm1.append(m_vec)
                    s_harm2.append(n_vec)

#returns the spherical harmonics particle which are neighbours to the particle the spherical harmonics particle which are neighbours to the particle    
                               
        param1 = (sum(s_harm1))/ (len(s_harm1))
        
        param2 = (sum(s_harm2))/ (len(s_harm2))
 
#param1 and param2 are both 1d arrays (lengths 13 and 9), which contain the q6 and q4 set of vectors respectively, for one particle    
    
        s = np.append(param1,param2)
        
#s is a 1d array with 22 compnents   
    
        s1 = s.tolist()

#s1 iss in list form        
        
        parameters.append(s1)

#parameters will have 500 1d arrays containing the 22 componet descriptor fo each particle in the configuration        


    with open('solid_T1.0_rho1.00_loc_inv.txt', 'ab') as filehandle:

        # store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
       pickle.dump(parameters, filehandle)
  
    return parameters


#This function repeats this process over the number of configurations. 

#timings were done using the time.time() function

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
        
        
        
        p = obtain_parameters(coordinates)
        
        print("Processed a configuration")
        
        m.append(sc)
        

        
    xyz.close()        
    return m

t= read_xyz("solid_T1.0_rho1.00.xyz")

print(t)

