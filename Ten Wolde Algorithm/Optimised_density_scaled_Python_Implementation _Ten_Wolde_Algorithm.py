import scipy as sc
import numpy as np
import pickle
import scipy.special


#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop


# density scaling is applied to L and neighbour_dist parameters


rho = 1.0
L = 7.93701/(rho)**(1/3)
l = 6

neighbour_dist = 1.5/(rho)**(1/3)
neighbour_dist_sqrd = neighbour_dist**2

def reduce_vector(i,j,L):
    
# constantly converting lists to numpy arrays is slow. Instead I've stored
# the lists as numpy arrays right from the start.
    
    #vi = np.array(i) 
    #vj = np.array(j)            
    r = j - i
    d = r/L
                
# This is faster than a loop as it can be vectorized in numpy
    d=d-np.floor(d+0.5)  

#this corrects for the atoms clsoe to edge of box           
        
    r = d*L
     
    return r
   
#creates reduced vector                    
  
def obtain_parameters(coordinates):

    
    mrange = range(-l,l+1)

    parameters = []


#gives list of all spherical harmonics for each particle 


    for a, i in enumerate(coordinates):
        s_harm = []
             
    
        for b, j in enumerate(coordinates):
            
            if a != b:
                
                r = reduce_vector(i,j,L)

                if np.dot(r,r) <= neighbour_dist_sqrd:
                    
                    the = np.arctan2(r[1],r[0])
                    phi = np.arccos(r[2]/np.linalg.norm(r))
                    
                    m_vec = np.array([sc.special.sph_harm(m,l,the,phi) for m in mrange ])
                    
                    s_harm.append(m_vec)

#returns the spherical harmonics particle which are neighbours to the particle the spherical harmonics particle which are neighbours to the particle    
                               
        param = (sum(s_harm)) / (len(s_harm))
        s = param.tolist()
        parameters.append(s)


    with open('misalignedq30_test_q6.txt', 'ab') as filehandle:

           
        pickle.dump(parameters, filehandle)

# stores the data as binary data stream     
    
    return parameters

#This function obtains the q6 vector for each particle in the configuration


def normalised_parameters(parameters):
    

    norm_p = [[i / np.linalg.norm(j) for i in j] for j in parameters]
    norm_array=np.array([np.array(n) for n in norm_p])
    
    return norm_array

#function normalises each q6 vector and returns it as a array of arrays

def phase_finder(norm_param,coordinates,L):
    
    a = 0
    solid_count = 0 
    
    for a, i in enumerate(norm_param):
        
        ci = coordinates[a]
        count = 0
        
        for b, j in enumerate(norm_param):
            
            cj = coordinates[b]
        
            if a != b:
                
                r = reduce_vector(ci,cj,L)
                distance_sqd = np.dot(r,r)
                
                if distance_sqd<=neighbour_dist_sqrd:
                    
                    if np.vdot(j, i)>0.5:
                        count +=1

        
        if count > 7:
            solid_count +=1
    
    return solid_count

#This function checks if each particle has 7 conncetions.
#If true, the particle is set as solid-like and appended to list
#If false, the particle is set a liquid-like

def read_xyz(filename):

    m=[]
    
    xyz = open(filename)
    i = 0
 #   for a in range(1):    
    while(True):
    #while i<1:
        
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
        norm = normalised_parameters(p)   
        sc = phase_finder(norm,coordinates,L)
        
        print("Processed a configuration")
        
        m.append(sc)
        
        i+=1
    xyz.close()        
    return m

#This function repeats the first 3 functions over the number of configurations. 
#Given that there is only 1 configuration, the while condition is set to execute only once.
#timings were done using the time.time() function

t= read_xyz("q30.xyz")

print(t)

