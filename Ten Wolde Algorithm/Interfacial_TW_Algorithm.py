import scipy as sc
import numpy as np
import pickle
import scipy.special


#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

rho = 0.965
L_x = 7.8559306
L_y = 7.8559306
L_z = 67.14331

#Periodic boundary conditions were not cubic, so this was specified in the components of L

l = 6
neighbour_dist = 1.5/(rho)**(1/3)
neighbour_dist_sqrd = neighbour_dist**2


def reduce_vector(i,j,L_x,L_y,L_z):
    
# constantly converting lists to numpy arrays is slow. Instead I've stored
# the lists as numpy arrays right from the start.
    
    #vi = np.array(i) 
    #vj = np.array(j)            
    o = j - i
    
    d_x = o[0]/L_x
    d_y = o[1]/L_y
    d_z = o[2]/L_z
    
    d = np.array([d_x,d_y,d_z])
    
                
# This is faster than a loop as it can be vectorized in numpy
    d=d-np.floor(d+0.5)  

#this corrects for the atoms clsoe to edge of box           
        
    r_x = d[0]*L_x
    r_y = d[1]*L_y
    r_z = d[2]*L_z     
    
    r = np.array([r_x,r_y,r_z])
    
    return r
   
#creates reduced vector                    
  
def obtain_parameters(coordinates):

    
    mrange = range(-l,l+1)

    parameters = []


#gives list of all spherical harmonics for each particle 


    for a, i in enumerate(coordinates):
        s_harm = []
        
#moved b = 0 outside bracket,put back if it messes thinbgs        
    
        for b, j in enumerate(coordinates):
            
            if a != b:
                
                r = reduce_vector(i,j,L_x,L_y,L_z)

                if np.dot(r,r) <= neighbour_dist_sqrd:
                    
                    the = np.arctan2(r[1],r[0])
                    phi = np.arccos(r[2]/np.linalg.norm(r))
                    
                    m_vec = np.array([sc.special.sph_harm(m,l,the,phi) for m in mrange ])
                    
                    s_harm.append(m_vec)

#returns the spherical harmonics particle which are neighbours to the particle the spherical harmonics particle which are neighbours to the particle    
                               
        param = (sum(s_harm)) / (len(s_harm))
        s = param.tolist()
        parameters.append(s)


    with open('traj.txt', 'ab') as filehandle:
# store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
        pickle.dump(parameters, filehandle)
            
    
    return parameters



def normalised_parameters(parameters):
    norm_p = [[i / np.linalg.norm(j) for i in j] for j in parameters]
    norm_array=np.array([np.array(n) for n in norm_p])
    return norm_array

#function normalises each parameter and returns it as a array of arrays

def phase_finder(norm_param,coordinates,L_x,L_y,L_z):
    
    a = 0
    solid_count = 0 
    
    for a, i in enumerate(norm_param):
        
        ci = coordinates[a]
        count = 0
        
        for b, j in enumerate(norm_param):
            
            cj = coordinates[b]
        
            if a != b:
                
                r = reduce_vector(ci,cj,L_x,L_y,L_z)
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
    w = 1
    xyz = open(filename)
    
 #   for a in range(1):    
    while (True):
    
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
        sc = phase_finder(norm,coordinates,L_x,L_y,L_z)
        
        print("Processed a configuration")
        
        w += 1
        
        m.append(sc)
        
  
    xyz.close()        
    return m

#This function repeats the first 3 functions over the number of configurations. 
#Given that there is only 1 configuration, the while condition is set to execute only once.
#timings were done using the time.time() function

t= read_xyz("traj.xyz")
print(t)

