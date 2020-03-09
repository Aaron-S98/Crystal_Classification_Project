
import scipy as sc
import numpy as np
from scipy import special
import pickle

#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

L = 7.93701
l = 6
neighbour_dist = 1.5
neighbour_dist_sqrd = neighbour_dist**2



def reduce_vector(i,j,L):
    
    vi = np.array(i) 
    vj = np.array(j)            
    r = vj - vi
    d = r/L
                
    for x in range(len(d)):
        d[x]=d[x]-np.floor(d[x]+0.5)  

#this corrects for the atoms clsoe to edge of box           
        
        r = d*L
    return r
   
#creates reduced vector                    
  
def obtain_parameters(coordinates):

    a = 0
    parameters = []


#gives list of all spherical harmonics for each particle 

    for a in range(len(coordinates)):
        i = coordinates[a]
        s_harm = []
        b = 0
        
#moved b = 0 outside bracket,put back if it messes thinbgs        
        
        for b in range(len(coordinates)): 
        
            j = coordinates[b]
            
            if i != j:
                
                r = reduce_vector(i,j,L)

                the = np.arctan2(r[1],r[0])
                phi = np.arccos(r[2]/np.linalg.norm(r))
                m_vec = np.empty(2*l+1,dtype=np.complex128)
                k=0        
#Using tan2 function, +0 and -0 are distinct floating point numbers
        
                for m in range (-l,l+1):
                    m_vec[k] = sc.special.sph_harm(m,l,the,phi)
                    k += 1
                    
                if np.dot(r,r) <= neighbour_dist_sqrd:
                    s_harm.append(m_vec)
                    b +=1
                    if b >= len(coordinates):
                        break
                    j = coordinates[b]
                else:
                    b += 1
                    if b >= len(coordinates):
                        break
                    j = coordinates[b]
#returns the spherical harmonics particle which are neighbours to the particle the spherical harmonics particle which are neighbours to the particle    
                            
            else:
                b += 1
                if b >= len(coordinates):
                    break
                j = coordinates[b]
                                
#if statement acts as a check, as b cannot be larger than coordinatesement acts as a check, as b cannot be larger than coordinates   
                                
        param = (sum(s_harm)) / (len(s_harm))
        s = param.tolist()
        parameters.append(s)
        a+=1
    with open('solid.parameters.txt', 'ab') as filehandle:
# store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
        pickle.dump(parameters, filehandle)
            
    
    return parameters

#f


def normalised_parameters(parameters):
    norm_p = [[i / np.linalg.norm(j) for i in j] for j in parameters]
    norm_array=np.array([np.array(n) for n in norm_p])
    return norm_array

#function normalises each parameter and returns it as a array of arrays

def phase_finder(norm_param,coordinates,L):
    
   
    
    a = 0
    solid_count = 0 
    
    for a in range(len(norm_param)):
        
        i = norm_param[a]
        ci = coordinates[a]
        count = 0
        b = 0   
        
        
        for b in range(len(norm_param)):
            
            j = norm_param[b]
            cj = coordinates[b]
        
            if a != b:
                r = reduce_vector(ci,cj,L)
                distance_sqd = np.dot(r,r)
                
                if distance_sqd<=neighbour_dist_sqrd:
                    
                    if np.vdot(j, i)>0.5:
                        count +=1
                        b += 1
                    
                    if b >= len(norm_param):
                        break
                    j = norm_param[b]
                    cj = coordinates[b]
                
                else:
                    b += 1
                    if b >= len(norm_param):
                        break
                    j = norm_param[b]
                    cj = coordinates[b]
            else:
                b += 1
                if b >= len(norm_param):
                    break
                j = coordinates[b]
                cj = coordinates[b]
        
        if count > 7:
            solid_count +=1
            a += 1
        else:
            a+=1  
    
    return solid_count


def read_xyz(filename):

    m=[]
    
    xyz = open(filename)
    
    a = 0
    
    try:
    
        while(True):
    
            atoms = []
            coordinates = []
    
            n_atoms = int(xyz.readline())
            title = xyz.readline()
    
            for line in xyz:
            
                atom,x,y,z = line.split()
                atoms.append(atom)
                coordinates.append([float(x), float(y), float(z)])
                if len(coordinates) == n_atoms:
                    break
        
            p = obtain_parameters(coordinates)
            norm = normalised_parameters(p)   
            sc = phase_finder(norm,coordinates,L)
        
        
            a += 1
            m.append(sc)
        
  
        xyz.close()        
        return m
    
    except EOFError:
        pass

t= read_xyz("solid.xyz")
print(t)

