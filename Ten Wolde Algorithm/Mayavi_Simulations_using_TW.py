
import scipy as sc
import numpy as np
import pickle
import scipy.special
from mayavi import mlab

#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

L = 7.93701
l = 6
neighbour_dist = 1.15
neighbour_dist_sqrd = neighbour_dist**2
min_count = 7
threshold_dist = 0.5


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
  
def obtain_parameters(coordinates):

    
    mrange = range(-l,l+1)

    parameters = []
    neighbours = []

#gives list of all spherical harmonics for each particle 

# DQ using enumerate is more "pythonic" than that you were doing....
    for a, i in enumerate(coordinates):
        s_harm = []
        neighbours_a = []
        
#moved b = 0 outside bracket,put back if it messes thinbgs        
    
        for b, j in enumerate(coordinates):
            
            if a != b:
                
                r = reduce_vector(i,j,L)

                if np.dot(r,r) <= neighbour_dist_sqrd:
                    
                    neighbours_a.append(b)

#   Appends everyneighbour of an atom to a list 
                                    
                    the = np.arctan2(r[1],r[0])
                    phi = np.arccos(r[2]/np.linalg.norm(r))
                    
                    m_vec = np.array([sc.special.sph_harm(m,l,the,phi) for m in mrange ])
                    
                    s_harm.append(m_vec)

#   Returns the spherical harmonics particle which are neighbours to the particle the spherical harmonics particle which are neighbours to the particle    
                               
        param = (sum(s_harm)) / (len(s_harm))
        s = param.tolist()
        parameters.append(s)

        neighbours.append(neighbours_a)
        
#   Appends list of neighbours of each particle to a list 
        
        
    #with open('d.solid.parameters.txt', 'ab') as filehandle:
# store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
        #pickle.dump(parameters, filehandle)
            
    
    return parameters, neighbours

#f


def normalised_parameters(parameters):
    norm_p = [[i / np.linalg.norm(j) for i in j] for j in parameters]
    norm_array=np.array([np.array(n) for n in norm_p])
    return norm_array

#function normalises each parameter and returns it as a array of arrays

def phase_finder(norm_param,coordinates,neighbours,threshold_dist,min_count):
    
    solid_count = 0 
    liquid_coordin = []
    solid_coordin = []
    
    for a, i in enumerate(norm_param):
        
        n_a = neighbours[a]
        count = 0
        
        for b, q in enumerate(n_a):
            
            j = norm_param[q]
        
                    
            if np.vdot(j, i)>threshold_dist:
                count +=1

        
        if count > min_count:
            
            solid_count +=1
            solid_coordin.append(coordinates[a])
            
        else:
            liquid_coordin.append(coordinates[a])
            
            
    return solid_count,liquid_coordin,solid_coordin


def read_xyz(filename):

    m=[]
    
    xyz = open(filename)
    
 #   for a in range(1):    
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
        
        
        
        p,neighbours = obtain_parameters(coordinates)
        norm = normalised_parameters(p)   
        sc,liquid_coordin,solid_coordin = phase_finder(norm,coordinates,neighbours,threshold_dist,min_count)
        
        print("Processed a configuration")
        
        
        
        m.append(sc)
        
        x_solid = []
        y_solid = []
        z_solid = []
        
        for u,v,w in solid_coordin:
            x_solid.append(u)
            y_solid.append(v)
            z_solid.append(w)
        
        x_liquid = []
        y_liquid = []
        z_liquid = []
  
        for f,g,h in liquid_coordin:
            x_liquid.append(f)
            y_liquid.append(g)
            z_liquid.append(h)
    
        edge_x = [L/2]
        edge_y = [L/2]
        edge_z = [L/2]
    
        mlab.points3d(edge_x, edge_y, edge_z, scale_factor=L,color=(0.5,0.5,0.5), mode='cube',resolution=12,opacity=0.4)
    
    #produces a single transparent cube defining the boundaries crystal sturcture
    
        mlab.points3d(x_liquid, y_liquid, z_liquid, scale_factor=0.98,color=(0.1,0.1,0.9), mode='sphere',resolution=12,opacity=1)
        mlab.points3d(x_solid, y_solid, z_solid, scale_factor=0.98,color=(0.9,0.1,0.1), mode='sphere',resolution=12,opacity=1)
        mlab.show()
    
    xyz.close()        
    return m

t= read_xyz("liquid.xyz")
print(t)

