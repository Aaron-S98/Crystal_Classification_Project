import scipy as sc
import numpy as np
import pickle
import scipy.special
import itertools
import py3nj as py
import time

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
A1 = (4*np.pi)/(2*l+1)
A2 = (4*np.pi)/(2*l2+1)
A1sqrt = A1**(1/2)
A2sqrt = A2**(1/2)

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
    nrange = range(-l2,l2+1)
    parameters = []
    
# gives list of all spherical harmonics for each particle 

# DQ using enumerate is more "pythonic" than that you were doing....
    for a, i in enumerate(coordinates):
        s_harm1 = []
        s_harm2 = []
        
        conc = []
# moved b = 0 outside bracket,put back if it messes things        
    
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
#there are 500 q6m values in param1 and param2

        local_inv_q6 = A1sqrt*np.linalg.norm(param1)
        local_inv_q4 = A2sqrt*np.linalg.norm(param2)
        
       #calculate w_descriptor      

        result1 = [seq for i in range(len(mrange), 0, -1) for seq in itertools.combinations(range(-l,l+1), 3) if sum(seq) == 0]
        result2 = [seq for i in range(len(nrange), 0, -1) for seq in itertools.combinations(range(-l2,l2+1), 3) if sum(seq) == 0]

        #produces w(i) matrix as list of tuples, however there are repeating tuples in the list.

        out1 = set([i for i in result1])
        out2 = set([i for i in result2])

        #remove duplicates from list of tuples
        
        

        w_l = list(out1)
        w_l2 = list(out2)

        #coverts set into list of unique tuples
       
        list_of_w_l = []
        list_of_w_l2 = []

 

        for i, j in enumerate(w_l):
    
            m1 = j[0]
            m2 = j[1]
            m3 = j[2]
    
            qlm1 = param1[m1+6]
            qlm2 = param1[m2+6]
            qlm3 = param1[m3+6]
    
    #corrects for index position with m values as it starts from m=-l
    
            x = py.wigner3j(l*2, l*2, l*2, j[0]*2, j[1]*2, j[2]*2)
    
            x1 = x*qlm1*qlm2*qlm3
    
            list_of_w_l.append(x1)
    
# appends each combination of each wl(i)*qlm1*qlm2*qlm3 matrix to a list


        for i, j in enumerate(w_l2):
    
            m1 = j[0]
            m2 = j[1]
            m3 = j[2]
    
            ql2m1 = param2[m1+4]
            ql2m2 = param2[m2+4]
            ql2m3 = param2[m3+4]
    
    #corrects for index position with m values as it starts from m=-l
    
            y = py.wigner3j(l2*2, l2*2, l2*2, j[0]*2, j[1]*2, j[2]*2)
            y1 = y*ql2m1*ql2m2*ql2m3
    
            list_of_w_l2.append(y1)
    
  
# appends each combination of each wl2(i)*qlm1*qlm2*qlm3 matrix to a list

        w_l_loc = sum(list_of_w_l)
        w_l2_loc = sum(list_of_w_l2)

# w_l(i) parameter
        
        w_l_hat_re_im = (w_l_loc)/(np.linalg.norm(param1))**3
        w_l2_hat_re_im = (w_l2_loc)/(np.linalg.norm(param2))**3
        
# get both real and imaginary parts of the w hat loc invaraint descriptor.        
        
        w_l_hat = np.real(w_l_hat_re_im)
        w_l2_hat = np.real(w_l2_hat_re_im)


        li = np.append(local_inv_q6,local_inv_q4)
        
        li_list = li.tolist()
        
        q6q4w6w4 = li_list + [w_l_hat, w_l2_hat]
        
        #appends all 4 descriptors to 1 list for each particle
        
        parameters.append(q6q4w6w4)
        
    
#appends all descriptors to 1 list with 4 descriptors for each of the 500 particles


    with open('q6q4w6w4_loc_inv_misaligned_quench.txt', 'ab') as filehandle:

        # store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
        pickle.dump(parameters, filehandle)
  
    return parameters

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
        
        
        
        p = obtain_parameters(coordinates)
        
        print("Processed a configuration")
        
        m.append(sc)
        
        
    time2 = time.time()
    
    print(time2-time1)
    
    xyz.close()        
    return m

t= read_xyz("solid.xyz")

