"""

Created on Fri Nov  8 02:08:45 2019

 

@author: Aaron Sinclair

 

Reads xyz format file containing sets of co-ordinates defining multiple crystal lattices

Each atom has close neighbours (atoms within a certain distance) and a quantum q6 vector.

For each atom, calculates the number of close neighbours and their q6 vector

For each lattice, calculates the number of atoms with >= 7 close neighbours and a q6 vector > 0.5

Lattices with >= 95% consituent atoms with close neighbours can be classified as solid

Latices with <= 5% consituent atoms with close neighbours can be classified as liquid

Lattices may fail to be identified as solid or liquid if the data is incorrect or the chosen length scale



"""

 

import scipy as sc

import numpy as np

from scipy import special

import pickle

 

#function to read the xyz file

#this returns 500 xyz coordinates in a list of lists with 3 variables in each

 

L = 7.93701     # ?? Define

l = 6           # ?? Define

neighbour_dist = 1.5   # ?? Define

neighbour_dist_sqrd = neighbour_dist**2

 

 

def reduce_vector(i,j,L):

    """ Summary: ???

   

    Parameters

    ----------

    i : ???

        ???

    j : ??

        ???

    i : List tuple

        ???

       

    Returns

    ----------

    ???

        ????

    """   

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

    """ Summary: Obtains q6 vector for each atom3

   

    Parameters

    ----------

    coordinates : List tuple

        (x,y,z) co-ordinates for atomes in lattice

    Returns

    ----------

    List

        ????

    """

   

    a = 0

    parameters = []

 

 

#gives list of all spherical harmonics for each particle

 

    for a in range(len(coordinates)):

                  

        i = coordinates[a]

        s_harm = []

        b = 0

               

        

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

#returns the spherical harmonics particle which are neighbours to the particle    

                            

            else:

                b += 1

                if b >= len(coordinates):

                    break

                j = coordinates[b]

                               

        #if statement acts as a check, as b cannot be larger than coordinatesement  

                                

        param = (sum(s_harm)) / (len(s_harm))

        s = param.tolist()

        parameters.append(s)

        a+=1

       

    # store the data as binary data stream

    # ab instead of wb, so it doesnt overwrite parameters       

    with open('solid.parameters.txt', 'ab') as filehandle:

 

          

        pickle.dump(parameters, filehandle)

           

    

    return parameters

 

 

 

def normalised_parameters(parameters):

    """ Summary: returns normalised q6 vectors

   

    Parameters

    ----------

    parameters : q6 vectors for each atom 

        ???

       

    Returns

    ----------

    normalised q6 vectors

        ????

    """ 

    

    norm_p = [[i / np.linalg.norm(j) for i in j] for j in parameters]

    norm_array=np.array([np.array(n) for n in norm_p])

    return norm_array

 

#function normalises each parameter and returns it as a array of arrays

 

def phase_finder(norm_param,coordinates,L):

    """ Counts number atoms with >7 near neighbours atoms in the same crystal

    

    Atoms is considered near neightbour if within a threshold distance

   

    Parameters

    ----------

    norm_param : normailised q6 vector

        

    coordinates : xyz coordinates of the atoms in each lattice

        

    L : subatomic length constant

        

       

    Returns

    ----------

    The number of atoms in the lattice that can be considered solid

          

    """

   

    print ("Phase_finder")           

def read_xyz(filename):
    
    """ Summary: Finds how many particles/500 can be considered solid

   

    Parameters

    ----------

    filename : String

        Path to file in xyz format defining multiple crystal lattices

       

    Returns

    ----------

    number of atoms that can be considered solid in each lattice

          

    """   


    
    n_lattices = 0

    m=[]

   

    xyz = open(filename)

   

    # Iterate through the co-ordinate sets in xyz file

   

    for line in xyz:

       

        n_lattices += 1

 

        # Start reading data for next lattice

        # Assume next line contains number of atoms in the lattice

       

        atoms = []

        coordinates = []

   

        n_atoms = int(line)

       

        # Assume next line is the title

       

        title = xyz.readline().rstrip('\n')

        print("title = {}, solid atom count = {}".format(title,n_atoms))

 

   

        # Read the next n_atoms lines from the file - assume these are the co-ordinates for the atoms the lattice

       

        for count in range(n_atoms):

           

            line = xyz.readline()

           

            atom,x,y,z = line.split()

            atoms.append(atom)

            coordinates.append([float(x), float(y), float(z)])

       

        p = obtain_parameters(coordinates)

 

        norm = normalised_parameters(p)

 

        sc = phase_finder(norm,coordinates,L)

        

        m.append((title, sc))

 

    xyz.close()       

    print( "{} lattices processed".format(n_lattices))

   

    return m

 

#t= read_xyz("solid.xyz")

t= read_xyz("solid.xyz")

 

print(t)
