import scipy as sc
import numpy as np
import pickle
import scipy.special
import ase.spacegroup as spg
import ase.neighborlist as nl
from ase.io import read
#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

rho = 1
L = 7.93701/(rho)**(1/3)
l = 6
neighbour_dist = 4.25/(rho)**(1/3)
neighbour_dist_sqrd = neighbour_dist**2
                  
  
def obtain_parameters(file):

    atomsGeTe=read(file)
    symbols = atomsGeTe.get_chemical_symbols()

    FirstAtom, SecondAtom, vects = nl.neighbor_list(['i','j','D'], atomsGeTe, 4.25, self_interaction=False)

    cell = atomsGeTe.get_cell()
    newvects = nl.mic(vects, cell, pbc=[True, True, True])

    
    descriptors = []


# gives list of all spherical harmonics for each particle 

# DQ using enumerate is more "pythonic" than that you were doing....
    for j in range(len(symbols)):
        
        indices = [c for c,atom in enumerate(FirstAtom) if atom==j]
        vectors = np.array([newvects[i] for i in indices])
        
        
        ord1_fn_1 = []
        ord1_fn_2 = []
       
        ord2_fn_1 = []
        ord2_fn_2 = []
        ord2_fn_3 = []
        
        ord3_fn_1 = []
        ord3_fn_2 = []
        ord3_fn_3 = []
        ord3_fn_4 = []
        
        ord4_fn_1 = []
        ord4_fn_2 = []
        ord4_fn_3 = []
        ord4_fn_4 = []
        ord4_fn_5 = []
        
        ord5_fn_1 = []
        ord5_fn_2 = []
        ord5_fn_3 = []
        ord5_fn_4 = []
        ord5_fn_5 = []
        ord5_fn_6 = []
        
        b_ord4_fn_1 = []
        
# moved b = 0 outside bracket,put back if it messes things        

        for a, i in enumerate(vectors):         
                
                mag_vector = np.linalg.norm(i)
                
                x,y,z = i[0]/mag_vector,i[1]/mag_vector,i[2]/mag_vector
                    
                o1_f1 = 0.5*x**1 + 0.866025*y**1
                o1_f2 = 1*z**1
                     
                o2_f1 = 0.540062*x**2 + -0.801784*x**1*y**1 + 0.0771517*y**2 + -0.617213*z**2
                o2_f2 = 0.92582*x**1*y**1 + 0.534522*y**2 + -0.534522*z**2
                o2_f3 = 0.707107*x**1*z**1 + 1.22474*y**1*z**1 

                o3_f1 = 0.53619*x**3 + 0.121136*x**2*y**1 + -1.32882*x**1*y**2 + 0.121136*y**3 + -0.279751*x**1*z**2 + -0.484544*y**1*z**2
                o3_f2 = 0.312772*x**2*y**1 + 0.722315*x**1*y**2 + 0.312772*y**3 + -0.722315*x**1*z**2 + -1.25109*y**1*z**2
                o3_f3 = 1.12916*x**2*z**1 + -1.15045*x**1*y**1*z**1 + 0.464948*y**2*z**1 + -0.531369*z**3
                o3_f4 = 1.78227*x**1*y**1*z**1 + 1.02899*y**2*z**1 + -0.342997*z**3

                o4_f1 = + 0.285044*x**4 + 0.542539*x**3*y**1 + -0.432264*x**2*y**2 + -0.97657*x**1*y**3 + 0.15975*y**4 + -1.278*x**2*z**2 + 1.30209*x**1*y**1*z**2 + -0.526235*y**2*z**2 + 0.300706*z**4
                o4_f2 = + 1.19161*x**3*y**1 + -0.893343*x**2*y**2 + -0.63434*x**1*y**3 + 0.16087*y**4 + 0.893343*x**2*z**2 + -1.67181*x**1*y**1*z**2 + -0.0718782*y**2*z**2 + -0.136911*z**4
                o4_f3 = + 1.14953*x**3*z**1 + 0.48431*x**2*y**1*z**1 + -2.33014*x**1*y**2*z**1 + 0.48431*y**3*z**1 + -0.372822*x**1*z**3 + -0.645746*y**1*z**3
                o4_f4 = + 0.518321*x**2*y**2 + 0.598506*x**1*y**3 + 0.172774*y**4 + -0.518321*x**2*z**2 + -1.79552*x**1*y**1*z**2 + -1.55496*y**2*z**2 + 0.345547*z**4
                o4_f5 = + 0.854242*x**2*y**1*z**1 + 1.97279*x**1*y**2*z**1 + 0.854242*y**3*z**1 + -0.657596*x**1*z**3 + -1.13899*y**1*z**3                     
                    
                o5_f1 = + 0.240391*x**5 + -0.509292*x**4*y**1 + -0.876962*x**3*y**2 + 1.23302*x**2*y**3 + -0.077379*x**1*y**4 + -0.0589707*y**5 + -1.52695*x**3*z**2 + -0.643317*x**2*y**1*z**2 + 3.09516*x**1*y**2*z**2 + -0.643317*y**3*z**2 + 0.247613*x**1*z**4 + 0.428878*y**1*z**4
                o5_f2 = + 0.96686*x**4*y**1 + 0.964265*x**3*y**2 + -1.72842*x**2*y**3 + -0.727203*x**1*y**4 + 0.234432*y**5 + -0.964265*x**3*z**2 + -0.615905*x**2*y**1*z**2 + 1.47042*x**1*y**2*z**2 + -0.615905*y**3*z**2 + 0.237062*x**1*z**4 + 0.410603*y**1*z**4
                o5_f3 = + 0.900562*x**4*z**1 + 0.400687*x**3*y**1*z**1 + -0.0495722*x**2*y**2*z**1 + -2.00344*x**1*y**3*z**1 + 0.437888*y**4*z**1 + -1.7846*x**2*z**3 + 1.60275*x**1*y**1*z**3 + -0.859252*y**2*z**3 + 0.264385*z**5
                o5_f4 = + 0.17967*x**3*y**2 + 0.518662*x**2*y**3 + 0.419229*x**1*y**4 + 0.103732*y**5 + -0.17967*x**3*z**2 + -1.55599*x**2*y**1*z**2 + -3.05439*x**1*y**2*z**2 + -1.55599*y**3*z**2 + 0.598899*x**1*z**4 + 1.03732*y**1*z**4
                o5_f5 = + 3.13679*x**3*y**1*z**1 + -2.06432*x**2*y**2*z**1 + -1.33807*x**1*y**3*z**1 + 0.519245*y**4*z**1 + 0.688106*x**2*z**3 + -1.79872*x**1*y**1*z**3 + -0.350385*y**2*z**3 + -0.0337721*z**5
                o5_f6 = + 1.77394*x**2*y**2*z**1 + 2.04837*x**1*y**3*z**1 + 0.591312*y**4*z**1 + -0.591312*x**2*z**3 + -2.04837*x**1*y**1*z**3 + -1.77394*y**2*z**3 + 0.236525*z**5
                     
                b_o4_f1 = + 0.365148*x**4 + -1.09545*x**2*y**2 + 0.365148*y**4 + -1.09545*x**2*z**2 + -1.09545*y**2*z**2 + 0.365148*z**4
                    
                   # symmetry adapted functions will be calculated for one neighour
                    
                ord1_fn_1.append(o1_f1)
                ord1_fn_2.append(o1_f2)
                    
                ord2_fn_1.append(o2_f1)
                ord2_fn_2.append(o2_f2)
                ord2_fn_3.append(o2_f3)
                    
                ord3_fn_1.append(o3_f1)
                ord3_fn_2.append(o3_f2)
                ord3_fn_3.append(o3_f3)
                ord3_fn_4.append(o3_f4)
                    
                ord4_fn_1.append(o4_f1)
                ord4_fn_2.append(o4_f2)
                ord4_fn_3.append(o4_f3)
                ord4_fn_4.append(o4_f4)
                ord4_fn_5.append(o4_f5)
                    
                ord5_fn_1.append(o5_f1)
                ord5_fn_2.append(o5_f2)
                ord5_fn_3.append(o5_f3)
                ord5_fn_4.append(o5_f4)
                ord5_fn_5.append(o5_f5)
                ord5_fn_6.append(o5_f6)
                    
                b_ord4_fn_1.append(b_o4_f1)
       
                    # this will append each set of values to a list 

#This returns the sum of all the neighbours for each symmetry adapted function                               
        
        s11 = sum(ord1_fn_1) 
        s12 = sum(ord1_fn_2) 
       
        s21 = sum(ord2_fn_1)
        s22 = sum(ord2_fn_2)
        s23 = sum(ord2_fn_3)
        
        s31 = sum(ord3_fn_1) 
        s32 = sum(ord3_fn_2) 
        s33 = sum(ord3_fn_3) 
        s34 = sum(ord3_fn_4)
        
        s41 = sum(ord4_fn_1) 
        s42 = sum(ord4_fn_2) 
        s43 = sum(ord4_fn_3) 
        s44 = sum(ord4_fn_4) 
        s45 = sum(ord4_fn_5)
        
        s51 = sum(ord5_fn_1) 
        s52 = sum(ord5_fn_2) 
        s53 = sum(ord5_fn_3) 
        s54 = sum(ord5_fn_4) 
        s55 = sum(ord5_fn_5) 
        s56 = sum(ord5_fn_6)
    
        sb1 = sum(b_ord4_fn_1)
        

#param1 and param2 are both 1d arrays (lengths 13 and 9), which contain the q6 and q4 set of vectors respectively, for one particle    
    
        component_1 = np.sqrt(s11**2 + s12**2)
        component_2 = np.sqrt(s21**2 + s22**2 + s23**2)
        component_3 = np.sqrt(s31**2 + s32**2 + s33**2 + s34**2)
        component_4 = np.sqrt(s41**2 + s42**2 + s43**2 + s44**2 + s45**2)
        component_5 = np.sqrt(s51**2 + s52**2 + s53**2 + s54**2 + s55**2 + s56**2)
        component_6 = np.sqrt(sb1**2)
        
        
        descr = [component_1,component_2,component_3,component_4,component_5,component_6]

        
# appending all of the descriptors for          
        
        descriptors.append(descr)

#parameters will have 500 1d arrays containing the 22 componet descriptor fo each particle in the configuration        


   # with open('symmetry_adapted_GeTe_quenched.txt', 'ab') as filehandle:

        # store the data as binary data stream
#ab instead of wb, so it doesnt overwrite parameters
           
     #   pickle.dump(descriptors, filehandle)
            
       
  
    return descriptors

def read_xyz(filename):


    j = 1
    xyz = open(filename)
    
  
    while j<=1:
    
        atoms = []
        coordinates = []
    
        try:
            n_atoms = int(xyz.readline())
        except:
            print("Reached end of file")
            break
       
        p = obtain_parameters(filename)

        
        print("Processed a configuration")
        
        j += 1
        
    xyz.close()        
    

    return p

t= read_xyz("quenched_99.xyz")


