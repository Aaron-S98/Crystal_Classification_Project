import numpy as np 
import pickle
import scipy as sc

a = []
b = []
c = []

with open("quench_locinv_run1-20.txt", 'rb') as fr:
    try:
       x = 0
       while True:
            
            a.append(pickle.load(fr))
            x += 1
    except EOFError:
        pass
    fr.close()

with open("solid_locinv.txt", 'rb') as fr:
    try:
        y = 0
        while True:
           
            b.append(pickle.load(fr))
            y += 1
    except EOFError:
        pass

with open("liquid_locinv.txt", 'rb') as fr:
    try:
        z = 0
        while True:
           
            c.append(pickle.load(fr))
            z += 1
    except EOFError:
        pass    
    
                            
#Gives parameters of each atom for 101 configurations as the training data
        
A = flatten = [item for sublist in a for item in sublist]
B = flatten = [item for sublist in b for item in sublist]
C = flatten = [item for sublist in c for item in sublist]

#flattens list of list into one list



#converts list into 1d array
quenched_input = []

solid_input = []
liquid_input = []

for element in A:
    
    component_A = []
    for config in element:
    
        A_real = config.real
        A_imag = config.imag
        component_A.append(A_real)
        component_A.append(A_imag)
    
    quenched_input.append(component_A)    

print (len(quenched_input))

for element in B:
    
    component_B = []
    for config in element:
    
        B_real = config.real
        B_imag = config.imag
        component_B.append(B_real)
        component_B.append(B_imag)
        
    solid_input.append(component_B)    
        
print (len(solid_input))

for element in C:
    
    component_C = []
    for config in element:
    
        C_real = config.real
        C_imag = config.imag
        component_C.append(C_real)
        component_C.append(C_imag)
        
    liquid_input.append(component_C)    
        
print (len(liquid_input))


#seperates real and complex parts and appends them to the 2D
    
output_quenched = np.array(quenched_input)
output_solid = np.array(solid_input)
output_liquid = np.array(liquid_input)
   
#converts into 1d array

train = []
test = []

quenched_train = output_quenched[0:10000]
sol_train = output_solid[0:5000]
liq_train = output_liquid[0:50000]

train.append(quenched_train)
train.append(sol_train)
train.append(liq_train)

#   This slices a 3rd of the quenched, solid and liquid arrays for training


sol_test = output_solid[5000:1000]
liq_test = output_liquid[12500:25000]

quenched_test = output_quenched[12500:25000]


test.append(quenched_test)


test.append(sol_test)
test.append(liq_test)

#   This slices half the liquid and solid arrays for testing

train_list = [item for sublist in train for item in sublist]
test_list = [item for sublist in test for item in sublist]

#breaks up parameters from a list of list to one long list

train_input = np.array(train_list)
test_input = np.array(test_list) 

print(len(train_input))
print(len(test_input))

#np.save('train_input_liquid_locinv.npy',train_input)
#np.save('test_input_conc_q4_q6_multiple_tempertures_quenched_solid_liquid.npy',test_input)

#creating labels

a = np.ones(55000)
b = np.zeros(50000)


train_labels = np.append(a,b)

print(len(train_labels))

m = np.ones(17000)
n = np.zeros(17000)

test_labels = np.append(m,n)

#np.save('train_input_liquid_locinv_labels.npy',train_labels)
#np.save('test_input_conc_q4_q6_multiple_tempertures_quenched_solid_liquid_labels.npy',test_labels)


print(len(test_labels))

