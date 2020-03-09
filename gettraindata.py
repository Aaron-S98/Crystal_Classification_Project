import numpy as np 
import pickle

a = []
b = []


with open("conc_q4_q6_liquid.txt", 'rb') as fr:
    try:
        x = 0
        while True:
            
            a.append(pickle.load(fr))
            x += 1
    except EOFError:
        pass
    fr.close()

with open("conc_q4_q6_solid.txt", 'rb') as fr:
    try:
        y = 0
        while True:
           
            b.append(pickle.load(fr))
            y += 1
    except EOFError:
        pass

#Gives parameters of each atom for 101 configurations as the training data
        
A = flatten = [item for sublist in a for item in sublist]
B = flatten = [item for sublist in b for item in sublist]

#flattens list of list into one list



#converts list into 1d array
solid_input = []
liquid_input = []

for element in A:
    
    component_A = []
    for config in element:
    
        A_real = config.real
        A_imag = config.imag
        component_A.append(A_real)
        component_A.append(A_imag)
    
    liquid_input.append(component_A)    

print (len(liquid_input))

for element in B:
    
    component_B = []
    for config in element:
    
        B_real = config.real
        B_imag = config.imag
        component_B.append(B_real)
        component_B.append(B_imag)
        
    solid_input.append(component_B)    
        
print (len(solid_input))

    
#seperates real and complex parts and appends them to the 2D
    
output_liquid = np.array(liquid_input)
output_solid = np.array(solid_input)
   
#converts into 1d array

train = []
test = []

liq_train = output_liquid[0:25500]
sol_train = output_solid[0:25500]
train.append(liq_train)
train.append(sol_train)

#   This slices half the liquid and solid arrays for training

liq_test = output_liquid[25500:]
sol_test = output_solid[25500:]
test.append(liq_test)
test.append(sol_test)

#   This slices half the liquid and solid arrays for testing
train_list = [item for sublist in train for item in sublist]
test_list = [item for sublist in test for item in sublist]

train_input = np.array(train_list)
test_input = np.array(test_list) 

print(len(train_input))
print(len(test_input))

np.save('"train_conc_q4_q6_solid_liquid',train_input)
np.save('test_conc_q4_q6_solid_liquid',test_input)

#creating labels

a = np.zeros(25500)
b = np.ones(25500)

train_labels = np.append(a,b)

m = np.zeros(25000)
n = np.ones(25000)

test_labels = np.append(m,n)

np.save('train_conc_q4_q6_solid_liquid_labels',train_labels)
np.save('test_conc_q4_q6_solid_liquid_labels',test_labels)




