import numpy as np

# sigmoid function to normalize inputs
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# sigmoid derivatives to adjust synaptic weights
def sigmoid_derivative(x):
    return x * (1 - x)

# Import input layers from the time history data files
import pandas as pd
f1=open('sa.txt','r')
df1 = pd.read_csv(f1, sep="\s+", header=0, names=['x','y', 'acc'])
#time=df.loc[:,"time"]
acc=df1.loc[:5,"acc"] 
print acc

# Import training output from a txt file    
f2=open('pid.txt','r')
df2 = pd.read_csv(f2, sep="\s+", header=0, names=['x','y', 'pid'])
#time=df.loc[:,"time"]
pid=df2.loc[:5,"pid"] 
# input dataset
#training_inputs = np.array([[0,0,1],
#                            [1,1,1],
#                            [1,0,1],
#                            [0,1,1]])

# output dataset
#training_outputs = np.array([[0,1,1,0]]).T
training_inputs = np.array([acc]).T
#print training_inputs
             
training_outputs = np.array([pid]).T
print training_outputs
# seed random numbers to make calculation
np.random.seed(1)

# initialize weights randomly with mean 0 to create weight matrix, synaptic weights
synaptic_weights = 2 * np.random.random((1,len(acc))) - 1

print('Random starting synaptic weights: ')
print(synaptic_weights)

# Iterate 10,000 times
for iteration in range(10000):

    # Define input layer
    input_layer = training_inputs
    # Normalize the product of the input layer with the synaptic weights
    outputs = sigmoid(np.dot(input_layer, synaptic_weights))

    # how much did we miss?
    error = training_outputs - outputs

    # multiply how much we missed by the
    # slope of the sigmoid at the values in outputs
    adjustments = error * sigmoid_derivative(outputs)

    # update weights
    synaptic_weights += np.dot(input_layer.T, adjustments)

print('Synaptic weights after training: ')
print(synaptic_weights)

print("Output After Training:")
print(outputs)
