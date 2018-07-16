# Artificial Neural Network Project

This is an artificial neural network class implemented in C++ over 3 days in my Easter holiday.

## Prerequisites

To use the project on your computer you will need:

```
C++ Compiler that supports C++11/C++14
```
## Setup

1) Ensure you have a working C++11+ compiler.
2) Download the [Neural Network Folder](https://github.com/DMarke99/CPP-Projects/tree/master/NeuralNetwork).
3) Execute by compiling [main.cpp](https://github.com/DMarke99/CPP-Projects/blob/master/NeuralNetwork/main.cpp)

## Documentation

### Classes
#### Matrix

Defines a mathematical matrix.

##### Private Attributes:

r: The number of rows in the matrix.  
c: The number of colums in the matrix.  
vals: The values stored in the matrix.

##### Methods

Matrix: Initialiser for the matrix class.  

row: Accessor for number of rows.  
col: Accessor for number of columns.  
operator(): Accessor for reference to values.  
val: Accessor for values.  
operator<<: Output operator for matrices.  

identity: The identity matrix.  
one: The matrix of ones.  
rand: A matrix of random variables taken from the Uniform distribution.  
randn: A matrix of random variables taken from the Gaussian distribution.  

isSquare: Boolean determining whether a matrix is square.  
isOrthogonal: Boolean determining whether a matrix is orthogonal.  
isSymmetric: Boolean determining whether a matrix is symmetrical.  
isAntisymmetric: Boolean determining whether a matrix is antisymmetrical.  

transpose: Returns the transpose of a matrix.  
inv: Returns the inverse of a matrix.  
det: Returns the determinant of a matrix.  
trace: Returns the trace of a matrix.  
norm: Returns the Frobenius norm of a matrix.  
map: Applies a function elementwise to a matrix and returns the result.  

getRow: Returns a row/rows from a matrix.  
getCol: Returns a column/columns from a matrix.  
omitRow: Deletes a particular row from a matrix and returns the result.  
omitCol: Deletes a particular column from a matrix and returns the result.  

concat: Concatenates a set of matrices and returns the resultant matrix.  
dot: Elementwise multiplication for matrices.

operator+: Addition for matrices.  
operator-: Subtraction for matrices.  
operator*: Scalar and matrix multiplication.  
operator/: Scalar and matrix division.  
operator==: Equality for matrices.  
import: Imports a matrix from a csv file.  

#### NeuralNetwork

Defines a recurrent neural network class.

##### Private Attributes:

inputs: The number of nodes in the first layer of the ANN.  
outputs: The number of nodes in the output layer of the ANN.  
W: The weights associated with each pair of adjacent layers.  
layers: The layer scheme of the ANN.  
activation: The activation function used in the hidden layers. Choice between Sigmoid, ReLU, leaky ReLU and linear.  
classify: A boolean determining whether the ANN is used as a classifier. If a classifier, softmax is used in final layer.  

##### Methods

NeuralNetwork: Initialiser for the neural network class.  

fit: Fits the neural network to input data.
predict: Generates a prediction based on current weights.

## Contributers

This project was created entirely by myself over 3 days in my Easter break. The greatest challenge I overcame was implementing the feed forward and backpropogation in the ANN, and calculating the derivatives/gradients using matrices.




