//
//  NeuralNetwork.cpp
//  NeuralNetwork
//
//  Created by Diamor Marke on 25/03/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "NeuralNetwork.hpp"

/*Activation functions and their derivatives
 x: the input into the function
 derivative: boolean to determine whether to return the derivative
 */

//fast sigmoid function; f(x) = x/(1+|x|)
double sigmoid(const double& x, const bool& derivative){
    if (derivative){return 1/((1 + abs(x))*(1 + abs(x)));}
    else {return x/(1 + abs(x));}
}

//ReLU; f(x) = max(x, 0)
double ReLU(const double& x, const bool& derivative){
    if (derivative){return (x >= 0) ? x : 0;}
    else {return (x >= 0) ? 1 : 0;}
}

//leaky ReLU; f(x) = x for x >= 0, 0.01 * x otherwise
double leakyReLU(const double& x, const bool& derivative){
    if (derivative){return (x >= 0) ? 1 : 0.01;}
    else {return (x >= 0) ? x : 0.01 * x;}
}

//linear: f(x) = x
double linear(const double& x, const bool& derivative){
    return (derivative ? 1 : x);
}

/*softmax: Defines the softmax function for output layer row-wise
 softmax(V) = exp(V)/sum(exp(V))
 softmax maps from R^n to [0, 1]^n with the constraints:
 sum(softmax(V)) = 1 and each component is greater than 0
 X: the matrix on which softmax is applied
 */
Matrix softmax(Matrix X){
    
    //maps exp to the input values
    X = X.map([](double x)->double{return exp(x);});
    
    //calculates the sum of each row
    Matrix sum = X * Matrix::one(X.col(), 1);
    
    //divides the value in each row by the sum of the row
    for (int i = 1; i <= X.row(); ++i){
        for (int j = 1; j <= X.col(); ++j){
            X(i,j) = X.val(i,j)/sum.val(i,1);
        }
    }
    
    return X;
}

/*NeuralNetwork; Initialisation of the Neural Network Object.
 inputs: The number of input nodes.
 outputs: The number of output nodes.
 hidden_layers: The numbers of nodes in hidden layers
 activation: The activation function being used.
 classify: Boolean determining whether the neural network is a classifier.
*/
NeuralNetwork::NeuralNetwork(const int& inputs, const int& outputs, deque<int> hidden_layers, function<double(double,bool)> activation, bool classify){
    
    //input and output layers must have a positive integer input
    if (inputs <= 0 or outputs <= 0){
        throw InvalidInputError();
    }
    
    //Copies initial parameters
    this->inputs = inputs;
    this->outputs = outputs;
    this->activation = activation;
    this->classify = classify;
    
    //Initialises hidden layer weights
    hidden_layers.push_back(outputs);
    hidden_layers.push_front(inputs);
    
    //Initialises hidden weights of layer i with values from Uniform(-1/sqrt(ni), 1/sqrt(ni)), where ni is the number of nodes in layer i;
    for (int i = 0; i < hidden_layers.size() - 1; ++i){
        W.push_back((2 * Matrix::rand(hidden_layers[i] + 1, hidden_layers[i+1]) - Matrix::one(hidden_layers[i] + 1, hidden_layers[i+1])) / (sqrt(hidden_layers[i] + 1)));
    }
    
    //Saves the total layer scheme
    this->layers = hidden_layers;
}

/*fit; Fits the Neural Network to the input data.
 X: The matrix of independent variables
 y: The matrix of dependent variables
 epochs: The number of epochs the data is trained for
 batch_size: The number of values trained before updating the weights
 learning_rate: The rate at which the weights are updated.
 */
void NeuralNetwork::fit(const Matrix& X, const Matrix& y, const int& epochs, const int& batch_size, const double& learning_rate){
    
    //Matrixes must be compatible with input parameters and each others
    if (inputs != X.col()  or outputs != y.col() or X.row() != y.row()){
        throw IncompatibleMatrixError();
    }
    
    //Initialises vector that will store gradients
    vector<Matrix> zeros;
    
    for (int i = 0; i < layers.size()-1; ++i){
        zeros.push_back(Matrix(layers[i]+1, layers[i+1]));
    }
    
    //Performs feedforward and backpropogation for a set number of epochs
    for (int epoch = 1; epoch <= epochs; ++epoch){
        
        //Sets the initial gradients to 0 of the given dimensions
        vector<Matrix> dW = zeros;
        
        //Sets initial parameters
        //Activated contains all the activated layers
        deque<Matrix> activated;
        
        //Preactivated contains all the preactivated layers
        deque<Matrix> preactivated;
        
        //Iterates over all the rows
        for (int row = 1; row <= X.row(); ++row){
            
            //Clears activated and preactivated layers
            activated.clear();
            preactivated.clear();
            
            //y_pred is the current prediction for y
            //it is also used to store the current value as it is being fed forward
            Matrix y_pred = X.getRow(row);
            activated.push_back(X.getRow(row));

            //Feeds the input value forward
            for (int i = 0; i < zeros.size() - 1; ++i){
                
                //Preactivated = Activated * Weights
                y_pred = (Matrix::concat({Matrix::one(1), y_pred},1) * W[i]);
                preactivated.push_back(y_pred);
                
                //Activated = Activation_function(Preactivated)
                y_pred = y_pred.map([this](double x){return activation(x, false);});
                activated.push_back(y_pred);
            }
            
            //Feeds forward to the final layer
            //Preactivated = Activated * Weights
            y_pred = (Matrix::concat({Matrix::one(1), y_pred},1) * W[zeros.size() - 1]);
            preactivated.push_back(y_pred);
            
            //Activated = Activation_function(Preactivated)
            if (classify){
                
                //Uses softmax if model is a classifier
                y_pred = softmax(y_pred);
            }
            else{
                
                //Uses activation function if model isn't a classifier
                y_pred = y_pred.map([this](double x){return activation(x, false);});
            }
            activated.push_back(y_pred);
            
            //Calculates the derivative of the loss function and starts backpropogation of the final layer
            Matrix loss;
            if (classify){
                loss = (y_pred - y.getRow(row))*(-1);
                
                //Calculates da/dz
                Matrix da = Matrix(y_pred.col());
                
                for (int i = 1; i <= da.row(); ++i){
                    for (int j = i+1; j <= da.row(); ++j){
                        da(i,j) = da(j,i) = -y_pred.val(1,i) * y_pred(1,j);
                    }
                }
                
                for (int i = 1; i <= da.row(); ++i){
                    da(i,i) = (1 - y_pred.val(1,i)) * y_pred(1,i);
                }
                
                //dL/dz = dL/da * da/dz
                loss = loss * da;
            }
            else{
                loss = y.getRow(row) - y_pred;
                
                //dL/dz = dL/da * da/dz
                loss = Matrix::dot(loss, preactivated[zeros.size() - 1].map([this](double x){return activation(x, true);}));
            }
            int i = zeros.size() - 1;
            
            //dL/dw = dL/dz * dz/dw
            dW[i] = dW[i] + (Matrix::concat({Matrix::one(1),activated[i]},1).transpose() * loss);
            
            //dL/da = dL/dz * dz/da
            loss = loss * W[i].transpose();
            
            //Backpropogates through the network generating gradients at all the remaining layers
            for (int i = zeros.size() - 2; i >= 0; --i){
                
                //dL/dz = dL/da * da/dz
                loss = Matrix::dot(loss.omitCol(1), preactivated[i].map([this](double x){return activation(x, true);}));
                
                //dL/dw = dL/dz * dz/dw
                dW[i] = dW[i] + (Matrix::concat({Matrix::one(1),activated[i]},1).transpose() * loss);
                
                //dL/da = dL/dz * dz/da
                loss = loss * W[i].transpose();
            }
            
            //Updates the weights after every set number of values
            if (row % batch_size == 0){
                for (int i = 0; i < layers.size()-1; ++i){
                    W[i] = W[i] + learning_rate * dW[i];
                }
                
                //Resets the gradient
                dW = zeros;
            }
        }
        
        //Generates an error report every number of epochs
        if ((epoch % 1) == 0){
            cout << "Epoch " << epoch << ": \t Error - " << (y - predict(X)).norm()/ y.row() / y.col() << endl;
        }
    }
}


/*predict; Generates a prediction for y given the current weights.
 X: The matrix of values used to predict y
 */
Matrix NeuralNetwork::predict(Matrix X){
    
    //X must have the correct number of columns
    if (inputs != X.col()){
        throw IncompatibleMatrixError();
    }
    
    
    //Feeds input value forward through network
    for (int i = 0; i < layers.size() - 2; ++i){
        X = (Matrix::concat({Matrix::one(X.row(), 1), X},1) * W[i]).map([this](double x){return activation(x, false);});
    }
    
    //Feeds forward through final layer
    if (classify){
        
        //Uses softmax if model is a classifier
        X = softmax(Matrix::concat({Matrix::one(X.row(), 1), X},1) * W[layers.size() - 2]);
    }else{
        
        //Uses activation function if not a classifier
        X = (Matrix::concat({Matrix::one(X.row(), 1), X},1) * W[layers.size() - 2]).map([this](double x){return activation(x, false);});
    }
    return X;
}
