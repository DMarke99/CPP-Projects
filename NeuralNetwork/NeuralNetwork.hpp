//
//  NeuralNetwork.hpp
//  NeuralNetwork
//
//  Created by Diamor Marke on 25/03/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef NeuralNetwork_hpp
#define NeuralNetwork_hpp

#include <deque>
#include "matrix.cpp"

//Activation functions
double sigmoid(const double& x, const bool& derivative);
double ReLU(const double& x, const bool& derivative);
double leakyReLU(const double& x, const bool& derivative);
double linear(const double& x, const bool& derivative);

//Output layer activation for classification models/probability
Matrix softmax(Matrix X);

//NeuralNetwork object that encapsulates a Neural Network
class NeuralNetwork{
private:
    
    //Neural Network Attributes
    int inputs;
    int outputs;
    vector<Matrix> W;
    deque<int> layers;
    function<double(double,bool)> activation;
    bool classify;
    
public:
    
    //Initialisation of the neural network object
    NeuralNetwork(const int& input, const int& output, deque<int> hidden_layers={}, function<double(double,bool)> activation=leakyReLU, bool classify=true);
    
    //Fits the neural network to input data X and output data y
    void fit(const Matrix& X, const Matrix& y, const int& epochs=100, const int& batch_size=32, const double& learning_rate=0.0001);
    
    //Using the current weights, generates a prediction for y given input X
    Matrix predict(Matrix X);
};
#endif /* NeuralNetwork_hpp */
