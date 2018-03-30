//
//  main.cpp
//  NeuralNetwork
//
//  Created by Diamor Marke on 23/03/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "NeuralNetwork.cpp"

int main(int argc, const char * argv[]) {
    
    /*Creates example X and y to run on neural net
     X is constant with random error
     y is a softmaxed constant with random error
     Ideally the neural net will detect the constant component of X and use it to predict the constant component of y
     */
    Matrix X = Matrix::randn(1000,8) + 0.5 * Matrix::one(1000, 8);
    Matrix y = 0.05 * Matrix::randn(1000,3) + 0.6 * Matrix::one(1000, 3);
    y = softmax(y);
    
    NeuralNetwork Regressor = NeuralNetwork(8, 3, {8, 5}, leakyReLU, true);
    Regressor.fit(X,y,100);
    
    return 0;
}
