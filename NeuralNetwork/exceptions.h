//
//  exceptions.h
//  NeuralNetwork
//
//  Created by Diamor Marke on 23/03/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef exceptions_h
#define exceptions_h
#include <exception>
using namespace std;

//IncompatibleMatrixError; Called whenever matrix/matrices are incompatible for a particular operation.
struct IncompatibleMatrixError : public exception {
    const char * what () const throw () {
        return "Matrix/Matrices incompatible for current operation";
    }
};

//NotInvertibleError; Called whenever matrix/matrices operations attempt to invert a non-invertible matrix
struct NotInvertibleError : public exception {
    const char * what () const throw () {
        return "Matrix is not invertible";
    }
};

//InvalidInputError; Called whenever inputs parameters don't make sense in the context of the method
struct InvalidInputError : public exception {
    const char * what () const throw () {
        return "Input parameters are invalid";
    }
};

//InconsistentInputError; Called when erronous data is imported
struct InconsistentInputError : public exception {
    const char * what () const throw () {
        return "File input is inconsistent; contains missing or extra columns in certain rows";
    }
};

#endif /* exceptions_h */
