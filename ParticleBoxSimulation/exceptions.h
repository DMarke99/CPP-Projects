//
//  exceptions.h
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 05/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef exceptions_h
#define exceptions_h
#include <exception>
using namespace std;

struct IncompatibleMatrixError : public exception {
    const char * what () const throw () {
        return "Matrix/Matrices incompatible for current operation";
    }
};

struct DimensionError : public exception {
    const char * what () const throw () {
        return "Dimensions of operand(s) are not compatible";
    }
};

struct NotInvertibleError : public exception {
    const char * what () const throw () {
        return "Matrix is not invertible";
    }
};

struct InvalidInputError : public exception {
    const char * what () const throw () {
        return "Input parameters are invalid";
    }
};

struct InconsistentInputError : public exception {
    const char * what () const throw () {
        return "File input is inconsistent; contains missing or extra columns in certain rows";
    }
};

#endif /* exceptions_h */
