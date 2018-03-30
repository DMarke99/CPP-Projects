//
//  matrix.hpp
//  NeuralNetwork
//
//  Created by Diamor Marke on 23/03/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <fstream>
#include "exceptions.h"
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;
//Matrix class that encapsulates a mathematical Matrix
class Matrix{
private:
    
    //Attributes
    int r;
    int c;
    vector<double> vals;
public:
    
    //Initialisers
    Matrix(const int& row, const int& col);
    Matrix(const int& n) : Matrix(n, n){};
    Matrix(): Matrix(3, 3){};
    Matrix(const Matrix& M);
    
    //Destructor
    ~Matrix(){};
    
    //Accessors for matrix dimensions
    int row() const {return r;};
    int col() const {return c;};
    
    //Special matrices
    static Matrix identity(const int& n);
    static Matrix one(const int& n){return one(n,n);};
    static Matrix one(const int& row, const int& col);
    static Matrix rand(const int& n){return rand(n,n);};
    static Matrix rand(const int& row, const int& col);
    static Matrix randn(const int& n){return randn(n,n);};
    static Matrix randn(const int& row, const int& col);
    
    //Accessors for values in matrix
    double& operator()(const int& i, const int& j); //used for accessing reference to value in array 
    double val(const int& i, const int& j) const ; //used for accessing value in array
    
    //Display of matrix
    friend ostream& operator<<(ostream& os, const Matrix& M);
    
    //Matrix properties
    bool isSquare() const;
    bool isOrthogonal() const;
    bool isSymmetric() const;
    bool isAntisymmetric() const;
    
    //Matrix functions
    Matrix transpose() const;
    Matrix inv() const;
    double det() const;
    double trace() const;
    double norm() const;
    Matrix map(const function<double(double)>& func) const;
    
    //Selection operations
    Matrix getRow(const int& i) const;
    Matrix getRow(const vector<int>& i) const;
    Matrix getCol(const int& i) const;
    Matrix getCol(const vector<int>& i) const;
    Matrix omitRow(const int& i) const;
    Matrix omitCol(const int& i) const;
    
    //Special matrix operators
    static Matrix concat(const vector<Matrix>& Matrices, const int& axis=0);
    static Matrix dot(const Matrix& A, const Matrix& B);
};

//Matrix Operations
Matrix operator+(const Matrix& A, const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator*(const double& k, const Matrix& M);
Matrix operator*(const Matrix& M, const double& k){return k * M;};
Matrix operator/(const Matrix& M, const double& k);
Matrix operator/(const Matrix& A, const Matrix& B);
bool operator==(const Matrix& A, const Matrix& B);

//Import; imports a matrix from a csv file at a specified filepath
Matrix import(const string& filepath, const char& delimiter, const bool& headers);

#endif /* matrix_hpp */
