//
//  linalg.hpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 05/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#ifndef linalg_hpp
#define linalg_hpp

#include "exceptions.h"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <math.h>
#include <random>
#include <stdio.h>
#include <vector>


using namespace std;

class Matrix{
private:
    int r;
    int c;
    vector<double> vals;
public:
    Matrix(const int& row, const int& col);                   //initialisers
    Matrix(const int& n) : Matrix(n, n){};
    Matrix(): Matrix(3, 3){};
    Matrix(const Matrix& M);
    Matrix(const vector<vector<double>>& arr);
    
    int row() const {return r;};                              //accessors
    int col() const {return c;};
    
    static Matrix identity(const int& n);                     //identity matrix
    static Matrix rotate_x(const double& theta);              //rotation in x-axis matrix
    static Matrix rotate_y(const double& theta);              //rotation in y-axis matrix
    static Matrix rotate_z(const double& theta);              //rotation in z-axis matrix
    
    double& operator()(const int& i, const int& j);           //non-constant accessor for values
    double val(const int& i, const int& j) const ;            //constant accessor for values
    
    friend ostream& operator<<(ostream& os, const Matrix& M); //output stream operator for matrices
    
    bool isSquare() const;                                    //matrix properties
    bool isOrthogonal() const;
    bool isSymmetric() const;
    bool isAntisymmetric() const;
    
    Matrix transpose() const;                                 //standard matrix functions
    Matrix inv() const;
    double det() const;
    double trace() const;
    double norm() const;
};

class Vector{
private:
    int n;
    vector<double> vals;
public:
    Vector(const int& n);                                      //initialisers
    Vector(): Vector(3){};
    Vector(const Vector& v);
    Vector(const vector<double>& v);
    Vector(const double arr[]);
    
    static Vector rand();
    
    int dim() const {return n;};
    double& operator()(const int& i);                          //non-constant accessor for values
    double val(const int& i) const ;                           //constant accessor for values
    
    friend ostream& operator<<(ostream& os, const Vector& v);  //output stream operator for vectors
    double norm() const;
};

Matrix operator+(const Matrix& A, const Matrix& B); //matrix addition
Matrix operator-(const Matrix& A, const Matrix& B); //matrix subtraction
Vector operator+(const Vector& A, const Vector& B); //vector addition
Vector operator-(const Vector& A, const Vector& B); //vector subtraction

Matrix operator*(const Matrix& A, const Matrix& B); //matrix multiplication
Matrix operator*(const double& k, const Matrix& M); //scalar matrix multiplication
Matrix operator*(const Matrix& M, const double& k);
Vector operator*(const double& k, const Vector& v); //scalar vector multiplication
Vector operator*(const Vector& v, const double& k);
double operator*(const Vector& u, const Vector& v); //vector dot product
Vector operator*(const Matrix& A, const Vector& v); //left matrix multiplication to column vector
Vector operator*(const Vector& v, const Matrix& A); //right matrix multiplication to row vector
Vector operator^(const Vector& u, const Vector& v); //vector product in R3

Matrix operator/(const Matrix& M, const double& k); //scalar matrix division
Vector operator/(const Vector& v, const double& k); //scalar vector division
Matrix operator/(const Matrix& A, const Matrix& B); //matrix division (multiplication of B by inverse of A)
bool operator==(const Matrix& A, const Matrix& B); //matrix equality operator
bool operator==(const Vector& u, const Vector& v); //vector equality operator
#endif /* linalg_hpp */
