//
//  matrix.cpp
//  ParticleBoxSimulation
//
//  Created by Diamor Marke on 05/06/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "linalg.hpp"

// Matrices
// Matrix Initialisers
Matrix::Matrix(const int& row, const int& col): r(row), c(col){
    if ((row<=0) or (col<=0)){
        throw invalid_argument("Invalid arguments");
    }
    
    vals = vector<double>(row*col,0);
}

Matrix::Matrix(const Matrix &M){
    this->r = M.row();
    this->c = M.col();
    this->vals = M.vals;
}

Matrix::Matrix(const vector<vector<double>>& arr){
    if (arr.size() == 0 or arr[1].size() == 0){
        throw invalid_argument("Invalid arguments");
    }
    int r = arr.size();
    int c = arr[1].size();
    for (vector<double> row : arr){
        if (row.size() != c){
            throw invalid_argument("Invalid arguments");
        }
    }
    
    this->r = r;
    this->c = c;
    
    for (int i = 0; i < r; ++i){
        for (int j = 0; j < c; ++j){
            this->vals.push_back(arr[i][j]);
        }
    }
}

//Special Matrices
Matrix Matrix::identity(const int& n){
    Matrix M = Matrix(n);
    for(int i = 1; i <= n; ++i){
        M(i,i) = 1;
    }
    return M;
}

Matrix Matrix::rotate_x(const double& theta){
    return Matrix({{1, 0, 0}, {0, cos(theta), -sin(theta)}, {0, sin(theta), cos(theta)}});
}

Matrix Matrix::rotate_y(const double& theta){
    return Matrix({{cos(theta), 0, -sin(theta)}, {0, 1, 0}, {sin(theta), 0, cos(theta)}});
}

Matrix Matrix::rotate_z(const double& theta){
    return Matrix({{cos(theta), -sin(theta), 0}, {sin(theta), cos(theta), 0}, {0, 0, 1}});
}

// Matrix Accessors
double& Matrix::operator()(const int& i, const int& j){
    if (i < 1 or i > row() or j < 1 or j > col()){
        throw InvalidInputError();
    }
    
    return vals[(i-1)*col()+(j-1)];
}

double Matrix::val(const int& i, const int& j) const {
    if (i < 1 or i > row() or j < 1 or j > col()){
        throw InvalidInputError();
    }
    
    return vals[(i-1)*col()+(j-1)];
};

// Matrix Booleans
bool Matrix::isSquare() const {
    return this->row() == this->col();
}

bool Matrix::isOrthogonal() const {
    Matrix A = this->transpose() * *this;
    return (A == identity(A.row()));
}

bool Matrix::isSymmetric() const {
    if (!isSquare()){
        return false;
    }
    
    for (int i = 1; i <= this->row(); ++i){
        for (int j = i+1; j <= this->row(); ++j){
            if (this->val(i,j) != this->val(j,i)){
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isAntisymmetric() const {
    if (!isSquare()){
        return false;
    }
    
    for (int i = 1; i <= this->row(); ++i){
        for (int j = i; j <= this->row(); ++j){
            if (this->val(i,j) != -this->val(j,i)){
                return false;
            }
        }
    }
    return true;
}

// Matrix Operations
Matrix Matrix::transpose() const {
    Matrix M = Matrix(this->col(), this->row());
    for(int i = 1; i <= row(); ++i){
        for(int j = 1; j <= col(); ++j){
            M(j,i) = this->val(i, j);
        }
    }
    return M;
}

Matrix Matrix::inv() const {
    if (!isSquare()){
        throw IncompatibleMatrixError();
    }
    
    int n = this->row();
    Matrix M = Matrix(*this);
    Matrix I = identity(n);
    
    //Finds inverse by performing Gaussian Elimination on the augmented matrix (M|I)
    for (int i = 1; i <= n; ++i){
        
        //Finds the first non-zero term in column i under row i
        int j = i;
        while (M(j,i) == 0){
            ++j;
            
            //If there is no non-zero value then the matrix is not invertible
            if (j>n){
                throw NotInvertibleError();
            }
        }
        
        //Swaps rows i and j if i != j
        if(j != i){
            double temp;
            
            //Swaps rows
            for (int k = 1; k <= this->row(); ++k){
                temp = M(j,k);
                M(j,k) = M(i,k);
                M(i,k) = temp;
                
                temp = I(j,k);
                I(j,k) = I(i,k);
                I(i,k) = temp;
            }
        }
        
        //Reduces all elements in row i by a factor of M(i,i)
        double factor = M(i,i);
        for (int k = 1; k <= n; ++k){
            M(i,k) = M(i,k)/factor;
            I(i,k) = I(i,k)/factor;
        }
        
        //Cancels terms all terms in column i under row i by using EROs
        for (int k = i+1; k <= n; ++k){
            factor = -M(k,i);
            for (int l = 1; l <= n; ++l){
                M(k,l) = M(k,l) + factor * M(i,l);
                I(k,l) = I(k,l) + factor * I(i,l);
            }
        }
        
        //Cancels terms all terms in row i to the right of column i by using ECOs
        for (int k = i+1; k <= n; ++k){
            factor = -M(i,k);
            for (int l = 1; l <= n; ++l){
                M(l,k) = M(l,k) + factor * M(l,i);
                I(l,k) = I(l,k) + factor * I(l,i);
            }
        }
        
    }

    return I;
}

double Matrix::det() const {
    if (!isSquare()){
        throw IncompatibleMatrixError();
    }
    
    double temp;
    double det = 1;
    int n = this->row();
    Matrix M = Matrix(*this);
    
    //Finds determinant by converting the matrix into an upper triangular matrix U
    //The determinant of U is the product of the terms on the leading diagonal
    //Multiplied by the sign of the permutation used to swap rows
    for (int i = 1; i <= n; ++i){
        
        //Finds the first non-zero term in column i under row i
        int j = i;
        while (M(j,i) == 0){
            ++j;
            
            //If there is no non-zero value then the determinant is 0
            if (j>n){
                return 0;
            }
        }
        
        //Swaps rows i and j if i != j
        if(j != i){
            
            //Accounts for sign of the permutation
            det = -det;
            
            //Swaps rows
            for (int k = 1; k <= this->row(); ++k){
                temp = M(j,k);
                M(j,k) = M(i,k);
                M(i,k) = temp;
            }
        }
        
        //Cancels terms all terms in column i under row i by using an ERO
        for (int k = i+1; k <= this->row(); ++k){
            double factor = -M(k,i)/M(i,i);
            for (int l = i; l <= n; ++l){
                M(k,l) = M(k,l) + factor * M(i,l);
            }
        }
    }
    
    //Calculates the product of the terms on the leading diagonal
    for (int i = 1; i <= n; ++i){
        det = det * M(i,i);
    }
    
    return det;
}

double Matrix::trace() const {
    if (!isSquare()){
        throw IncompatibleMatrixError();
    }
    
    double tr = 0;
    
    for (int i = 1; i <= this->row(); ++i){
        tr = tr + this->val(i,i);
    }
    
    return tr;
}

double Matrix::norm() const {
    double res = 0;
    for (int i = 1; i <= row(); ++i){
        for (int j = 1; j <= col(); ++j){
            res += val(i,j)*val(i,j);
        }
    }
    return sqrt(res);
}

ostream& operator<<(ostream& os, const Matrix& M){
    for(int i = 1; i <= M.row(); ++i){
        for(int j = 1; j <= M.col(); ++j){
            os << setw(5) << M.val(i,j);
        }
        os << endl;
    }
    return os;
}

// Vectors
// Vector Initialisers
Vector::Vector(const int& n){
    if (n<=0){
        throw invalid_argument("Invalid arguments");
    }
    
    this->n = n;
    vals = vector<double>(n,0);
}

Vector::Vector(const Vector& v){
    this->n = v.dim();
    this->vals = v.vals;
}

Vector::Vector(const vector<double>& v){
    if (v.size()<=0){
        throw invalid_argument("Invalid arguments");
    }
    
    this->n = (int) v.size();
    this->vals = v;
}

Vector::Vector(const double arr[]){
    if (sizeof(arr)<=0){
        throw invalid_argument("Invalid arguments");
    }

    this->vals = vector<double>(arr, arr+sizeof(arr)/sizeof(arr[0]));
    this->n = (int) this->vals.size();
    
}

//generates a spherically uniform random vector in R3
Vector Vector::rand(){
    random_device rand_dev;
    mt19937 generator(rand_dev());
    std::uniform_real_distribution<double>  distr(-1, 1);
    
    double u = distr(generator);
    double v = distr(generator);
    
    double theta = 2 * M_PI * v;
    return Vector({sqrt(1-u*u)*cos(theta),sqrt(1-u*u)*sin(theta),u});
}

// Vector Accessors
double& Vector::operator()(const int& i){
    if (i < 1 or i > dim()){
        throw InvalidInputError();
    }
    
    return vals[i-1];
}

double Vector::val(const int& i) const {
    if (i < 1 or i > dim()){
        throw InvalidInputError();
    }
    
    return vals[i-1];
};

ostream& operator<<(ostream& os, const Vector& v){
    for(int i = 1; i <= v.dim(); ++i){
        os << setw(5) << v.val(i);
    }
    os << endl;
    return os;
}

// Operators
Matrix operator+(const Matrix& A, const Matrix& B){
    if ((A.col() != B.col()) or (A.row() != B.row())){
        throw DimensionError();
    }
    
    Matrix M = Matrix(A.row(), A.col());
    
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            M(i,j) = A.val(i,j) + B.val(i,j);
        }
    }
    return M;
}

Matrix operator-(const Matrix& A, const Matrix& B){
    if ((A.col() != B.col()) or (A.row() != B.row())){
        throw DimensionError();
    }
    
    Matrix M = Matrix(A.row(), A.col());
    
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            M(i,j) = A.val(i,j) - B.val(i,j);
        }
    }
    return M;
}

Vector operator+(const Vector& u, const Vector& v){
    if (u.dim() != v.dim()){
        throw DimensionError();
    }
    
    Vector w = Vector(v.dim());
    
    for (int i = 1; i <= u.dim(); ++i){
        w(i) = u.val(i) + v.val(i);
    }
    return w;
}

Vector operator-(const Vector& u, const Vector& v){
    if (u.dim() != v.dim()){
        throw DimensionError();
    }
    
    Vector w = Vector(v.dim());
    
    for (int i = 1; i <= u.dim(); ++i){
        w(i) = u.val(i) - v.val(i);
    }
    return w;
}

Matrix operator*(const double& k, const Matrix& M){
    Matrix N = Matrix(M.row(), M.col());
    
    for (int i = 1; i <= M.row(); ++i){
        for (int j = 1; j <= M.col(); ++j){
            N(i,j) = k * M.val(i,j);
        }
    }
    return N;
}

Matrix operator*(const Matrix& A, const Matrix& B){
    if (A.col() != B.row()){
        throw DimensionError();
    }
    
    Matrix M = Matrix(A.row(), B.col());
    
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            for (int k = 1; k <= B.col(); ++k){
                M(i,k) = M(i,k) + A.val(i,j) * B.val(j,k);
            }
        }
    }
    return M;
}

Matrix operator*(const Matrix& M, const double& k){
    return k * M;
};

Vector operator*(const double& k, const Vector& v){
    Vector w = Vector(v.dim());
    
    for (int i = 1; i <= v.dim(); ++i){
        w(i) = k * v.val(i);
    }
    return w;
}

Vector operator*(const Vector& v, const double& k){
    return k * v;
}

double operator*(const Vector& u, const Vector& v){
    if (u.dim() != v.dim()){
        throw DimensionError();
    }
    
    double res = 0;
    
    for (int i = 1; i <= u.dim(); ++i){
        res += u.val(i) * v.val(i);
    }
    return res;
}

Vector operator*(const Matrix& A, const Vector& v){
    if (A.col() != v.dim()){
        throw DimensionError();
    }
    
    Vector b = Vector(A.row());
    
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
        b(i) += A.val(i,j) * v.val(j);
        }
    }
    return b;
}

Vector operator*(const Vector& v, const Matrix& A){
    if (A.row() != v.dim()){
        throw DimensionError();
    }
    
    Vector b = Vector(A.col());
    
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            b(i) += A.val(i,j) * v.val(i);
        }
    }
    return b;
}

Vector operator^(const Vector& u, const Vector& v){
    if (u.dim() != v.dim() or v.dim() != 3){
        throw DimensionError();
    }
    
    return Vector({u.val(2)*v.val(3)-u.val(3)*v.val(2),
                   u.val(3)*v.val(1)-u.val(1)*v.val(3),
                   u.val(1)*v.val(2)-u.val(2)*v.val(1)});
}

Vector operator/(const Vector& v, const double& k){
    return v * (1/k);
}

Matrix operator/(const Matrix& M, const double& k){
    return M * (1/k);
}

Matrix operator/(const Matrix& A, const Matrix& B){
    if (!B.isSquare() or !A.isSquare() or (B.row() != A.row())){
        throw IncompatibleMatrixError();
    }

    int n = B.row();
    Matrix M = B;
    Matrix I = A;
    
    //Divides by performing Gaussian Elimination on the augmented matrix (A|B)
    for (int i = 1; i <= n; ++i){
        
        //Finds the first non-zero term in column i under row i
        int j = i;
        while (M(j,i) == 0){
            ++j;
            
            //If there is no non-zero value then the matrix is not invertible
            if (j>n){
                throw NotInvertibleError();
            }
        }
        
        //Swaps rows i and j if i != j
        if(j != i){
            double temp;
            
            //Swaps rows
            for (int k = 1; k <= n; ++k){
                temp = M(j,k);
                M(j,k) = M(i,k);
                M(i,k) = temp;
                
                temp = I(j,k);
                I(j,k) = I(i,k);
                I(i,k) = temp;
            }
        }
        
        //Reduces all elements in row i by a factor of M(i,i)
        double factor = M(i,i);
        for (int k = 1; k <= n; ++k){
            M(i,k) = M(i,k)/factor;
            I(i,k) = I(i,k)/factor;
        }
        
        //Cancels terms all terms in column i under row i by using EROs
        for (int k = i+1; k <= n; ++k){
            factor = -M(k,i);
            for (int l = 1; l <= n; ++l){
                M(k,l) = M(k,l) + factor * M(i,l);
                I(k,l) = I(k,l) + factor * I(i,l);
            }
        }
        
        //Cancels terms all terms in row i to the right of column i by using ECOs
        for (int k = i+1; k <= n; ++k){
            factor = -M(i,k);
            for (int l = 1; l <= n; ++l){
                M(l,k) = M(l,k) + factor * M(l,i);
                I(l,k) = I(l,k) + factor * I(l,i);
            }
        }
        
    }
    
    //Returns the result of division
    return I;
}

bool operator==(const Matrix& A, const Matrix& B){
    if ((A.col() != B.col()) or (A.row() != B.row())){
        return false;
    }
    
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            if (A.val(i,j) == B.val(i,j)){
                return false;
            }
        }
    }
    
    return true;
}

bool operator==(const Vector& u, const Vector& v){
    if (u.dim() != v.dim()){
        return false;
    }
    
    for (int i = 1; i <= u.dim(); ++i){
        if (u.val(i) != v.val(i)){
            return false;
        }
    }
    
    return true;
}
