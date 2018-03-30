//
//  matrix.cpp
//  NeuralNetwork
//
//  Created by Diamor Marke on 23/03/2018.
//  Copyright © 2018 Diamor Marke. All rights reserved.
//

#include "matrix.hpp"

using namespace std;

//Initialisers

/*Standard initialiser;
 row: the number of rows in the matrix.
 col: the number of columns in the matrix.
*/
Matrix::Matrix(const int& row, const int& col): r(row), c(col){
    if ((row<0) or (col<0)){
        throw invalid_argument("Invalid arguments");
    }
    vals = vector<double>(row*col,0);
}

/*Copy constructor; Used for creating deep copies.
 M: the matrix being copied.
 */
Matrix::Matrix(const Matrix &M){
    this->r = M.row();
    this->c = M.col();
    this->vals = M.vals;
}

//Special matrices

/*Identity Matrix; Returns the Identity Matrix.
 n: the dimension of the square matrix.
 */
Matrix Matrix::identity(const int& n){
    Matrix M = Matrix(n);
    for(int i = 1; i <= n; ++i){
        M(i,i) = 1;
    }
    return M;
}

/*One Matrix; Returns a matrix of 1s.
 row: the number of rows in the matrix.
 col: the number of columns in the matrix.
 */
Matrix Matrix::one(const int& row, const int& col){
    Matrix M = Matrix(row, col);
    for(int i = 1; i <= row; ++i){
        for(int j = 1; j <= col; ++j){
            M(i,j) = 1;
        }
    }
    return M;
}

/*Random Matrix; Returns a matrix of random uniformly distributed numbers.
 row: the number of rows in the matrix.
 col: the number of columns in the matrix.
 */
Matrix Matrix::rand(const int& row, const int& col){
    default_random_engine generator;
    uniform_real_distribution<double> uniform(0,1);
    Matrix M = Matrix(row, col);
    for(int i = 1; i <= row; ++i){
        for(int j = 1; j <= col; ++j){
            M(i,j) = uniform(generator);
        }
    }
    return M;
}

/*Random Matrix; Returns a matrix of random normally distributed numbers.
 row: the number of rows in the matrix.
 col: the number of columns in the matrix.
 */
Matrix Matrix::randn(const int& row, const int& col){
    default_random_engine generator;
    normal_distribution<double> normal(0,1);
    Matrix M = Matrix(row, col);
    for(int i = 1; i <= row; ++i){
        for(int j = 1; j <= col; ++j){
            M(i,j) = normal(generator);
        }
    }
    return M;
}

//isSquare; Boolean to determine whether matrix is square.
bool Matrix::isSquare() const {
    return this->row() == this->col();
}

/*isOrthogonal; Boolean to determine whether matrix is orthogonal.
 Matrix M is orthogonal iff M^-1 = M.transpose.
 */
bool Matrix::isOrthogonal() const {
    Matrix A = this->transpose() * *this;
    return (A == identity(A.row()));
}

/*isSymmetric; Boolean to determine whether matrix is symmetric.
 Matrix M is orthogonal iff M(i,j) = M(j,i) for all i,j.
 */
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

/*isAntisymmetric; Boolean to determine whether matrix is symmetric.
 Matrix M is orthogonal iff M(i,j) = -M(j,i) for all i,j.
 */
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

//Accessors

/*operator(); Container-like accessor of matrix values.
 i: the row index of the value being accessed.
 j: the column index of the value being accessed.
 */
double& Matrix::operator()(const int& i, const int& j){
    if (i < 1 or i > row() or j < 1 or j > col()){
        throw InvalidInputError();
    }
    return vals[(i-1)*col()+(j-1)];
}

/*vals; Constant accessor of matrix values.
 i: the row index of the value being accessed.
 j: the column index of the value being accessed.
 */
double Matrix::val(const int& i, const int& j) const {
    if (i < 1 or i > row() or j < 1 or j > col()){
        throw InvalidInputError();
    }
    return vals[(i-1)*col()+(j-1)];
};

//transpose; Returns the transpose of a matrix.
Matrix Matrix::transpose() const {
    Matrix M = Matrix(this->col(), this->row());
    for(int i = 1; i <= row(); ++i){
        for(int j = 1; j <= col(); ++j){
            M(j,i) = this->val(i, j);
        }
    }
    return M;
}

/*operator>>; Output stream for matrices
 Enables C++-like displaying of matrix M via cout << M.
 */
ostream& operator<<(ostream& os, const Matrix& M){
    for(int i = 1; i <= M.row(); ++i){
        for(int j = 1; j <= M.col(); ++j){
            os << M.val(i,j) << "\t";
        }
        os << endl;
    }
    return os;
}

//Matrix Operations

/*operator+; Matrix Addition.
 A: The left matrix in addition.
 B: The right matrix in addition.
 */
Matrix operator+(const Matrix& A, const Matrix& B){
    
    //Matrices must be same dimensions
    if ((A.col() != B.col()) or (A.row() != B.row())){
        throw IncompatibleMatrixError();
    }
    
    Matrix M = Matrix(A.row(), A.col());
    
    //Performs addition elementwise
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            M(i,j) = A.val(i,j) + B.val(i,j);
        }
    }
    return M;
}

/*operator-; Matrix Subtraction.
 A: The left matrix in subtraction.
 B: The right matrix in subtraction.
 */
Matrix operator-(const Matrix& A, const Matrix& B){
    
    //Matrices must be same dimensions
    if ((A.col() != B.col()) or (A.row() != B.row())){
        throw IncompatibleMatrixError();
    }
    
    Matrix M = Matrix(A.row(), A.col());
    
    //Performs subtraction elementwise
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            M(i,j) = A.val(i,j) - B.val(i,j);
        }
    }
    return M;
}

/*operator*; Scalar Matrix Multiplication.
 k: The scalar in multiplication.
 M: The matrix in multiplication.
 */
Matrix operator*(const double& k, const Matrix& M){
    Matrix N = Matrix(M.row(), M.col());
    
    //Multiplies each element by k elementwise
    for (int i = 1; i <= M.row(); ++i){
        for (int j = 1; j <= M.col(); ++j){
            N(i,j) = k * M.val(i,j);
        }
    }
    return N;
}

/*operator*; Scalar Matrix Division.
 k: The scalar in division.
 M: The matrix in division.
 */
Matrix operator/(const Matrix& M, const double& k){
    return M * (1/k);
}

/*operator*; Matrix Multiplication.
 A: The left matrix in multiplication.
 B: The right matrix in multiplication.
 */
Matrix operator*(const Matrix& A, const Matrix& B){
    
    //Matrices are compatible when A.col = B.row
    if (A.col() != B.row()){
        throw IncompatibleMatrixError();
    }
    
    Matrix M = Matrix(A.row(), B.col());
    
    //Performs multiplication algorithm
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            for (int k = 1; k <= B.col(); ++k){
                M(i,k) = M(i,k) + A.val(i,j) * B.val(j,k);
            }
        }
    }
    return M;
}

/*operator/; Matrix Division.
 Finds the unique solution X to the equation XB = A for invertible B.
 A: The left matrix in division.
 B: The right matrix in division.
 */
Matrix operator/(const Matrix& A, const Matrix& B){
    
    //Method can only divide two square matrices of the same dimensions;
    if (!B.isSquare() or !A.isSquare() or (B.row() != A.row())){
        throw IncompatibleMatrixError();
    }
    
    //Generates initial parameters
    int n = B.row();
    
    //Copies matrices
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

//operator==; Defines equality for matrices.
bool operator==(const Matrix& A, const Matrix& B){
    
    //Matrices must have the same dimensions to be equal
    if ((A.col() != B.col()) or (A.row() != B.row())){
        return false;
    }
    
    //Matrices must have the same elements to be
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            if (A.val(i,j) == B.val(i,j)){
                return false;
            }
        }
    }
    
    return true;
}

//det; Returns the determinant of a matrix.
double Matrix::det() const {
    
    //Determinant is only defined for square matrices
    if (!isSquare()){
        throw IncompatibleMatrixError();
    }
    
    //Generates initial parameters
    double temp;
    double det = 1;
    int n = this->row();
    
    //Deep copies the current matrix
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

//inv; Returns the inverse of a matrix.
Matrix Matrix::inv() const {
    
    //Inverse is only defined for square matrices
    if (!isSquare()){
        throw IncompatibleMatrixError();
    }
    
    //Generates initial parameters
    int n = this->row();
    
    //Copies current matrix
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
    
    //Returns the inverse matrix
    return I;
}

//trace; Returns the trace of a matrix.
double Matrix::trace() const {
    
    //Trace is only defined for square matrices
    if (!isSquare()){
        throw IncompatibleMatrixError();
    }
    
    //generates initial parameters
    double tr = 0;

    //sums the terms on the diagonal
    for (int i = 1; i <= this->row(); ++i){
        tr = tr + this->val(i,i);
    }
    
    return tr;
}

/*norm; Returns the Frobenius norm of a matrix.
 M.norm() = sum(M(i,j)^2).
*/
double Matrix::norm() const {
    double res = 0;
    for (int i = 1; i <= row(); ++i){
        for (int j = 1; j <= col(); ++j){
            res += val(i,j)*val(i,j);
        }
    }
    return res;
}

/*map; Defines map for matrices.
 Applies a function elementwise to a matrix and returns the result.
 func: the function being applied elementwise.
*/
Matrix Matrix::map(const function<double(double)>& func) const {
    Matrix M = Matrix(*this);
    for (int i = 1; i <= M.row(); ++i){
        for (int j = 1; j <= M.col(); ++j){
            M(i,j) = func(M.val(i,j));
        }
    }
    return M;
}

/*getRow; Selects a particular row from the current matrix.
 i: the index of the row being selected.
 */
Matrix Matrix::getRow(const int& i) const {
    Matrix M = Matrix(1,this->col());
    for (int j = 1; j <= this->col(); ++j){
        M(1,j) = this->val(i,j);
    }
    return M;
}

/*getRow; Selects a set of rows from the current matrix.
 rows: A vector of indexes of rows being selected.
 */
Matrix Matrix::getRow(const vector<int>& rows) const {
    Matrix M = Matrix(rows.size(),this->col());
    for (int j = 1; j <= this->col(); ++j){
        for (int k = 1; k <= rows.size(); ++k){
            M(k,j) = this->val(rows[k-1],j);
        }
    }
    return M;
}

/*getCol; Selects a particular column from the current matrix.
 i: the index of the column being selected.
 */
Matrix Matrix::getCol(const int& i) const {
    Matrix M = Matrix(this->row(),1);
    for (int j = 1; j <= this->row(); ++j){
        M(j,1) = this->val(j,i);
    }
    return M;
}

/*getCol; Selects a set of columns from the current matrix.
 cols: A vector of indexes of column being selected.
 */
Matrix Matrix::getCol(const vector<int>& cols) const {
    Matrix M = Matrix(this->row(),cols.size());
    for (int j = 1; j <= this->row(); ++j){
        for (int k = 1; k <= cols.size(); ++k){
            M(j,k) = this->val(j, cols[k-1]);
        }
    }
    return M;
}

/*omitRow; Omits a particular row from the current matrix.
 i: the index of the row being omitted.
 */
Matrix Matrix::omitRow(const int& i) const {
    Matrix M = Matrix(this->row() - 1,this->col());
    for (int j = 1; j <= this->col(); ++j){
        for (int k = 1; k <= this->col() - 1; ++k){
            M(k,j) = this->val(k + ((k >= i) ? 1 : 0),j);
        }
    }
    return M;
}

/*omitCol; Omits a particular column from the current matrix.
 i: the index of the column being omitted.
 */
Matrix Matrix::omitCol(const int& i) const {
    Matrix M = Matrix(this->row(),this->col() - 1);
    for (int j = 1; j <= this->row(); ++j){
        for (int k = 1; k <= this->col() - 1; ++k){
            M(j,k) = this->val(j,k + ((k >= i) ? 1 : 0));
        }
    }
    return M;
}


/*concat; Concatenates a collection of matrices.
 Matrices: The vector of matrices being combined.
 axis: The axis on which the matrices are combined.
 0 if combined along the column and 1 if combined along the row.
 */
Matrix Matrix::concat(const vector<Matrix> &Matrices, const int &axis){
    
    //If vector is empty return the 0x0 matrix
    if (Matrices.size() == 0){
        return Matrix::identity(0);
    }
    
    //1 and 0 are the only valid values for axis
    else if(axis == 1 or axis == 0){
        
        //Gets the length of the axis perpendicular to the direction of concatenation
        int axis_len = axis == 0 ? Matrices[0].col() : Matrices[0].row();
        
        //calculates the concatenated length of the matrix
        int length = axis == 0 ? Matrices[0].row() : Matrices[0].col();
        for (int i = 1; i < Matrices.size(); ++i){
            
            //The lengths on the perpendicular axis must all be the same
            if ((axis == 0 and Matrices[i].col() != axis_len) or (axis == 1 and Matrices[i].row() != axis_len)){
                throw IncompatibleMatrixError();
            }
            
            //Adds the length of the current matrix to the running total
            length += axis == 0 ? Matrices[i].row() : Matrices[i].col();
        }
        
        //Generates the blank resultant matrix
        Matrix result = axis == 0 ? Matrix(length, axis_len) : Matrix(axis_len, length);
        
        //Index is the current index along the concatenated axis being copied
        int index = 0;
        
        //Iterates over all the matrices
        for (int i = 0; i < Matrices.size(); ++i){
            
            //Other is the index along the concatenated axis in the current matrix being copied
            int other = axis == 0 ? Matrices[i].row() : Matrices[i].col();
            
            //Iterates over values in current matrix
            for (int k = 1; k <= other; ++k){
                
                //Increments index
                ++index;
                for (int j = 1; j <= axis_len; ++j){
                    
                    //Copies values from current matrix into resultant matrix
                    if (axis == 0){
                        result(index,j) = Matrices[i].val(k,j);
                    }
                    else{
                        result(j,index) = Matrices[i].val(j,k);
                    }
                }
            }
        }
        
        //Returns resultant matrix
        return result;
    }
    //If axis isn't 1 or 0 throw an exception
    else{
        throw InvalidInputError();
    }
}

/*dot; Returns the elementwise product of two matrices
 A: The left matrix in the elementwise product.
 B: The right matrix in the elementwise product.
 */
Matrix Matrix::dot(const Matrix& A, const Matrix& B){
    
    //Matrices must have the same dimensions to be dotted
    if ((A.col() != B.col()) or (A.row() != B.row())){
        throw IncompatibleMatrixError();
    }
    
    Matrix M = Matrix(A.row(), A.col());
    
    //Calculates the elementwise product.
    for (int i = 1; i <= A.row(); ++i){
        for (int j = 1; j <= A.col(); ++j){
            M(i,j) = A.val(i,j) * B.val(i,j);
        }
    }
    return M;
}

/*import; Imports a csv file into a matrix.
 filepath: The filepath of the csv file being imported.
 delimiter: The delimiter for values in a particular row.
 header: A boolean determining whether there are headers to omit.
 */
Matrix import(const string& filepath, const char& delimiter, const bool& headers){
    ifstream ip(filepath);
    
    if(!ip.is_open()) cout << "ERROR" << endl;
    string input;
    vector<string> lines;
    
    while (ip.good()){
        getline(ip,input,'\n');
        lines.push_back(input);
    }
    ip.close();
    
    //Removes first row if there is a header
    if (headers){
        lines.erase(lines.begin());
    }
    
    
    if (lines.size() == 0){
        return Matrix::identity(0);
    }
    else{
        //Determines matrix width
        stringstream ss(lines[0]);
        vector<string> tokens;
        
        while (getline(ss,input,',')){
            tokens.push_back(input);
        }
        
        int cols = tokens.size();
        
        //Creates resultant matrix
        Matrix res = Matrix(lines.size(), cols);
        
        //Iterates over every row
        for (int i = 1; i <= lines.size(); ++i){
            
            //Splits lines by using stringstreams
            stringstream ss(lines[i-1]);
            vector<string> line_tokens;
            
            //Splits line and stores value in a matrix
            while (getline(ss,input,delimiter)){
                line_tokens.push_back(input);

            }
            
            //All columns must have the same number of elements
            if (tokens.size() != cols){
                throw InconsistentInputError();
            }
            
            //the predicate for removal of characters in strings
            auto predicate = [](const char& c){
                string check = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890.";
                return (check.find(c) == string::npos);
            };
            
            //Converts values to doubles
            for (int j = 1; j <= cols; ++j){
                
                //Cleans up string before conversion
                line_tokens[j-1].erase(std::remove_if(line_tokens[j-1].begin(), line_tokens[j-1].end(), predicate), line_tokens[j-1].end());
                
                //Converts string into double
                stringstream(line_tokens[j-1]) >> res(i,j);
            }
        }

        return res;
    }    
}
