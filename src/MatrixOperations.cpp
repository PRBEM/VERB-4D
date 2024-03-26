#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include "MatrixOperations.h"

extern "C" {

int dgemv_(char *trans, long *m, long *n, double *alpha, 
            double *a, long *lda, double *x, long *incx, 
            double *beta, double *y, long *incy);

int dgemm_(char *transa, char *transb, long *m, long *n,
            long *k, double *alpha, double *a, long *lda, 
            double *b, long *ldb, double *beta, double *c__, 
            long *ldc);

int dgetrf_(long *m, long *n, double *a, long * lda, 
            long *ipiv, long *info);

int dgetri_(long *n, double *a, long *lda, long *ipiv, 
            double *work, long *lwork, long *info);

int dgesv_(long *n, long* nrhs, double* a, long* lda,
          long* ipiv, double* b, long* ldb, long* info);
}

Matrix1D<double> cat(const Matrix1D<double>& a, const Matrix1D<double>& b)
{
    if (!a.initialized && !b.initialized) {
        return a;
    } else if (!a.initialized && b.initialized) {
        return b;
    } else if (a.initialized && !b.initialized) {
        return a;
    } else {
        int size_total = a.size_q1 + b.size_q1;

        Matrix1D<double> result (size_total);

        for (int i = 0; i < a.size_q1; ++i)
            result[i] = a[i];

        for (int i = a.size_q1; i < size_total; ++i)
            result[i] = b[i - a.size_q1];

        return result;
    }
}

Matrix2D<double> cat(const Matrix2D<double>& a, const Matrix2D<double>& b)
{
    if (!a.initialized && !b.initialized) {
        return a;
    } else if (!a.initialized && b.initialized) {
        return b;
    } else if (a.initialized && !b.initialized) {
        return a;
    } else {
        if (a.size_q2 != b.size_q2) {
            std::cout << "Concatentated matrices are of different size\n";
            exit(EXIT_FAILURE);
        }

        int size_q2 = a.size_q2;

        int size_total = a.size_q1 + b.size_q1;

        Matrix2D<double> result (size_total, size_q2);

        for (int i = 0; i < a.size_q1; ++i)
            for (int j = 0; j < size_q2; ++j)
                    result[i][j] = a[i][j];

        for (int i = a.size_q1; i < size_total; ++i)
            for (int j = 0; j < size_q2; ++j)
                    result[i][j] = b[i - a.size_q1][j];

        return result;
    }
}

Matrix3D<double> cat(const Matrix3D<double>& a, const Matrix3D<double>& b)
{
    if (!a.initialized && !b.initialized) {
        return a;
    } else if (!a.initialized && b.initialized) {
        return b;
    } else if (a.initialized && !b.initialized) {
        return a;
    } else {
        if ((a.size_q2 != b.size_q2) || (a.size_q3 != b.size_q3)) {
            std::cout << "Concatentated matrices are of different size\n";
            exit(EXIT_FAILURE);
        }

        int size_q2 = a.size_q2;
        int size_q3 = a.size_q3;

        int size_total = a.size_q1 + b.size_q1;

        Matrix3D<double> result (size_total, size_q2, size_q3);

        for (int i = 0; i < a.size_q1; ++i)
            for (int j = 0; j < size_q2; ++j)
                for (int k = 0; k < size_q3; ++k)
                    result[i][j][k] = a[i][j][k];

        for (int i = a.size_q1; i < size_total; ++i)
            for (int j = 0; j < size_q2; ++j)
                for (int k = 0; k < size_q3; ++k)
                    result[i][j][k] = b[i - a.size_q1][j][k];

        return result;
    }
}

Matrix2D<double> transpose(const Matrix2D<double>& A) {
    Matrix2D<double> result(A.size_q2, A.size_q1);
    for (auto i = 0; i < result.size_q1; ++i) {
        for (auto j = 0; j < result.size_q2; ++j) {
            result[i][j] = A[j][i];
        }
    }
    return result;
}

Matrix1D<double> log10(const Matrix1D<double>& A) {
   Matrix1D<double> result {A.size_q1};
   for (auto i = 0; i < result.size_q1; ++i) {
       if (A[i] < 1e-31) result[i] = -31;
       else result[i] = std::log10(A[i]);
   }
   return result;
}

Matrix2D<double> log10(const Matrix2D<double>& A) {
   Matrix2D<double> result {A.size_q1, A.size_q2};
    for (auto i = 0; i < result.size_q1; ++i) {
        for (auto j = 0; j < result.size_q2; ++j) {
            if (A[i][j] < 1e-31) result[i][j] = -31;
            else result[i][j] = std::log10(A[i][j]);
        }
    }
    return result;
}

Matrix1D<double> pow(double base, const Matrix1D<double>& A) {
    Matrix1D<double> result {A.size_q1};
    for (auto i = 0; i < result.size_q1; ++i) {
        result[i] = std::pow(base, A[i]); 
    }
    return result;
}

Matrix2D<double> diag(const Matrix1D<double>& d) {
    Matrix2D<double> result {d.size_q1, d.size_q1};   
    for (auto i = 0; i < result.size_q1; ++i) {
        for (auto j = 0; j < result.size_q2; ++j) {
            result[i][j] = (i == j) ? d[i] : 0.;
        }
    }
    return result;
}


template<typename Comparator>
Matrix1D<bool> compare(
    const Matrix1D<double>& A, double number, 
    Comparator comparator
) {        
    Matrix1D<bool> result {A.size_q1};
    for(auto i = 0; i < A.size_q1; ++i) {
        result[i] = comparator(A[i], number);
    }
    return result;
}

Matrix1D<bool> operator==(const Matrix1D<double>& A, double number) {        
    return compare(A, number, std::equal_to<double>{});
}

Matrix1D<bool> operator>(const Matrix1D<double>& A, double number) {        
    return compare(A, number, std::greater<double>{});
}

Matrix1D<bool> operator>=(const Matrix1D<double>& A, double number) {        
    return compare(A, number, std::greater_equal<double>{});
}

Matrix1D<bool> operator<(const Matrix1D<double>& A, double number) {        
    return compare(A, number, std::less<double>{});
}

Matrix1D<bool> operator<=(const Matrix1D<double>& A, double number) {        
    return compare(A, number, std::less_equal<double>{});
}

Matrix1D<bool> operator==(double number, const Matrix1D<double>& A) {        
    return A == number;
}

Matrix1D<bool> operator>(double number, const Matrix1D<double>& A) {        
    return A < number;
}

Matrix1D<bool> operator>=(double number, const Matrix1D<double>& A) {        
    return A <= number;
}

Matrix1D<bool> operator<(double number, const Matrix1D<double>& A) {        
    return A > number;
}

Matrix1D<bool> operator<=(double number, const Matrix1D<double>& A) {        
    return A >= number;
}

template<typename Operation>
Matrix1D<bool> logicalOperation(
    const Matrix1D<bool>& a, const Matrix1D<bool>& b,
    Operation operation
) {
    if (a.size_q1 != b.size_q1) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "a.size_q1 = " << a.size_q1 << ", ";
        std::cout << "b.size_q1 = " << b.size_q1 << std::endl;; 
        exit(EXIT_FAILURE);
    }
    
    Matrix1D<bool> result {a.size_q1};

    for (auto i = 0; i < a.size_q1; ++i) {
        result[i] = operation(a[i], b[i]);
    }

    return result;
}

Matrix1D<bool> operator&&(const Matrix1D<bool>& a, const Matrix1D<bool>& b) {
    return logicalOperation(a, b, std::logical_and<bool>{});
}

Matrix1D<bool> operator||(const Matrix1D<bool>& a, const Matrix1D<bool>& b) {
    return logicalOperation(a, b, std::logical_or<bool>{});
}

/* Matrix1D<bool> operator&&(const Matrix1D<bool>& a, const Matrix1D<bool>& b) {
    if (a.size_q1 != b.size_q1) {
        cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        cout << "Matrix sizes are different. ";
        cout << "a.size_q1 = " << a.size_q1 << ", ";
        cout << "b.size_q1 = " << b.size_q1 << endl;; 
        exit(EXIT_FAILURE);
    }
    
    Matrix1D<bool> result {a.size_q1};

    for (auto i = 0; i < a.size_q1; ++i) {
        result[i] = (a[i] && b[i]);
    }

    return result;
} */

template<typename Comparator>
Matrix2D<bool> compare(
    const Matrix2D<double>& A, double number, Comparator comparator
) {        
    Matrix2D<bool> result {A.size_q1, A.size_q2};
    for(auto i = 0; i < A.size_q1; ++i) {
        for (auto j = 0; j < A.size_q2; ++j) {
            result[i][j] = comparator(A[i][j], number);
        }
    }
    return result;
}

Matrix2D<bool> operator==(const Matrix2D<double>& A, double number) {        
    return compare(A, number, std::equal_to<double>{});
}

Matrix2D<bool> operator>(const Matrix2D<double>& A, double number) {        
    return compare(A, number, std::greater<double>{});
}

Matrix2D<bool> operator>=(const Matrix2D<double>& A, double number) {        
    return compare(A, number, std::greater_equal<double>{});
}

Matrix2D<bool> operator<(const Matrix2D<double>& A, double number) {        
    return compare(A, number, std::less<double>{});
}

Matrix2D<bool> operator<=(const Matrix2D<double>& A, double number) {        
    return compare(A, number, std::less_equal<double>{});
}

Matrix2D<bool> operator==(double number, const Matrix2D<double>& A) {        
    return A == number;
}

Matrix2D<bool> operator>(double number, const Matrix2D<double>& A) {        
    return A < number;
}

Matrix2D<bool> operator>=(double number, const Matrix2D<double>& A) {        
    return A <= number;
}

Matrix2D<bool> operator<(double number, const Matrix2D<double>& A) {        
    return A > number;
}

Matrix2D<bool> operator<=(double number, const Matrix2D<double>& A) {        
    return A >= number;
}

Matrix2D<bool> operator&&(const Matrix2D<bool>& A, const Matrix2D<bool>& B) {
    if (A.size_q1 != B.size_q1 || A.size_q2 != B.size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "A.size_q1 = " << A.size_q1 << ", ";
        std::cout << "B.size_q1 = " << B.size_q1 << ", "; 
        std::cout << "A.size_q2 = " << A.size_q2 << ", ";
        std::cout << "B.size_q2 = " << B.size_q2 << std::endl; 
        exit(EXIT_FAILURE);
    }
    
    Matrix2D<bool> result {A.size_q1, A.size_q2};

    for (auto i = 0; i < A.size_q1; ++i) {
        for (auto j = 0; j < A.size_q2; ++j) {                
            result[i][j] = (A[i][j] && B[i][j]);
        }
    }

    return result;
}

Matrix1D<double> operator-(
    const Matrix1D<double>& a, const Matrix1D<double>& b
) {
    if (a.size_q1 != b.size_q1) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "a.size_q1 = " << a.size_q1 << ", ";
        std::cout << "b.size_q1 = " << b.size_q1 << std::endl;
        exit(EXIT_FAILURE);
    }

    Matrix1D<double> result {a.size_q1};
    for (auto i = 0; i < result.size_q1; ++i) {
        result[i] = a[i] - b[i];
    } 

    return result;
} 

Matrix1D<double> operator*(
    const Matrix2D<double>& A, const Matrix1D<double>& b
) {
    if (A.size_q2 != b.size_q1) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix/vector size mismatch. ";
        std::cout << "A.size_q2 = " << A.size_q2 << ", ";
        std::cout << "b.size_q1 = " << b.size_q1 << std::endl; 
        exit(EXIT_FAILURE);
    }

    Matrix1D<double> answer {A.size_q1};

    //transpose matrix A as lapack expects column major matrices 
    char trans {'T'};
    long m = static_cast<long>(A.size_q2);
    long n = static_cast<long>(A.size_q1);
    double alpha {1.};
    double *a {A[0]};
    long lda {m};
    double *x {b.matrix_array};
    long incx {1};
    double beta {0.};
    double *y {answer.matrix_array};
    long incy {1};

/*     int dgemv_(char *trans, long *m, long *n, double *alpha, 
                double *a, long *lda, double *x, long *incx, 
                double *beta, double *y, long *incy); */
    dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

    return answer;
} 

Matrix2D<double> operator*(
    const Matrix2D<double>& A, const Matrix2D<double>& B
) {

    if (A.size_q2 != B.size_q1) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix/matrix size mismatch. ";
        std::cout << "A.size_q2 = " << A.size_q2 << ", ";
        std::cout << "B.size_q1 = " << B.size_q1 << std::endl; 
        exit(EXIT_FAILURE);
    }
    //Lapack uses column major matrices (while we use row major) and therefore
    //returns the transpose of the product of A and B
    //to avoid an additional transpose in the end we compute B^T * A^T instead
    Matrix2D<double> answer {A.size_q1, B.size_q2};

    //transpose matrices A and B as lapack expects column major matrices 
    char transa {'N'};
    char transb {'N'};
    long m = static_cast<long>(B.size_q2);
    long n = static_cast<long>(A.size_q1);
    long k = static_cast<long>(B.size_q1);
    double alpha {1.0};
    double *a {A[0]};
    long lda = static_cast<long>(A.size_q2);
    double *b {B[0]};
    long ldb = static_cast<long>(B.size_q2);
    double beta {0.};
    double *c {answer[0]};
    long ldc {m};
    
    dgemm_(&transb, &transa, &m, &n, &k, &alpha, b, &ldb, a, &lda, &beta, 
            c, &ldc);

    return answer;
}


Matrix2D<double> operator^(Matrix2D<double> A, const int power) {
    if (A.size_q1 != A.size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix/matrix size mismatch. ";
        std::cout << "A.size_q1 = " << A.size_q2 << ", ";
        std::cout << "A.size_q2 = " << A.size_q2 << std::endl;
        exit(EXIT_FAILURE);
    }

    if (power == 1) {
        // nothing to do
        return A;
    }

    // we precalculate powers of 2 for efficiency (A^2, A^4, A^16, ...)
    std::bitset<32> power_bitset(power);
    //std::cout << "power: " << power << std::endl;
    //std::cout << "bitset: " << power_bitset << std::endl;
    //std::cout << "bitset count: " << power_bitset.count() << std::endl;
    int cache_size = power_bitset[0] ? power_bitset.count()-1 : power_bitset.count();

    std::vector<Matrix2D<double>> cache_matrixes(cache_size);
    for (int i = 0; i < cache_matrixes.size(); i++) {
        cache_matrixes[i] = Matrix2D<double>(A.size_q1, A.size_q2);
    }

    char transa{'N'};
    char transb{'N'};
    long m{A.size_q2};
    long n{A.size_q1};
    long k{A.size_q1};
    double alpha{1.0};
    long lda{A.size_q2};
    long ldb{A.size_q2};
    double beta{0.};
    long ldc{m};

    double* a;
    double* b;
    double* c;

    Matrix2D<double> A_copy = A; // A acts as a temporary
    int cache_matrix_cursor = 0;

    for (int i = 1; i < power_bitset.size(); i++) {
        //std::cout << "Iteration: " << i << std::endl;

        if (i == 1) {
            a = A[0];
            c = cache_matrixes[cache_matrix_cursor][0];
        } else {
            a = cache_matrixes[cache_matrix_cursor-1][0];
            c = cache_matrixes[cache_matrix_cursor][0];
        }

        dgemm_(&transb, &transa, &m, &n, &k, &alpha, a, &ldb, a, &lda, &beta,
            c, &ldc);

        // preparation for while loop
        if (power_bitset[i] == false and i > 1) { // use A as temporary; is already done for i == 1
            //std::cout << "A as temporary" << std::endl;
            A = cache_matrixes[cache_matrix_cursor];

            c = A[0];
            a = cache_matrixes[cache_matrix_cursor][0];
        }

        while(power_bitset[i] == false) {
            //std::cout << "bit false: " << i << std::endl;

            // switch a and c
            double* tmp = a;
            a = c;
            c = tmp;

            dgemm_(&transb, &transa, &m, &n, &k, &alpha, a, &ldb, a, &lda, &beta,
                c, &ldc);

            i++;
        }

        // if result ended up being in A, copy it back to the cache matrix
        if (c != cache_matrixes[cache_matrix_cursor][0]) {
            //std::cout << "Result in A, copying" << std::endl;
            cache_matrixes[cache_matrix_cursor] = A;
        }

        cache_matrix_cursor++;
        if (cache_matrix_cursor == cache_matrixes.size()) {
            break; // all filled
        }

    }

    Matrix2D<double> answer{A.size_q1, A.size_q2};

    if (cache_matrixes.size() == 1) {
        answer = cache_matrixes[0];
        c = answer[0]; // mark that answer is holding the result
    } else {
        for (int i = 1; i < cache_matrixes.size(); i++) {
            //std::cout << "Iteration: " << i << std::endl;

            if (i == 1) {
                a = cache_matrixes[i-1][0];
                b = cache_matrixes[i][0];
                c = answer[0];
            } else if (i % 2 == 0) {
                a = answer[0];
                b = cache_matrixes[i][0];
                c = A[0]; // we can safely reuse A here; it's data is not important
            } else {
                a = A[0];
                b = cache_matrixes[i][0];
                c = answer[0];
            }

            dgemm_(&transb, &transa, &m, &n, &k, &alpha, a, &ldb, b, &lda, &beta,
                c, &ldc);

        }
    }

    if (c == A[0]) { // result in A
        //std::cout << "result in A" << std::endl;
        if (power % 2 == 1) {
            a = A[0];
            b = A_copy[0];
            c = answer[0];

            dgemm_(&transb, &transa, &m, &n, &k, &alpha, a, &ldb, b, &lda, &beta,
            c, &ldc);

            return answer;
        }

        return A;

    } else { // result in answer
        //std::cout << "result in answer" << std::endl;
        if (power % 2 == 1) {
            a = answer[0];
            b = A_copy[0];
            c = A[0];

            dgemm_(&transb, &transa, &m, &n, &k, &alpha, a, &ldb, b, &lda, &beta,
            c, &ldc);

            return A;
        }

        return answer;

    }
}


Matrix2D<double> abtrans(
    const Matrix2D<double>& A, const Matrix2D<double>& B
) {

    if (A.size_q2 != B.size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix/matrix size mismatch. ";
        std::cout << "A.size_q2 = " << A.size_q2 << ", ";
        std::cout << "B.size_q2 = " << B.size_q2 << std::endl; 
        exit(EXIT_FAILURE);
    }

    Matrix2D<double> answer {A.size_q1, B.size_q1};

    //transpose matrices A and B as lapack expects column major matrices 
    char transa {'N'};
    char transb {'T'};
    long m = static_cast<long>(B.size_q1);
    long n = static_cast<long>(A.size_q1);
    long k = static_cast<long>(B.size_q2);
    double alpha {1.0};
    double *a {A[0]};
    long lda {k};
    double *b {B[0]};
    long ldb = static_cast<long>(B.size_q2);
    double beta {0.};
    double *c {answer[0]};
    long ldc {m};

    dgemm_(&transb, &transa, &m, &n, &k, &alpha, b, &ldb, a, &lda, &beta, 
            c, &ldc);

    return answer;
} 

Matrix2D<double> operator+(
    const Matrix2D<double>& A, const Matrix2D<double>& B) {

    if (A.size_q1 != B.size_q1 || A.size_q2 != B.size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "A.size_q1 = " << A.size_q1 << ", ";
        std::cout << "B.size_q1 = " << B.size_q1 << ", "; 
        std::cout << "A.size_q2 = " << A.size_q2 << ", ";
        std::cout << "B.size_q2 = " << B.size_q2 << std::endl; 
        exit(EXIT_FAILURE);
    }

    Matrix2D<double> result {A.size_q1, A.size_q2};
    for (size_t i = 0; i < result.size_q1; ++i) {
        for (size_t j = 0; j < result.size_q2; ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
} 

Matrix2D<double> operator-(
    const Matrix2D<double>& A, const Matrix2D<double>& B) {

    if (A.size_q1 != B.size_q1 || A.size_q2 != B.size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix sizes are different. ";
        std::cout << "A.size_q1 = " << A.size_q1 << ", ";
        std::cout << "B.size_q1 = " << B.size_q1 << ", "; 
        std::cout << "A.size_q2 = " << A.size_q2 << ", ";
        std::cout << "B.size_q2 = " << B.size_q2 << std::endl; 
        exit(EXIT_FAILURE);
    }

    Matrix2D<double> result {A.size_q1, A.size_q2};

    for (size_t i = 0; i < result.size_q1; ++i) {
        for (size_t j = 0; j < result.size_q2; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
} 


Matrix2D<double> inv(Matrix2D<double> A) {

    if (A.size_q1 != A.size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Matrix is not a square one. ";
        std::cout << "A.size_q1 = " << A.size_q1 << ", ";
        std::cout << "A.size_q2 = " << A.size_q2 << std::endl;
        exit(EXIT_FAILURE);
    }


    long m = static_cast<long>(A.size_q1);
    long n = static_cast<long>(A.size_q1);   
    double *a {A[0]};
    long lda {m};
    std::unique_ptr<long[]> ipiv {new long[m]};
    long info;
    long lwork {n};
    std::unique_ptr<double[]> work {new double[lwork]};

/*     int dgetrf_(long *m, long *n, double *a, long * lda, 
                long *ipiv, long *info);

    int dgetri_(long *n, double *a, long *lda, long *ipiv, 
                double *work, long *lwork, long *info); */
    dgetrf_(&m, &n, a, &lda, ipiv.get(), &info);

    dgetri_(&n, a, &lda, ipiv.get(), work.get(), &lwork, &info);

    return A;
}
// solve the system of transposed matrices A^T * X^T = B^T equivalent to X * A = B
Matrix2D<double>& trans_solve(const Matrix2D<double>& A, Matrix2D<double>& B) {
    long nrhs = static_cast<long>(B.size_q1);
    long n = static_cast<long>(A.size_q1);   
    long m = static_cast<long>(A.size_q2);   
    double *a {A[0]};
    double *b {B[0]};
    long lda {n};
    long ldb = static_cast<long>(B.size_q2);
    std::unique_ptr<long[]> ipiv {new long[m]};
    long info {0};
    long lwork {n};

/*     int dgetrf_(long *m, long *n, double *a, long * lda, 
                long *ipiv, long *info);

    int dgetri_(long *n, double *a, long *lda, long *ipiv, 
                double *work, long *lwork, long *info); */
    dgesv_(&n, &nrhs, a, &lda, ipiv.get(), b, &ldb, &info);
    return B;
}

Matrix1D<double> toMatrix1D(const Matrix2D<double>& A) {
    Matrix1D<double> result(A.size_q1 * A.size_q2);
    for (size_t i = 0; i < A.size_q1; ++i) {
        for (size_t j = 0; j < A.size_q2; ++j) {
            result[i * A.size_q2 + j] = A[i][j];
        }
    } 
    return result;
}

Matrix2D<double> toMatrix2D(
    const Matrix1D<double>& A, size_t size_q1, size_t size_q2
) {
    if (A.size_q1 != size_q1 * size_q2) {
        std::cout << "Error! In " << __FILE__ << ", line " << __LINE__ << ": ";
        std::cout << "Sizes mismatch:";
        std::cout << "A.size_q1 = " << A.size_q1 << ", ";
        std::cout << "size_q1 = " << size_q1 << ", "; 
        std::cout << "size_q2 = " << size_q2 << std::endl; 
        exit(EXIT_FAILURE);
    }

    Matrix2D<double> result(size_q1, size_q2);
    for (auto i = 0; i < size_q1; ++i) {
        for (auto j = 0; j < size_q2; ++j) {
            result[i][j] = A[i * size_q2 + j];
        }
    }

    return result;
}