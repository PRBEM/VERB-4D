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
    long m {A.size_q2};
    long n {A.size_q1};
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
    long m {B.size_q2};
    long n {A.size_q1};
    long k {B.size_q1};
    double alpha {1.0};
    double *a {A[0]};
    long lda {A.size_q2};
    double *b {B[0]};
    long ldb {B.size_q2};
    double beta {0.};
    double *c {answer[0]};
    long ldc {m};
    
    dgemm_(&transb, &transa, &m, &n, &k, &alpha, b, &ldb, a, &lda, &beta, 
            c, &ldc);

    return answer;
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
    long m {B.size_q1};
    long n {A.size_q1};
    long k {B.size_q2};
    double alpha {1.0};
    double *a {A[0]};
    long lda {k};
    double *b {B[0]};
    long ldb {B.size_q2};
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

    for (auto i = 0; i < result.size_q1; ++i) {
        for (auto j = 0; j < result.size_q2; ++j) {
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

    for (auto i = 0; i < result.size_q1; ++i) {
        for (auto j = 0; j < result.size_q2; ++j) {
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


    long m {A.size_q1};
    long n {A.size_q1};   
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

Matrix1D<double> toMatrix1D(const Matrix2D<double>& A) {
    Matrix1D<double> result(A.size_q1 * A.size_q2);
    for (auto i = 0; i < A.size_q1; ++i) {
        for (auto j = 0; j < A.size_q2; ++j) {
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