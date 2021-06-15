#include "test_runner.h"
#include "../src/Matrix.h"
#include "../src/MatrixOperations.h"
#include "MatrixOperationsTests.h"
#include "TestHelper.h"

#include <cmath>
#include <initializer_list>

void TestCat1D() {
    Matrix1D<double> a = initialize({1., 2., 3., 4., 5.});
    Matrix1D<double> b = initialize({6., 7., 8., 9., 10., 11., 12.});
    
    {
        Matrix1D<double> expected = initialize(
            {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.});
        auto result = cat(a, b);
        TestEqual(result, expected);
    }

    {
        Matrix1D<double> expected = initialize(
            {6., 7., 8., 9., 10., 11., 12., 1., 2., 3., 4., 5.});
        auto result = cat(b, a);
        TestEqual(result, expected);
    }
}

void TestCat2D() {
    Matrix2D<double> a = initialize(
        {{1., 2., 3., 4., 5.},
         {5., 6., 7., 8., 9.}}
    );
    Matrix2D<double> b = initialize(
        {{11., 12., 13., 14., 15.},
         {15., 16., 17., 18., 19.}, 
         {20., 21., 22., 23., 24.}}
    );

    {
        Matrix2D<double> expected = initialize(
            {{1., 2., 3., 4., 5.},
             {5., 6., 7., 8., 9.},
             {11., 12., 13., 14., 15.},
             {15., 16., 17., 18., 19.}, 
             {20., 21., 22., 23., 24.}}
        );
        auto result = cat(a, b);
        TestEqual(result, expected);
    }

    {
        Matrix2D<double> expected = initialize(
            {{11., 12., 13., 14., 15.},
             {15., 16., 17., 18., 19.}, 
             {20., 21., 22., 23., 24.},
             {1., 2., 3., 4., 5.},
             {5., 6., 7., 8., 9.}}
        );
        auto result = cat(b, a);
        TestEqual(result, expected);
    }
}

void TestCat3D() {
    Matrix3D<double> a = initialize(
        {
         {{11., 12., 13., 14., 15.},
          {15., 16., 17., 18., 19.}},
         
         {{21., 22., 23., 24., 25.},
          {25., 26., 27., 28., 29.}},
         
         {{31., 32., 33., 34., 35.},
          {35., 36., 37., 38., 39.}}}
    );
    Matrix3D<double> b = initialize(
        {
         {{41., 42., 43., 44., 45.},
          {45., 46., 47., 48., 49.}},
                                 
         {{51., 52., 53., 54., 55.},
          {55., 56., 57., 58., 59.}}}
    );

    {
        Matrix3D<double> expected = initialize(
            {
            {{11., 12., 13., 14., 15.},
             {15., 16., 17., 18., 19.}},
                                    
            {{21., 22., 23., 24., 25.},
             {25., 26., 27., 28., 29.}},
                                    
            {{31., 32., 33., 34., 35.},
             {35., 36., 37., 38., 39.}},
                                    
            {{41., 42., 43., 44., 45.},
             {45., 46., 47., 48., 49.}},
                                    
            {{51., 52., 53., 54., 55.},
             {55., 56., 57., 58., 59.}}}
        );

        auto result = cat(a, b);
        TestEqual(result, expected);
    }

    {
        Matrix3D<double> expected = initialize(
            {
            {{41., 42., 43., 44., 45.},
             {45., 46., 47., 48., 49.}},
            
            {{51., 52., 53., 54., 55.},
             {55., 56., 57., 58., 59.}},

            {{11., 12., 13., 14., 15.},
             {15., 16., 17., 18., 19.}},
            
            {{21., 22., 23., 24., 25.},
             {25., 26., 27., 28., 29.}},
            
            {{31., 32., 33., 34., 35.},
             {35., 36., 37., 38., 39.}}}
            
        );

        auto result = cat(b, a);
        TestEqual(result, expected);
    }
}

void TestCat() {
    TestCat1D();
    TestCat2D();
    TestCat3D();
}


void TestTranspose() {

    Matrix2D<double> a = initialize(
        {{1., 2., 3., 4., 5.},
         {5., 6., 7., 8., 9.}}
    );

    Matrix2D<double> expected = initialize(
        {{1., 5.},
         {2., 6.}, 
         {3., 7.}, 
         {4., 8.},
         {5., 9.}}
    );

    auto result = transpose(a);

    TestEqual(result, expected);
}

void TestLog10() {    
    {
        Matrix1D<double> a = initialize(
            {1., 2., 3., 4., 5.}
        );

        Matrix1D<double> expected(a.size_q1);
        for (auto i = 0; i < expected.size_q1; ++i) {
            expected[i] = std::log10(a[i]);
        }

        TestEqual(log10(a), expected);
    }

    {
        Matrix2D<double> a = initialize(
            {{1., 2., 3., 4., 5.},
             {5., 6., 7., 8., 9.}}
        );

        Matrix2D<double> expected(a.size_q1, a.size_q2);
        for (auto i = 0; i < expected.size_q1; ++i) {
            for (auto j = 0; j < expected.size_q2; ++j) {
                expected[i][j] = std::log10(a[i][j]);
            }
        }

        TestEqual(log10(a), expected);
    }
}

void TestPow() {    
    {
        Matrix1D<double> a = initialize(
            {1., 2., 3., 4., 5.}
        );

        Matrix1D<double> expected(a.size_q1);
        for (auto i = 0; i < expected.size_q1; ++i) {
            expected[i] = std::pow(10., a[i]);
        }

        TestEqual(pow(10., a), expected);
    }
}

void TestDiag() {
    {
        Matrix1D<double> a = initialize(
            {1., 2., 3., 4., 5.}
        );

        Matrix2D<double> expected = initialize(
            {
                {1., 0., 0., 0., 0.},   
                {0., 2., 0., 0., 0.},   
                {0., 0., 3., 0., 0.},   
                {0., 0., 0., 4., 0.},   
                {0., 0., 0., 0., 5.}
            }   
        );

        TestEqual(diag(a), expected);
    }
}

void TestOperators() {
    {
        Matrix1D<double> a = initialize(
            {1., 2., 1., 1., 5.}
        );

        Matrix1D<bool> expected = initialize(                
            {true, false, true, true, false}
        );
        TestEqual(a == 1., expected);
        TestEqual(1. == a, expected);

        expected = initialize(                
            {false, true, false, false, true}
        );
        TestEqual(a > 1., expected);
        TestEqual(a >= 2., expected);
        TestEqual(1. < a, expected);
        TestEqual(2. <= a, expected);

        expected = initialize(                
            {true, false, true, true, false}
        );
        TestEqual(a < 2., expected);
        TestEqual(2. > a, expected);

        expected = initialize(                
            {true, true, true, true, false}
        );
        TestEqual(a <= 2., expected);
        TestEqual(2. >= a, expected);
    }
    {
        Matrix1D<bool> a = initialize(                
            {true, true, false, false}
        );        
        Matrix1D<bool> b = initialize(                
            {true, false, true, false}
        );       

        Matrix1D<bool> expected = initialize(
            {true, false, false, false}
        ); 
        TestEqual(a && b, expected);

        expected = initialize(
            {true, true, true, false}
        ); 
        TestEqual(a || b, expected);
    }
    {
        Matrix2D<double> a = initialize(
            {
                {1., 2., 1., 1., 5.},
                {7., 8., 0., 1., -5.},
            }
        );

        Matrix2D<bool> expected = initialize(                
            {
                {true, false, true, true, false},
                {false, false, false, true, false}
            }
        );
        TestEqual(a == 1., expected);
        TestEqual(1. == a, expected);

        expected = initialize(                
            {
                {false, true, false, false, true},
                {true, true, false, false, false}
            }
        );
        TestEqual(a > 1., expected);
        TestEqual(a >= 2., expected);
        TestEqual(1. < a, expected);
        TestEqual(2. <= a, expected);

        expected = initialize(                
            {
                {true, false, true, true, false},
                {false, false, true, true, true}
            }
        );
        TestEqual(a < 2., expected);
        TestEqual(2. > a, expected);

        expected = initialize(                
            {
                {true, true, true, true, false},
                {false, false, true, true, true}
            }
        );
        TestEqual(a <= 2., expected);
        TestEqual(2. >= a, expected);
    }
    {
        Matrix2D<bool> a = initialize(
            {
                {true, true, false, false},
                {false, false, true, true}
            }
        );
        Matrix2D<bool> b = initialize(
            {
                {true, false, true, false},
                {false, true, false, true}
            }
        );

        Matrix2D<bool> expected = initialize(                
            {
                {true, false, false, false},
                {false, false, false, true}
            }
        );
        TestEqual(a && b, expected);
    }
    {
        Matrix1D<double> a = initialize(                
            {1., 2., 3., 4., 5.}
        );
        Matrix1D<double> b = initialize(                
            {-1., 2., -5., 8., 5.}
        );

        Matrix1D<double> expected = initialize(
            {2., 0., 8., -4., 0.}
        );
        TestEqual(a - b, expected);
    }
    {
        Matrix2D<double> A = initialize(                
            {
                {1., 2., 3., 4.},
                {5., 6., 7., 8.},
                {9., 1., 2., 3.}
            }
        );
        Matrix1D<double> x = initialize (
            {1., 2., 3., 4.}
        );
        Matrix1D<double> expected = initialize (                
            {
                1. + 2. * 2. + 3. * 3. + 4. * 4.,
                5. + 6. * 2. + 7. * 3. + 8. * 4.,
                9. + 1. * 2. + 2. * 3. + 3. * 4.,
            }
        );
        auto result = A * x;
        TestEqual(result, expected);
    }
    {
        Matrix2D<double> A = initialize(                
            {
                {1., 2., 3., 4.},
                {5., 6., 7., 8.},
                {9., 1., 2., 3.}
            }
        );
        Matrix2D<double> B = initialize(                
            {
                {10., 20., 30., 40.},
                {50., 60., 70., 80.},
                {90., 10., 20., 30.},
                {40., 50., 60., 70.}
            }
        );
        Matrix2D<double> expected = initialize (                
            {
                {540., 370., 470., 570.},
                {1300., 930., 1190., 1450.},
                {440., 410., 560., 710.}
            }
        );
        auto result = A * B;
        TestEqual(result, expected);
        
        Matrix2D<double> C = initialize(                
            {
                {10., 20., 30., 40.},
                {50., 60., 70., 80.},
                {40., 10., 20., 30.}
            });
        Matrix2D<double> expected_ACT = initialize (                
            {
                {300., 700., 240.},
                {700., 1740., 640.},
                {290., 890., 500.}
            }
        );
        result = abtrans(A, C);
        TestEqual(result, expected_ACT);
    }
        

    {
        Matrix2D<double> A = initialize(                
            {
                {1., 2., 3., 4.},
                {5., 6., 7., 8.},
                {9., 1., 2., 3.}
            }
        );
        Matrix2D<double> B = initialize(                
            {
                {10., 20., 30., 40.},
                {50., 60., 70., 80.},
                {90., 10., 20., 30.}
            }
        );
        Matrix2D<double> expected = initialize (                
            {
                {11., 22., 33., 44.},
                {55., 66., 77., 88.},
                {99., 11., 22., 33.}
            }
        );
        auto result = A + B;
        TestEqual(result, expected);
    }
    {
        Matrix2D<double> A = initialize(                
            {
                {-1., -2., -3., -4.},
                {-5., -6., -7., -8.},
                {-9., -1., -2., -3.}
            }
        );
        Matrix2D<double> B = initialize(                
            {
                {10., 20., 30., 40.},
                {50., 60., 70., 80.},
                {90., 10., 20., 30.}
            }
        );
        Matrix2D<double> expected = initialize (                
            {
                {11., 22., 33., 44.},
                {55., 66., 77., 88.},
                {99., 11., 22., 33.}
            }
        );
        auto result = B - A;
        TestEqual(result, expected);
    }
}

void TestMatrixInversion() {
    {
        Matrix2D<double> A = diag(initialize({1., 1., 1., 1.}));
        TestEqual(inv(A), A);
        TestEqual(inv(A * 2.), A / 2.);
    }
    {
        Matrix2D<double> A = initialize(
            {
                {4., 7., 0., 0.},
                {2., 6., 0., 0.},
                {0., 0., 4., 7.},
                {0., 0., 2., 6.}
            }
        ); 
        Matrix2D<double> expected = initialize(
            {
                {0.6, -0.7, 0., 0.},
                {-0.2, 0.4, 0., 0.},
                {0., 0., 0.6, -0.7},
                {0., 0., -0.2, 0.4}
            }
        );
        TestEqualTol(inv(A), expected);
    }
    {
        double a = 1./std::sqrt(2.);
        Matrix2D<double> A = initialize(
            {
                {a, -a},
                {a, a}
            }
        );
        Matrix2D<double> expected = initialize(
            {
                {a, a},
                {-a, a}
            }
        );
        TestEqualTol(inv(A), expected);
    }
}

void TestConversions() {
    {
        Matrix2D<double> A = initialize(                
            {
                {-1., -2., -3., -4.},
                {-5., -6., -7., -8.},
                {-9., -1., -2., -3.}
            }
        );
        Matrix1D<double> expected = initialize(                
            {
                -1., -2., -3., -4.,
                -5., -6., -7., -8.,
                -9., -1., -2., -3.
            }
        );
        TestEqual(toMatrix1D(A), expected);
    }
    {
        Matrix1D<double> a = initialize(                
            {
                -1., -2., -3., -4.,
                -5., -6., -7., -8.,
                -9., -1., -2., -3.
            }
        );
        Matrix2D<double> expected = initialize(                
            {
                {-1., -2., -3., -4.},
                {-5., -6., -7., -8.},
                {-9., -1., -2., -3.}
            }
        );
        TestEqual(toMatrix2D(a, 3, 4), expected);
    }
}

void TestMatrixOperations() {
    TestRunner tr;
    RUN_TEST(tr, TestCat);
    RUN_TEST(tr, TestTranspose);
    RUN_TEST(tr, TestLog10);
    RUN_TEST(tr, TestPow);
    RUN_TEST(tr, TestDiag);
    RUN_TEST(tr, TestOperators);
    RUN_TEST(tr, TestMatrixInversion);
    RUN_TEST(tr, TestConversions);
}
