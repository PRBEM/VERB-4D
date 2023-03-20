#include "../src/DataAssimilationHelper.h"
#include "../src/DataAssimilation.h"
#include "../src/Parameters.h"
#include "../src/MatrixOperations.h"
#include "TestHelper.h"
#include "test_runner.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace data_assimilation::internal;
using namespace data_assimilation;

void TestDataAssimilationAuxiliaryFunctions() {
    {
        ASSERT(str2bool("True"));
        ASSERT(str2bool("TRUE"));
        ASSERT(str2bool("true"));
        ASSERT(!str2bool("false"));
        ASSERT(!str2bool("False"));
        ASSERT(!str2bool("FALSE"));
    } 

    {
        string filename {"__ParametersTmp__.ini"};
        string parameterName {"some_parameter"};
        string parameterName2 {"some_parameter2"};
        ofstream file(filename);
        if (file.is_open()) {
            file.clear();
            file << parameterName << " = true" << endl;
            file << parameterName2 << " = false" << endl;
            file.close();
            ::Parameters parameters {filename};
            remove(filename.c_str());
            
            bool result;
            bool result2;
            getParameterBool(parameters, parameterName, result, true);
            getParameterBool(parameters, parameterName2, result2, true);
            
            ASSERT(result);
            ASSERT(!result2);
        } else {
            cerr << "File " << filename << " was not open " << endl;
        }

    }
}

void TestPMFDataCat() {
    std::vector<ProcessedMatFileData> data(2);
    data[0].MLT = initialize({1., 2., 3.});
    data[0].R = initialize({10., 20., 30.});
    data[0].V = initialize({
        {
            {100., 200., 300.},
            {500., 600., 700.},
            {510., 610., 710.}
        },
        {
            {800., 900., 110.},
            {210., 310., 410.},
            {220., 320., 420.}
        }
    });
    data[0].K = initialize({
        {1000., 2000., 3000.},
        {4000., 5000., 6000.},
        {4100., 5100., 6100.}
    });
    data[0].PSD = initialize({
        {
            {10000., 20000., 30000.},
            {40000., 50000., 60000.}
        },
        {
            {70000., 80000., 90000.},
            {11000., 21000., 31000.}
        }
    });
    data[1].MLT = initialize({4., 5.});
    data[1].R = initialize({40., 50.});
    data[1].V = initialize({
        {
            {190., 290., 390.},
            {590., 690., 790.},
            {890., 990., 880.}
        },
        {
            {890., 990., 190.},
            {290., 390., 490.},
            {590., 690., 790.}
        }
    });
    data[1].K = initialize({
        {1900., 2900., 3900.},
        {4900., 5900., 6900.}
    });
    data[1].PSD = initialize({
        {
            {19000., 29000., 39000.},
            {40000., 59000., 69000.}
        }
    });

    ProcessedMatFileData expected;
    expected.MLT = initialize({1., 2., 3., 4., 5.});
    expected.R = initialize({10., 20., 30., 40., 50.});
    expected.V = initialize({
        {
            {100., 200., 300.},
            {500., 600., 700.},
            {510., 610., 710.}
        },
        {
            {800., 900., 110.},
            {210., 310., 410.},
            {220., 320., 420.}
        },
        {
            {190., 290., 390.},
            {590., 690., 790.},
            {890., 990., 880.}
        },
        {
            {890., 990., 190.},
            {290., 390., 490.},
            {590., 690., 790.}
        }
    });
    expected.K = initialize({
        {1000., 2000., 3000.},
        {4000., 5000., 6000.},
        {4100., 5100., 6100.},
        {1900., 2900., 3900.},
        {4900., 5900., 6900.}
    });
    expected.PSD = initialize({
        {
            {10000., 20000., 30000.},
            {40000., 50000., 60000.}
        },
        {
            {70000., 80000., 90000.},
            {11000., 21000., 31000.}
        },
        {
            {19000., 29000., 39000.},
            {40000., 59000., 69000.}
        }
    });
    auto result = cat(data);
    TestEqualTol(result.MLT, expected.MLT);
    TestEqualTol(result.R, expected.R);
    TestEqualTol(result.V, expected.V);
    TestEqualTol(result.K, expected.K);
    TestEqualTol(result.PSD, expected.PSD);
}

void TestInterpolation() {
   // to do: test four corner interpolation
   ASSERT_EQUAL(0,0);     
}

void TestKalmanFilter() {
    Matrix1D<double> forecast = initialize(
        {1., 2., 3.}
    );
    Matrix2D<double> Pa = initialize(
        {
            {0.8 , 0.2, 0.01},
            {0.1 , 0.9, 0.25},
            {0.02, 0.5, 0.83}
        }
    );
    Matrix2D<double> M = initialize(
        {
            {4., 0., 0.},
            {1., 4., 1.},
            {0., 1., 4.}
        }
    );
    Matrix2D<double> Q = diag(forecast * 0.5);
    Matrix1D<double> obs = initialize(
        {1.15, 1.9}
    );
    Matrix2D<double> H = initialize(
        {
            {1., 0., 0.},
            {0., 1., 0.}
        }
    );
    Matrix2D<double> R = diag(obs * 0.5);
    runKalmanFilter(forecast, Pa, M, Q, obs, H, R);
    Matrix1D<double> forecast_expected = initialize(
        {
            1.141739752121648,
            1.907275963366603,
            2.892779016721860
        }
    );
    Matrix2D<double> Pa_expected = initialize( 
        {
            {0.548465083488248 ,  0.012711909576093,  -0.111720190195332},
            {0.009632627132196 ,  0.904750510910109,   0.539295640565006},
            {-0.122416645000452,  0.715218960054244,  10.106010915435636}
        }
    );
    PRINT(Pa);
    PRINT(Pa_expected);
    TestEqualTol(forecast, forecast_expected, 1e-14);
    TestEqualTol(Pa, Pa_expected, 1e-14);
}

void TestModelMatrixConvectionWithoutStabilityCheck() {
    double dP = M_PI * 2. / 3.;
    double dR = 0.5;
    double dt = 4.;
    Matrix2D<double> VP = initialize(
        {
            {1., 2., 3.},
            {4., 5., 6.},
            {-7., -8., -9.},
            {1., 2., 3.}
        }
    );
    Matrix2D<double> VR = initialize(
        {
            {10., 20., 30.},
            {-40., -50., -60.},
            {70., 80., 90.},
            {10., 20., 30.}
        }
    );
    Matrix2D<double> Loss = initialize(
        {
            {-1., -2., -3.},
            {-5., -6., -7.},
            {-8., -9., -4.},
            {-2., -1.,  0.}
        }
    ); 
    
    auto CR = VR * dt / dR;
    auto CP = VP * dt / dP;

    auto sign = [](double val) {
        return val == 0. 
            ? 0. 
            : (val > 0.) ? 1. : -1.;
    };

    Matrix2D<double> expected (CP.size_q1 * CP.size_q2, 
        CP.size_q1 * CP.size_q2);
    expected = 0.;
    expected[0][0] = 1.;
    expected[1][7] = 0.5 * std::abs(CP[0][1]) * (1. + sign(CP[0][1])); 
    expected[1][0] = 0.5 * std::abs(CR[0][1]) * (1. + sign(CR[0][1])); 
    expected[1][1] = 1. - std::abs(CR[0][1]) - std::abs(CP[0][1])
                     + Loss[0][1] * dt;
    expected[1][2] = 0.5 * std::abs(CR[0][1]) * (1. - sign(CR[0][1])); 
    expected[1][5] = 0.5 * std::abs(CP[0][1]) * (1. - sign(CP[0][1])); 
    expected[2][2] = 1.;
    expected[3][3] = 1.;
    expected[4][1] = 0.5 * std::abs(CP[1][1]) * (1. + sign(CP[1][1]));
    expected[4][3] = 0.5 * std::abs(CR[1][1]) * (1. + sign(CR[1][1]));
    expected[4][4] = 1. - std::abs(CR[1][1]) - std::abs(CP[1][1]) 
                     + Loss[1][1] * dt;
    expected[4][5] = 0.5 * std::abs(CR[1][1]) * (1. - sign(CR[1][1]));
    expected[4][7] = 0.5 * std::abs(CP[1][1]) * (1. - sign(CP[1][1]));
    expected[5][5] = 1.;
    expected[6][6] = 1.;
    expected[7][4] = 0.5 * std::abs(CP[2][1]) * (1. + sign(CP[2][1]));
    expected[7][6] = 0.5 * std::abs(CR[2][1]) * (1. + sign(CR[2][1]));
    expected[7][7] = 1. - std::abs(CR[2][1]) - std::abs(CP[2][1])
                     + Loss[2][1] * dt;
    expected[7][8] = 0.5 * std::abs(CR[2][1]) * (1. - sign(CR[2][1]));
    expected[7][10] = 0.5 * std::abs(CP[2][1]) * (1. - sign(CP[2][1]));
    expected[8][8] = 1.;
    expected[9][9] = 1.;
    expected[10][7] = 0.5 * std::abs(CP[3][1]) * (1. + sign(CP[3][1]));
    expected[10][9] = 0.5 * std::abs(CR[3][1]) * (1. + sign(CR[3][1]));
    expected[10][10] = 1. - std::abs(CR[3][1]) - std::abs(CP[3][1]) 
                       + Loss[3][1] * dt;
    expected[10][11] = 0.5 * std::abs(CR[3][1]) * (1. - sign(CR[3][1]));
    expected[10][1] = 0.5 * std::abs(CP[3][1]) * (1. - sign(CP[3][1]));
    expected[11][11] = 1.; 

    auto result = getModelMatrixConvection2DNoStabilityCheck(
        VP, VR, Loss, dt, dP, dR);

    TestEqualTol(expected, result);
}

void TestModelMatrixConvectionWithStabilityCheck() {
    double dP = M_PI * 2. / 3.;
    double dR = 0.5;
    double dt = 4.;
    Matrix2D<double> VP = initialize(
        {
            {1., 2., 3.},
            {4., 5., 6.},
            {-7., -8., -9.},
            {1., 2., 3.}
        }
    );
    Matrix2D<double> VR = initialize(
        {
            {10., 20., 30.},
            {-40., -50., -60.},
            {70., 80., 90.},
            {10., 20., 30.}
        }
    );
    Matrix2D<double> Loss = initialize(
        {
            {-1., -2., -3.},
            {-5., -6., -7.},
            {-8., -9., -4.},
            {-2., -1.,  0.}
        }
    ); 

    double maximumCourantNumber {0.5};
    auto splitStep = splitTimeStepCourantCondition(maximumCourantNumber, 
           dt, VP, VR, dP, dR); 

    const int& nt = splitStep.first;
    const double& dt_new = splitStep.second;

    ASSERT_EQUAL(nt, 720*2 + 1);
    ASSERT(std::abs(dt_new - dt / 1441.) < 1e-15);

    auto CR = VR * dt_new / dR;
    auto CP = VP * dt_new / dP;

    auto sign = [](double val) {
        return val == 0. 
            ? 0. 
            : (val > 0.) ? 1. : -1.;
    };

    auto expectedOneStep = getModelMatrixConvection2DNoStabilityCheck(
            VP, VR, Loss, dt_new, dP, dR);

    Matrix2D<double> expected(expectedOneStep);
    for (auto i = 0; i < nt; ++i) {
        expected = expected * expectedOneStep;
    }

    auto result = getModelMatrixConvection2D(
        VP, VR, Loss, dt, dP, dR);

    TestEqualTol(expected, result, 1e-12);
}

void TestConditionalMean() {
    Matrix1D<double> a = initialize(
        {1., 2., 4., 5., 6., 7., 8.}
    );
    Matrix1D<bool> flags = initialize(
        {true, false, false, true, false, true, true}
    );

    auto result = conditionalMean(a, flags); 
    auto expected = 5.25;
    TestEqualTol(result, expected);

    flags = false;
    result = conditionalMean(a, flags);
    // expected = FILLVAL;
    // ASSERT_EQUAL(result, expected);
    ASSERT(isnan(result));
}

void TestBin() {
    Matrix2D<double> P = initialize(
        {
            {10., 10., 10., 10., 10.},
            {20., 20., 20., 20., 20.},
            {30., 30., 30., 30., 30.},
            {40., 40., 40., 40., 40.}
        }
    ); 
    Matrix2D<double> R = initialize(
        {
            {1., 2., 3., 4., 5.},
            {1., 2., 3., 4., 5.},
            {1., 2., 3., 4., 5.},
            {1., 2., 3., 4., 5.}
        }
    );
   
    Observations data;

    data.R = initialize(
        {
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9,
         1.4, 1.3, 1.6, 2.1, 2.8, 3.2, 3.7, 4.1, 4.8, 4.9
        }
    );

    data.P = initialize(
        {
         11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 
         13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 
         17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 
         22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 
         29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 
         33., 33., 33., 33., 33., 33., 33., 33., 33., 33., 
         38., 38., 38., 38., 38., 38., 38., 38., 38., 38., 
         39., 39., 39., 39., 39., 39., 39., 39., 39., 39. 
        }
    );

    data.PSD = initialize(
        {
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.,
         10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
        }
    );

    auto result = bin(data, P, R);

    Matrix2D<double> expected = initialize(
        {
            {15., 35., 55., 75., 95.},
            {15., 35., 55., 75., 95.},
            {15., 35., 55., 75., 95.},
            {15., 35., 55., 75., 95.}
        }
    );

    TestEqualTol(result, expected, 1e-12);
    
    result = bin(data, P, R, "linear");
    TestEqualTol(result, expected, 1e-12);

    Matrix1D<double> logs = initialize(
        {
            std::pow(10., log10(10. * 20.)), 
            std::pow(10., log10(30. * 40.)), 
            std::pow(10., log10(50. * 60.)), 
            std::pow(10., log10(70. * 80.)), 
            std::pow(10., log10(90. * 100.)) 
        }
    );

    result = bin(data, P, R, "log10");
    expected = initialize(
        {
            {logs[0], logs[1], logs[2], logs[3], logs[4]},
            {logs[0], logs[1], logs[2], logs[3], logs[4]},
            {logs[0], logs[1], logs[2], logs[3], logs[4]},
            {logs[0], logs[1], logs[2], logs[3], logs[4]}
        }
    );
}

void TestConvertToObservationSpace() {
    Matrix2D<double> data = initialize(
        {
            {FILLVAL, 1., FILLVAL},
            {2., FILLVAL, 3.}
        }
    );

    ObservationSpace expected;
    expected.data = initialize({1., 2., 3.});
    expected.H = initialize(
        {
            {0., 1., 0., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0.},
            {0., 0., 0., 0., 0., 1.},
        }
    );

    auto result = convertToObservationSpace(data);

    TestEqualTol(result.data, expected.data);
    TestEqualTol(result.H, expected.H);
}

void TestDataAssimilationInternal() {
    TestRunner tr;
    RUN_TEST(tr, TestDataAssimilationAuxiliaryFunctions);
    RUN_TEST(tr, TestKalmanFilter);
    RUN_TEST(tr, TestModelMatrixConvectionWithoutStabilityCheck);
    RUN_TEST(tr, TestModelMatrixConvectionWithStabilityCheck);
    RUN_TEST(tr, TestConditionalMean);
    RUN_TEST(tr, TestBin);
    RUN_TEST(tr, TestConvertToObservationSpace);
    RUN_TEST(tr, TestInterpolation);
    RUN_TEST(tr, TestPMFDataCat);
}
