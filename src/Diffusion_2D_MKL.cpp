#include "Matrix.h"
#include "BoundaryConditionType.hpp"
#include "Logger.h"
#include <vector>
typedef Matrix2D<double> mat2d;
void initialize_sparse_values(
        const mat2d& v_grid, const mat2d& k_grid,
        BoundaryConditionType v_lower, BoundaryConditionType v_upper,
        BoundaryConditionType k_lower, BoundaryConditionType k_upper,
        const mat2d& Dvv, const mat2d& Dvk, const mat2d& Dkv, mat2d& Dkk, 
        const mat2d& jacobian, const mat2d& loss, double dt,
        std::vector<double>& values)
{
    int v_size = v_grid.size_q1;
    int k_size = k_grid.size_q2;
    values.clear();
    values.reserve(
        v_size * (k_size - 2) * 9
    );

    double v_offdiag_lower, v_offdiag_upper;
    double k_offdiag_lower, k_offdiag_upper;
    switch(v_lower)
    {
        case BoundaryConditionType::ConstantDerivative:
            v_offdiag_lower = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            v_offdiag_lower = 0;
            break;
        default:
            Logger::error << "Lower diffusion boundary condition in V must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }
    switch(v_upper)
    {
        case BoundaryConditionType::ConstantDerivative:
            v_offdiag_upper = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            v_offdiag_upper = 0;
            break;
        default:
            Logger::error << "Upper diffusion boundary condition in V must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }
    switch(k_lower)
    {
        case BoundaryConditionType::ConstantDerivative:
            k_offdiag_lower = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            k_offdiag_lower = 0;
            break;
        default:
            Logger::error << "Lower diffusion boundary condition in K must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }
    switch(k_upper)
    {
        case BoundaryConditionType::ConstantDerivative:
            k_offdiag_upper = 1;
            break;
        case BoundaryConditionType::ConstantValue:
            k_offdiag_upper = 0;
            break;
        default:
            Logger::error << "Upper diffusion boundary condition in K must be either constant value or constant derivative\n";
            exit(EXIT_FAILURE);
    }

    // lower K boundary, lower V boundary
    values.push_back(1);
    values.push_back(-v_offdiag_lower);
    values.push_back(-k_offdiag_lower);
    values.push_back(k_offdiag_lower * v_offdiag_lower);

    // lower K boundary, inner V
    for(int i = 1; i < v_size - 1; i++)
    {
        values.push_back(1);
        values.push_back(-k_offdiag_lower);
    }

    // lower K boundary, upper V boundary
    values.push_back(-v_offdiag_upper);
    values.push_back(1);
    values.push_back(k_offdiag_lower * v_offdiag_upper);
    values.push_back(-k_offdiag_lower);

    // inner K
    for(int j = 1; j < k_size - 1; j++)
    {
        // lower V boundary
        values.push_back(1);
        values.push_back(-v_offdiag_lower);

        // inner V
        for(int i = 1; i < v_size - 1; i++)
        {
            double dv_left  = v_grid[i][j] - v_grid[i-1][j];
            double dv_right = v_grid[i+1][j] - v_grid[i][j];
            double dk_left  = k_grid[i][j] - v_grid[i][j-1];
            double dk_right = k_grid[i][j+1] - v_grid[i][j];

            // f(i-1,j-1) coeff
            double coeff = 0.25 * (jacobian[i-1][j] * Dvk[i-1][j] + jacobian[i][j-1] * Dkv[i][j-1]) / (jacobian[i][j] * dv_left * dk_left);
            values.push_back(-coeff);

            // f(i-1,j) coeff
            double factor = 0.25 / (jacobian[i][j] * dv_left);
            coeff = factor * (
                2 * (jacobian[i][j]   * Dvv[i][j]   / dv_right + jacobian[i-1][j] * Dvv[i-1][j] / dv_left)
                  + (jacobian[i-1][j] * Dvk[i-1][j] / dk_right - jacobian[i-1][j] * Dvk[i-1][j] / dk_left)
                  + (jacobian[i][j]   * Dkv[i][j]   / dk_right - jacobian[i][j]   * Dkv[i][j]   / dk_left)
            );
            values.push_back(-coeff);
                
            // f(i-1,j+1) coeff
            coeff = -0.25 * (jacobian[i-1][j] * Dvk[i-1][j] + jacobian[i][j+1] * Dkv[i][j+1]) / (jacobian[i][j] * dv_left * dk_right);
            values.push_back(-coeff);

            // f(i,j-1) coeff
            factor = 0.25 / (jacobian[i][j] * dk_left);
            coeff = factor * (
                2 * (jacobian[i][j]   * Dkk[i][j]   / dk_right + jacobian[i][j-1] * Dkk[i][j-1] / dk_left)
                  + (jacobian[i][j]   * Dvk[i][j]   / dv_right - jacobian[i][j]   * Dvk[i][j]   / dv_left)
                  + (jacobian[i][j-1] * Dkv[i][j-1] / dv_right - jacobian[i][j-1] * Dkv[i][j-1] / dv_left)
            );
            values.push_back(-coeff);        
            
            // f(i,j) coeff, 1/dt and loss
            factor = -0.5 / jacobian[i][j];
            coeff = factor * (
                jacobian[i+1][j] * Dvv[i+1][j] / (dv_right * dv_right) + 2 * jacobian[i][j] * Dvv[i][j] / (dv_left * dv_right) + jacobian[i-1][j] * Dvv[i-1][j] / (dv_left * dv_left)
              + jacobian[i][j+1] * Dkk[i][j+1] / (dk_right * dk_right) + 2 * jacobian[i][j] * Dkk[i][j] / (dk_left * dk_right) + jacobian[i][j-1] * Dkk[i][j-1] / (dk_left * dk_left)
              - (1 / (dv_right * dk_right) - 1 / (dv_right * dk_left) - 1 / (dv_left * dk_right) + 1 / (dv_left * dk_left)) * (Dvk[i][j] + Dkv[i][j]) * jacobian[i][j] / 2
            ); 
            values.push_back(1/dt - loss[i][j] - coeff);

            // f(i,j+1) coeff
            factor = 0.25 / (jacobian[i][j] * dk_right);
            coeff = factor * (
                2 * (jacobian[i][j+1] * Dkk[i][j+1] / dk_right + jacobian[i][j]   * Dkk[i][j]   / dk_left)
                  +  jacobian[i][j]   * Dvk[i][j]   / dv_left  - jacobian[i][j]   * Dvk[i][j]   / dv_right
                  +  jacobian[i][j]   * Dkv[i][j]   / dv_left  - jacobian[i][j+1] * Dkv[i][j+1] / dv_right
            );
            values.push_back(-coeff);

            // f(i+1,j-1) coeff
            coeff = -0.25 * (jacobian[i+1][j] * Dvk[i+1][j] + jacobian[i][j-1] * Dkv[i][j-1]) / (jacobian[i][j] * dv_right * dk_left);
            values.push_back(-coeff);
            
            // f(i+1,j) coeff
            factor = 0.25 / (jacobian[i][j] * dv_right);
            coeff = factor * (
                2 * (jacobian[i+1][j] * Dvv[i+1][j] / dv_right + jacobian[i][j]   * Dvv[i][j]   / dv_left)
                  + (jacobian[i+1][j] * Dvk[i+1][j] / dk_left  - jacobian[i+1][j] * Dvk[i+1][j] / dk_right)
                  + (jacobian[i][j]   * Dkv[i][j]   / dk_left  - jacobian[i][j]   * Dkv[i][j]   / dk_right)
            );
            values.push_back(-coeff);
            
            // f(i+1,j+1) coeff
            coeff = 0.25 * (jacobian[i+1][j] * Dvk[i+1][j] + jacobian[i][j+1] * Dkv[i][j+1]) / (jacobian[i][j] * dv_right * dk_right);
            values.push_back(-coeff);
        }

        // upper V boundary
        values.push_back(-v_offdiag_lower);
        values.push_back(1);
    }
    // upper K boundary, lower V boundary
    values.push_back(-k_offdiag_upper);
    values.push_back(k_offdiag_upper * v_offdiag_lower);
    values.push_back(1);
    values.push_back(-v_offdiag_lower);

    // upper K boundary, inner V
    for(int i = 1; i < v_size - 1; i++)
    {
        values.push_back(-k_offdiag_upper);
        values.push_back(1);
    }

    // upper K boundary, upper V boundary
    values.push_back(k_offdiag_upper * v_offdiag_upper);
    values.push_back(-k_offdiag_upper);
    values.push_back(-v_offdiag_upper);
    values.push_back(1);
}