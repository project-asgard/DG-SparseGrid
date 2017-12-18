#include <iostream>
#include <vector>
#include "amgx_c.h"

int main(int argc, char** argv) {
    // Initialize AMGX
    AMGX_initialize();
    AMGX_initialize_plugins();    

    // Use default solver configuration
    AMGX_config_handle config;
    AMGX_config_create(&config, "");

    // Simple GPU resource handle
    AMGX_resources_handle resource;
    AMGX_resources_create_simple(&resource, config);

    // Create a solver based upon the config
    AMGX_solver_handle solver;
    AMGX_solver_create(&solver, resource, AMGX_mode_dDDI, config);

    // Create double precision device matrix/vectors for system: A*x = b
    
    // A, 4 x 4 identity matrix
    AMGX_matrix_handle A;
    AMGX_matrix_create(&A, resource, AMGX_mode_dDDI);
    std::vector<double> A_data = {1.0, 1.0, 1.0, 1.0};
    std::vector<int> A_row_ptr = {0,1,2,3,4};
    std::vector<int> A_column_indices = {0,1,2,3};
    AMGX_matrix_upload_all(A, 4, 4, 1, 1, A_row_ptr.data(), A_column_indices.data(), A_data.data(), 0);

    // b(right hand side)
    AMGX_vector_handle b;
    AMGX_vector_create(&b, resource, AMGX_mode_dDDI);
    std::vector<double> b_data = {1.0, 0.5, 0.25, 0.125};
    AMGX_vector_upload(b, 4, 1, b_data.data());

    // Initial guess for x is zero vector
    AMGX_vector_handle x;
    AMGX_vector_create(&x, resource, AMGX_mode_dDDI);
    AMGX_vector_set_zero(x, 4, 1);

    // setup and kickoff Solver
    AMGX_solver_setup(solver, A);
    AMGX_solver_solve(solver, b, x);

    // Download solution vector from solver
    std::vector<double> x_data(4);
    AMGX_vector_download(x, x_data.data());

    // Print results
    std::cout<<x_data[0]<<std::endl;
    std::cout<<x_data[1]<<std::endl;
    std::cout<<x_data[2]<<std::endl;
    std::cout<<x_data[3]<<std::endl;

    // Cleanup matrix/vectors
    AMGX_matrix_destroy(A);
    AMGX_vector_destroy(x);
    AMGX_vector_destroy(b);

    // Destroy solver components
    AMGX_solver_destroy(solver);
    AMGX_resources_destroy(resource);
    AMGX_config_destroy(config);

    // Finalize AMGX
    AMGX_finalize_plugins();
    AMGX_finalize();

    return 0;
}
