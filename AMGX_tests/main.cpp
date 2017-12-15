#include <iostream>
#include "amgx_c.h"

int main(int argc, char** argv) {
    // Initialize AMGX
    AMGX_initialize();
    AMGX_initialize_plugins();    

    // Use default solver configuration
    AMGX_config_handle config;
    AMGX_config_create(&config, "");

    // Simple GPU resource handle
    AMGX_resources_handle rsrc;
    AMGX_resources_create_simple(&rsrc, config);

    // Create a solver based upon the config
    AMGX_solver_handle solver;
    AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, config);

    //
    // Create double precision device matrix/vectors for system: matrix * soln = rhs
    // After creating on the host they will be sent to the solver
    // pinned memory could also be used
    //
    
    // matrix, 4x4 identity
    AMGX_matrix_handle matrix;
    AMGX_matrix_create(&matrix, rsrc, AMGX_mode_dDDI);
    double data[] = {1, 1, 1, 1};
    int col_ind[] = {0, 1, 2, 3};
    int row_ptr[] = {0, 1, 2, 3, 4};
    AMGX_matrix_upload_all(matrix, 4, 4, 1, 1, row_ptr, col_ind, data, 0);

    // rhs, set to all 1's
    AMGX_vector_handle rhs;
    AMGX_vector_create(&rhs, rsrc, AMGX_mode_dDDI);
    double rhs_data[] = {1,2,1,1};
    AMGX_vector_upload(rhs, 4, 1, rhs_data);

    // Initial guess is zero for simplicity
    AMGX_vector_handle soln;
    AMGX_vector_create(&soln, rsrc, AMGX_mode_dDDI);
    AMGX_vector_set_zero(soln, 4, 1);

    // setup and kickoff Solver
    AMGX_solver_setup(solver, matrix);
    AMGX_solver_solve(solver, rhs, soln);

    // Download solution vector from solver
    double soln_data[4];
    AMGX_vector_download(soln, soln_data);

    // Print results
    std::cout<<soln_data[0]<<std::endl;
    std::cout<<soln_data[1]<<std::endl;
    std::cout<<soln_data[2]<<std::endl;
    std::cout<<soln_data[3]<<std::endl;

    // Cleanup matrix/vectors
    AMGX_matrix_destroy(matrix);
    AMGX_vector_destroy(rhs);
    AMGX_vector_destroy(soln);

    // Destroy solver components
    AMGX_solver_destroy(solver);
    AMGX_resources_destroy(rsrc);
    AMGX_config_destroy(config);

    // Finalize AMGX
    AMGX_finalize_plugins();
    AMGX_finalize();

    return 0;
}
