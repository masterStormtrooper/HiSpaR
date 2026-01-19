// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <omp.h>

// Include HiSpa headers
#include "chromosome.h"

using namespace Rcpp;

//' Run HiSpa MCMC Analysis (Internal)
//' 
//' @description
//' Internal C++ function. Use hispa_analyze() from R instead.
//' All results are saved to output_dir by the Chromosome class.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.hispa_analyze_cpp)]]
bool hispa_analyze_cpp(
    std::string input_file,
    std::string output_dir,
    int mcmc_iterations = 6000,
    int num_clusters = 0,
    int mcmc_burn_in = 0,
    double mcmc_initial_sd = 0.1,
    double mcmc_sd_floor = 0.0001,
    double mcmc_sd_ceil = 0.3,
    bool use_cluster_init = false,
    int cluster_init_iterations = 1000,
    double cluster_initial_sd = 0.1,
    bool save_samples = false,
    int sample_interval = 50,
    bool verbose = true
) {
    try {
        // Make explicit copies of string parameters
        std::string input_file_copy(input_file);
        std::string output_dir_copy(output_dir);
        
        if (verbose) {
            Rcpp::Rcout << "HiSpa - Hi-C Spatial Analysis\n";
            Rcpp::Rcout << "==============================\n";
            Rcpp::Rcout << "Input file: " << input_file_copy << "\n";
            Rcpp::Rcout << "Output directory: " << output_dir_copy << "\n\n";
        }
        
        // Scope block to ensure Chromosome destructs before spdlog::shutdown()
        {
            // Initialize Chromosome object with output directory
            Chromosome my_chromosome(output_dir_copy);
            
            // Load contact matrix from file
            if (verbose) {
                Rcpp::Rcout << "Loading contact matrix from: " << input_file_copy << "\n";
            }
            bool load_result = my_chromosome.load_data_from_file(input_file_copy);
            if (!load_result) {
                Rcpp::stop("Failed to load contact matrix from file");
            }
            
            // Set MCMC options (using defaults from hispa_main.cpp)
            my_chromosome.set_skip_zero_contact_loci(true);
            my_chromosome.set_sample_from_prior(false);
            
            // --- 1. Pre-processing: Assign clusters ---
            if (verbose) {
                Rcpp::Rcout << "\n=== Pre-processing ===\n";
                Rcpp::Rcout << "Assigning clusters...\n";
            }
            
            if (num_clusters > 0) {
                my_chromosome.assign_clusters(num_clusters);
            } else {
                my_chromosome.assign_clusters();  // Auto-detect sqrt(n)
            }
            
            // Build cluster relationships (using default distance threshold = 2)
            if (verbose) {
                Rcpp::Rcout << "Building cluster relationships...\n";
            }
            my_chromosome.build_cluster_relationships_by_distance(2);
            
            // --- 2. Initialize structure ---
            if (verbose) {
                Rcpp::Rcout << "\n=== Structure Initialization ===\n";
            }
            
            if (use_cluster_init) {
                if (verbose) {
                    Rcpp::Rcout << "Initializing structure from individual cluster MCMC...\n";
                    Rcpp::Rcout << "Cluster MCMC iterations: " << cluster_init_iterations << "\n";
                    Rcpp::Rcout << "Cluster initial SD: " << cluster_initial_sd << "\n";
                }
                my_chromosome.initialize_structure_from_clusters(
                    cluster_init_iterations, 0, cluster_initial_sd);
                
                // --- 3. Assemble global structure ---
                if (verbose) {
                    Rcpp::Rcout << "Assembling global structure...\n";
                }
                my_chromosome.assemble_global_structure();
            } else {
                if (verbose) {
                    Rcpp::Rcout << "Using random position initialization...\n";
                }
                my_chromosome.initialize_positions();
            }
            
          
            // --- 4. Run main MCMC ---
            if (verbose) {
                Rcpp::Rcout << "\n=== Main MCMC Algorithm ===\n";
                Rcpp::Rcout << "Running MCMC with " << mcmc_iterations << " iterations...\n";
                Rcpp::Rcout << "Burn-in: " << mcmc_burn_in << "\n";
                Rcpp::Rcout << "Initial SD: " << mcmc_initial_sd << "\n";
                Rcpp::Rcout << "SD floor: " << mcmc_sd_floor << "\n";
                Rcpp::Rcout << "SD ceiling: " << mcmc_sd_ceil << "\n";
                if (save_samples) {
                    Rcpp::Rcout << "Saving samples: every " << sample_interval 
                               << " iterations after burn-in\n";
                }
            }
            
            my_chromosome.run_mcmc(
                mcmc_iterations,
                mcmc_burn_in,
                mcmc_initial_sd,
                mcmc_sd_floor,
                mcmc_sd_ceil,
                save_samples,
                sample_interval
            );
            
            // --- 5. Results saved by Chromosome class ---
            if (verbose) {
                Rcpp::Rcout << "\n=== Analysis Complete ===\n";
                Rcpp::Rcout << "All results saved: \n";
                Rcpp::Rcout << "  - final_positions.txt\n";
                Rcpp::Rcout << "  - initial_positions.txt\n";
                Rcpp::Rcout << "  - log_likelihood_trace.txt\n";
                Rcpp::Rcout << "  - block_timings.txt\n";
                Rcpp::Rcout << "  - mcmc_log.txt\n";
                Rcpp::Rcout << "DEBUG: About to exit scope block (destructor will be called)\n";
                Rcpp::Rcout.flush();
            }
        }  // Chromosome object goes out of scope and destructs here
        
        // // Shutdown spdlog to ensure proper cleanup of all loggers
        // spdlog::shutdown();
        
        // if (verbose) {
        //     Rcpp::Rcout << "DEBUG: After spdlog::shutdown()\n";
        //     Rcpp::Rcout.flush();
        // }       
        // spdlog::drop_all(); 
        // spdlog::shutdown();

        return true;
        
    } catch (const std::exception& e) {
        Rcpp::Rcerr << "C++ Exception: " << e.what() << std::endl;
        Rcpp::stop("Error in hispa_analyze: " + std::string(e.what()));
    } catch (...) {
        Rcpp::Rcerr << "Unknown C++ exception caught" << std::endl;
        Rcpp::stop("Unknown error in hispa_analyze");
    }
}
