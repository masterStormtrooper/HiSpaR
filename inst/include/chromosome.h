// =========================================================================
// chromosome.h - Header file for the Chromosome class
// =========================================================================
#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath> // For sqrt
#include <algorithm> // For std::min and std::max
#include <chrono> // For timing
#include <filesystem> // For directory creation
#include <sstream> // For string streams
#include <iomanip> // For time formatting
#include <armadillo> // Include the Armadillo library header
#include <omp.h> // For OpenMP
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

// Platform-specific includes for memory usage
#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>
#endif


/*
 * @struct GammaPrior
 * @brief Represents a Gamma(shape, rate) prior distribution for a specific distance type.
 * 
 * Parameterization: Gamma(shape, rate) where mean = shape/rate and variance = shape/rate^2
 * The log prior density is: (shape - 1) * log(x) - rate * x + shape * log(rate) - lgamma(shape)
 */
struct GammaPrior {
    double shape;  // Shape parameter (alpha)
    double rate;   // Rate parameter (beta)
    
    // Default constructor (uninformative prior)
    GammaPrior() : shape(0.001), rate(0.001) {}
    
    // Parameterized constructor
    GammaPrior(double shape_, double rate_) : shape(shape_), rate(rate_) {}
    
    /*
     * @brief Calculate the log prior density for a given distance value.
     * @param x The distance to evaluate (must be positive).
     * @return The log prior density, or -inf if x <= 0.
     */
    inline double log_density(double x) const {
        if (x <= 0) {
            return -std::numeric_limits<double>::infinity();
        }
        return (shape - 1.0) * std::log(x) - rate * x + shape * std::log(rate) - std::lgamma(shape);
    }
    
    /*
     * @brief Get the mean of the Gamma distribution.
     */
    inline double mean() const {
        return shape / rate;
    }
    
    /*
     * @brief Get the variance of the Gamma distribution.
     */
    inline double variance() const {
        return shape / (rate * rate);
    }
};


/*
 * @struct DistancePriors
 * @brief Contains gamma priors for different genomic distance separations between loci.
 * 
 * This structure holds gamma priors indexed by genomic distance (in number of loci).
 * For example:
 *   - distance 1: priors[1] represents adjacent loci (i, i+1)
 *   - distance 2: priors[2] represents loci separated by 1 (i, i+2)
 *   - distance k: priors[k] represents loci separated by k-1 (i, i+k)
 */
struct DistancePriors {
    std::map<int, GammaPrior> priors;  // Map from genomic distance to gamma prior
    std::vector<int> active_distances;  // Sorted list of distances with fitted priors
    
    // Default constructor (no priors)
    DistancePriors() {}
    
    /*
     * @brief Constructor that initializes empty priors for specified distances.
     * @param distances Vector of genomic distances (e.g., {1, 2, 5, 10}).
     */
    DistancePriors(const std::vector<int>& distances) {
        for (int d : distances) {
            if (d > 0) {
                priors[d] = GammaPrior();  // Initialize with uninformative prior
                active_distances.push_back(d);
            }
        }
        std::sort(active_distances.begin(), active_distances.end());
    }
    
    /*
     * @brief Check if a prior exists for a given distance.
     * @param distance The genomic distance to check.
     * @return True if a prior exists for this distance.
     */
    inline bool has_prior(int distance) const {
        return priors.find(distance) != priors.end();
    }
    
    /*
     * @brief Get the prior for a specific distance.
     * @param distance The genomic distance.
     * @return Reference to the GammaPrior for this distance.
     */
    inline const GammaPrior& get_prior(int distance) const {
        auto it = priors.find(distance);
        if (it != priors.end()) {
            return it->second;
        }
        static GammaPrior default_prior;
        return default_prior;
    }
    
    /*
     * @brief Set the prior for a specific distance.
     * @param distance The genomic distance.
     * @param prior The GammaPrior to set.
     */
    inline void set_prior(int distance, const GammaPrior& prior) {
        priors[distance] = prior;
        if (std::find(active_distances.begin(), active_distances.end(), distance) == active_distances.end()) {
            active_distances.push_back(distance);
            std::sort(active_distances.begin(), active_distances.end());
        }
    }
    
    /*
     * @brief Calculate the log prior for all relevant distances involving locus j.
     * @param positions The full position matrix (n x d).
     * @param pos_j The position of locus j (1 x d row vector).
     * @param j The index of the locus (0 to n-1).
     * @return Sum of log prior densities for all active distances involving j.
     */
    inline double log_prior_for_locus(const arma::mat& positions, const arma::rowvec& pos_j, arma::uword j) const {
        arma::uword n = positions.n_rows;
        double log_prior_sum = 0.0;
        
        for (int dist : active_distances) {
            const GammaPrior& prior = priors.at(dist);
            
            // Skip if this is an uninformative prior
            if (prior.shape <= 0.01 && prior.rate <= 0.01) {
                continue;
            }
            
            // Check backward distance (j-dist to j)
            if (j >= static_cast<arma::uword>(dist)) {
                double d = arma::norm(pos_j - positions.row(j - dist), 2);
                log_prior_sum += prior.log_density(d);
            }
            
            // Check forward distance (j to j+dist)
            if (j + dist < n) {
                double d = arma::norm(positions.row(j + dist) - pos_j, 2);
                log_prior_sum += prior.log_density(d);
            }
        }
        
        return log_prior_sum;
    }
    
    /*
     * @brief Check if any informative priors exist.
     * @return True if at least one prior is informative (shape > 0.01 and rate > 0.01).
     */
    inline bool has_informative_priors() const {
        for (const auto& pair : priors) {
            if (pair.second.shape > 0.01 && pair.second.rate > 0.01) {
                return true;
            }
        }
        return false;
    }
};


/*
 * @struct ClusterAdjacency
 * @brief Stores the indices of loci belonging to a cluster, its neighbors,
 * and its non-neighboring "stranger" clusters.
 */
struct ClusterAdjacency {
    arma::uvec self_indices;
    arma::uvec neighbor_indices;
    arma::uvec stranger_indices;
};


/*
 * @brief Convolve contact matrix by counting nonzero elements in upper triangular part of window.
 * Only processes upper triangular elements (excluding diagonal), counts nonzero elements
 * in the upper triangular portion of the window, then makes result symmetric.
 * @param contact_matrix Input contact matrix (arma::mat).
 * @param half_k Half window size (default=3).
 * @return Convoluted matrix with same shape as input (symmetric).
 */
arma::mat convolute_contacts(const arma::mat& contact_matrix, int half_k = 3);


/*
 * @class Chromosome
 * @brief Represents a chromosome and its associated Hi-C contact data.
 */
class Chromosome {
private:
    std::string chromosome_name;
    arma::mat contact_matrix;
    arma::uvec cluster_labels;
    arma::mat backbone_contact_matrix;
    std::vector<ClusterAdjacency> cluster_adjacencies;
    arma::mat position_matrix;
    arma::mat cluster_center_position_matrix;
    arma::mat pairwise_distance_matrix;
    arma::mat theta_matrix;  // Precomputed theta = exp(beta0 + beta1 * log(distance)) for convolution
    arma::mat convoluted_contact_matrix;
    double beta0;
    double beta1;
    
    // Prior distributions for different genomic distances between loci
    DistancePriors distance_priors;
    arma::mat prior_position_matrix;  // Positions used to fit the gamma prior
    arma::mat prior_contact_matrix;   // Contact matrix from prior data
    
    // MCMC options
    bool skip_zero_contact_loci;  // Whether to skip loci with no contacts during sampling
    bool sample_from_prior;  // Whether to sample new positions from prior positions

    // MCMC sample storage and tracking
    arma::vec mcmc_trace_beta0;
    arma::vec mcmc_trace_beta1;
    std::vector<arma::mat> mcmc_trace_cluster_centers;
    arma::vec mcmc_trace_log_likelihood;
    std::vector<double> mcmc_trace_block_durations;
    arma::mat best_position_matrix;
    double max_log_likelihood;
    
    // Initialization results
    std::vector<arma::mat> initial_cluster_structures;
    arma::mat backbone_structure;
    std::vector<double> initial_beta0s;
    std::vector<double> initial_beta1s;
    double backbone_beta0;
    double backbone_beta1;


    // Logger
    // std::shared_ptr<spdlog::logger> logger;

    // Private helper for logging memory
    void log_memory_usage();
    
    // Private MCMC sampling helpers
    void sample_beta_parameters(double& current_ll, double sd_b0, double sd_b1, 
                                int& accepted_b0, int& accepted_b1);
    void sample_cluster_centers(double& current_ll, double sd_center, int& accepted_center);
    void sample_locus_positions(double& current_ll, double sd_locus, int& accepted_locus);
    void sample_locus_positions_convoluted(double& current_ll, double sd_locus, int& accepted_locus, int k = 3);
    
    // Private helper for finding rotation matrix
    arma::mat rotation_solver(arma::vec sv, arma::vec dv) const;
    // Private helper to rotate a set of positions around their centroid
    // positions: n x d matrix (each row is a point)
    // rotation: d x d rotation matrix
    arma::mat rotate_positions(const arma::mat& positions, const arma::mat& rotation) const;
    
    // Private helper for likelihood calculation with explicit betas
    // double calculate_log_likelihood(const arma::mat& distances, const arma::mat& contacts, double b0, double b1) const;
    /*
     * @brief Optimized log-likelihood calculation using a direct memory loop.
     * @tparam T1 Matrix/subview type for distances.
     * @tparam T2 Matrix/subview type for contacts.
     * @param distances_expr An Armadillo object (e.g., mat, subview) of distances.
     * @param contacts_expr An Armadillo object of contacts.
     * @param b0 The beta0 parameter (intercept).
     * @param b1 The beta1 parameter (slope).
     * @return The calculated total log-likelihood.
     */
    template<typename T1, typename T2>
    double calculate_log_likelihood(
        const arma::Base<double, T1>& distances_expr, 
        const arma::Base<double, T2>& contacts_expr, 
        double b0, double b1) const 
    {
        const T1& distances = distances_expr.get_ref();
        const T2& contacts = contacts_expr.get_ref();
        // CRITICAL: Prevent out-of-bounds if matrices are misaligned
        if (distances.n_elem != contacts.n_elem) {
            std::cerr << "Error: Distance and contact matrices have different number of elements." << std::endl;
            return -arma::datum::inf;
        }
        double total_log_likelihood = 0.0;
        
        // const double* dist_mem = distances.memptr();
        // const double* cont_mem = contacts.memptr();
        // const arma::uword n_elem = distances.n_elem;

        // for (arma::uword i = 0; i < n_elem; ++i) {
        //     const double d = dist_mem[i];
            
        //     if (d <= 0) {
        //         continue;
        //     }

        //     const double c = cont_mem[i];
        //     const double log_lam = b0 + b1 * std::log(d);
        //     const double term = c * log_lam - std::exp(log_lam);
            
        //     if (std::isfinite(term)) {
        //         total_log_likelihood += term;
        //     }
        // }

        for (arma::uword i = 0; i < distances.n_elem; ++i) {
            const double d = distances(i);
            
            // Prevent log(0) or negative distances
            if (d <= 1e-12) continue;

            const double c = contacts(i);
            const double log_lam = b0 + b1 * std::log(d);
            const double term = c * log_lam - std::exp(log_lam);
        
            if (std::isfinite(term)) {
                total_log_likelihood += term;
            }
        }
        return total_log_likelihood;
    }


public:
    Chromosome(const std::string& name);
    // Constructor for sub-problems
    Chromosome(const std::string& name, const arma::mat& contacts);
    ~Chromosome();

    bool load_data_from_file(const std::string& filename);
    void assign_clusters();
    void assign_clusters(int num_clusters);
    bool assign_clusters(const arma::uvec& custom_labels);
    void calculate_backbone_contact_matrix(const std::vector<arma::uvec>& indices_by_cluster);
    void build_cluster_relationships(double k = 0.25);
    void build_cluster_relationships_by_distance(int window_size = 2);
    void initialize_positions();
    void initialize_structure_from_clusters(int sub_iterations, int sub_burn_in, double initial_sd);
    void scale_structure();
    void assemble_global_structure();
    
    /*
     * @brief Fit gamma priors for multiple genomic distances from a position matrix.
     * @param positions_file Path to the position matrix file (n x d, ASCII format).
     * @param distances Vector of genomic distances to fit (e.g., {1, 2, 5, 10} for adjacent, skip-1, skip-4, skip-9).
     * @return True if successfully fitted and stored the priors, false otherwise.
     */
    bool fit_distance_priors_from_file(const std::string& positions_file, const std::vector<int>& distances);
    
    /*
     * @brief Fit gamma priors for multiple genomic distances from a position matrix.
     * @param positions Position matrix (n x d).
     * @param distances Vector of genomic distances to fit (e.g., {1, 2, 5, 10}).
     * @return True if successfully fitted and stored the priors, false otherwise.
     */
    bool fit_distance_priors(const arma::mat& positions, const std::vector<int>& distances);
    
    // Legacy method for backward compatibility (fits only adjacent distances)
    bool fit_adjacent_distance_prior_from_file(const std::string& positions_file);
    bool fit_adjacent_distance_prior(const arma::mat& positions);
    
    arma::mat calculate_pairwise_distances(const arma::mat& pos1, const arma::mat& pos2) const;
    // double calculate_log_likelihood(const arma::mat& distances, const arma::mat& contacts) const;
    double calculate_log_likelihood(const arma::mat& distances, const arma::mat& contacts) const {
        return calculate_log_likelihood(distances, contacts, this->beta0, this->beta1);
    }

    /*
     * @brief Runs an MCMC simulation with adaptive proposal distributions.
     * @param iterations Total number of MCMC iterations.
     * @param burn_in Number of initial iterations to discard.
     * @param initial_sd The starting standard deviation for proposals.
     * @param sd_floor The minimum allowed standard deviation.
     * @param sd_ceiling The maximum allowed standard deviation.
     * @param save_samples Whether to save samples during MCMC.
     * @param sample_interval Save every kth sample after burn-in (default: 5).
     */
    void run_mcmc(int iterations, int burn_in, double initial_sd = 0.1, double sd_floor = 0.001, double sd_ceiling = 0.3, bool save_samples = false, int sample_interval = 5);
    void run_mcmc_convoluted(int iterations, int burn_in, double initial_sd = 0.1, double sd_floor = 0.001, double sd_ceiling = 0.3, bool save_samples = false, int sample_interval = 5, int k = 3);
    void run_mcmc_finetune(int iterations);
    // --- GETTERS ---
    const std::string& get_name() const;
    const arma::mat& get_contact_matrix() const;
    const arma::uvec& get_cluster_labels() const;
    const std::vector<ClusterAdjacency>& get_cluster_adjacencies() const;
    const arma::mat& get_backbone_matrix() const;
    const arma::mat& get_position_matrix() const;
    const arma::mat& get_cluster_center_position_matrix() const;
    const arma::mat& get_pairwise_distance_matrix() const;
    double get_beta0() const;
    double get_beta1() const;
    const arma::vec& get_mcmc_trace_beta0() const;
    const arma::vec& get_mcmc_trace_beta1() const;
    const std::vector<arma::mat>& get_mcmc_trace_cluster_centers() const;
    const arma::mat& get_best_position_matrix() const;
    const arma::vec& get_mcmc_trace_log_likelihood() const;
    const std::vector<double>& get_mcmc_trace_block_durations() const;
    const std::vector<arma::mat>& get_initial_cluster_structures() const;
    const arma::mat& get_backbone_structure() const;
    const std::vector<double>& get_initial_beta0s() const;
    const std::vector<double>& get_initial_beta1s() const;
    double get_backbone_beta0() const;
    double get_backbone_beta1() const;
    const DistancePriors& get_distance_priors() const;
    const GammaPrior& get_adjacent_distance_prior() const;  // Legacy accessor
    const arma::mat& get_prior_position_matrix() const;
    const arma::mat& get_prior_contact_matrix() const;
    bool get_skip_zero_contact_loci() const;
    bool get_sample_from_prior() const;
    const arma::mat& get_convoluted_contact_matrix() const;


    // --- SETTERS ---
    void set_beta0(double val);
    void set_beta1(double val);
    void set_distance_priors(const DistancePriors& priors);
    void set_distance_prior(int distance, double shape, double rate);
    void set_distance_prior(int distance, const GammaPrior& prior);
    void set_adjacent_distance_prior(double shape, double rate);  // Legacy setter
    void set_adjacent_distance_prior(const GammaPrior& prior);    // Legacy setter
    void set_position_matrix(const arma::mat& positions);
    void set_prior_contact_matrix(const arma::mat& contacts);
    bool load_prior_contact_matrix_from_file(const std::string& filename);
    void set_skip_zero_contact_loci(bool skip);
    void set_sample_from_prior(bool sample);
    void set_convoluted_contact_matrix(const arma::mat& convoluted_contacts);
    void compute_convoluted_contact_matrix(int k = 3);
};

#endif // CHROMOSOME_H
