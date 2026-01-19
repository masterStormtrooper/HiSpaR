// =========================================================================
// chromosome.cpp - Implementation file for the Chromosome class
// =========================================================================
#include "chromosome.h"
#include <set>
#include <Rcpp.h>


// =========================================================================
// Helper Function: convolute_contacts
// =========================================================================
/*
 * @brief Convolve contact matrix by counting nonzero elements in upper triangular part of window.
 * Only processes upper triangular elements (excluding diagonal), counts nonzero elements
 * in the upper triangular portion of the window, then makes result symmetric.
 * @param contact_matrix Input contact matrix (arma::mat).
 * @param half_k Half window size (default=3).
 * @return Convoluted matrix with same shape as input (symmetric).
 */
arma::mat convolute_contacts(const arma::mat& contact_matrix, int half_k) {
    arma::uword n = contact_matrix.n_rows;
    arma::mat result = arma::zeros<arma::mat>(n, n);
    
    // Only process upper triangular part (excluding diagonal)
    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = i + 1; j < n; j++) {
            // Define window boundaries with edge handling
            arma::uword i_start = (i >= static_cast<arma::uword>(half_k)) ? i - half_k : 0;
            arma::uword i_end = std::min(n, i + half_k + 1);
            arma::uword j_start = (j >= static_cast<arma::uword>(half_k)) ? j - half_k : 0;
            arma::uword j_end = std::min(n, j + half_k + 1);
            
            // Extract window
            arma::mat window = contact_matrix.submat(i_start, j_start, i_end - 1, j_end - 1);
            
            // Count only nonzero elements in the upper triangular part of the window
            int count = 0;
            for (arma::uword row = 0; row < window.n_rows; row++) {
                for (arma::uword col = 0; col < window.n_cols; col++) {
                    arma::uword global_row = i_start + row;
                    arma::uword global_col = j_start + col;
                    
                    // Only count if in upper triangular (global_col > global_row) and nonzero
                    if (global_col > global_row && window(row, col) != 0) {
                        count++;
                    }
                }
            }
            
            result(i, j) = count;
        }
    }
    
    // Make symmetric
    result = result + result.t();
    
    return result;
}

Chromosome::Chromosome(const std::string& name) : chromosome_name(name), beta0(0.0), beta1(0.0), skip_zero_contact_loci(true), sample_from_prior(false), max_log_likelihood(-arma::datum::inf) {
    std::cout << "Chromosome object '" << chromosome_name << "' created." << std::endl;
    try {
        std::filesystem::create_directory(chromosome_name);
        // std::string log_path = chromosome_name + "/mcmc_log.txt";
        // std::string logger_name = "logger_" + name;
        // std::replace(logger_name.begin(), logger_name.end(), '/', '_');
        // logger = spdlog::basic_logger_mt(logger_name, log_path, true); 
    } catch (const spdlog::spdlog_ex &ex) {
        std::cerr << "Log initialization failed: " << ex.what() << std::endl;
    }
}

Chromosome::Chromosome(const std::string& name, const arma::mat& contacts) : chromosome_name(name), contact_matrix(contacts), beta0(0.0), beta1(0.0), skip_zero_contact_loci(true), sample_from_prior(false), max_log_likelihood(-arma::datum::inf) {
    std::cout << "Chromosome object '" << chromosome_name << "' created for sub-problem." << std::endl;
    try {
        std::filesystem::create_directories(chromosome_name);
        // std::string log_path = chromosome_name + "/mcmc_log.txt";
        // std::string logger_name = "logger_" + name;
        // std::replace(logger_name.begin(), logger_name.end(), '/', '_');
        // logger = spdlog::basic_logger_mt(logger_name, log_path, true);
    } catch (const spdlog::spdlog_ex &ex) {
        std::cerr << "Log initialization failed for " << name << ": " << ex.what() << std::endl;
    }
    // For a sub-problem, all loci belong to a single cluster (label 0)
    arma::uvec single_cluster_labels(contact_matrix.n_rows, arma::fill::zeros);
    assign_clusters(single_cluster_labels);
}


Chromosome::~Chromosome() {
//    if (logger) {
//         std::string name = // logger->name();
//         // 1. Flush the logs
//         // logger->flush();
//         // 2. Clear the shared_ptr reference held by this class
//         logger.reset(); 
//         // 3. REMOVE from the global spdlog registry
//         spdlog::drop(name); 
//     }
    Rcpp::Rcout << "Chromosome object destroyed." << std::endl;
}

bool Chromosome::load_data_from_file(const std::string& filename) {
    Rcpp::Rcout << "Loading data from " << filename << std::endl;
    bool success = contact_matrix.load(filename, arma::raw_ascii);

    if (success) {
        Rcpp::Rcout << "Successfully loaded matrix with size: "
                  << contact_matrix.n_rows << " x " << contact_matrix.n_cols
                  << std::endl;
        cluster_labels.reset();
        cluster_adjacencies.clear();
        backbone_contact_matrix.reset();
        position_matrix.reset();
        cluster_center_position_matrix.reset();
        pairwise_distance_matrix.reset();
        max_log_likelihood = -arma::datum::inf;
    } else {
        Rcpp::Rcerr << "Error: Failed to load data from '" << filename << "'." << std::endl;
    }
    return success;
}

void Chromosome::assign_clusters() {
    if (contact_matrix.is_empty()) {
        Rcpp::Rcerr << "Error: Cannot assign clusters. Load data first." << std::endl;
        return;
    }
    arma::uword n_loci = contact_matrix.n_rows;
    int num_clusters = static_cast<int>(sqrt(n_loci));
    if (num_clusters == 0) num_clusters = 1;
    
    assign_clusters(num_clusters);
}

void Chromosome::assign_clusters(int num_clusters) {
    if (contact_matrix.is_empty()) {
        Rcpp::Rcerr << "Error: Cannot assign clusters. Load data first." << std::endl;
        return;
    }
    if (num_clusters <= 0) {
        Rcpp::Rcerr << "Error: Number of clusters must be positive." << std::endl;
        return;
    }
    arma::uword n_loci = contact_matrix.n_rows;
    arma::uword cluster_size = n_loci / num_clusters;
    cluster_labels.set_size(n_loci);
    Rcpp::Rcout << "Assigning " << n_loci << " loci into " << num_clusters << " clusters..." << std::endl;
    for (arma::uword i = 0; i < n_loci; ++i) {
        arma::uword label = i / cluster_size;
        if (label >= num_clusters) {
            label = num_clusters - 1;
        }
        cluster_labels(i) = label;
    }
    cluster_adjacencies.clear();
    backbone_contact_matrix.reset();
}

bool Chromosome::assign_clusters(const arma::uvec& custom_labels) {
    if (contact_matrix.is_empty()) {
        Rcpp::Rcerr << "Error: Cannot assign clusters. Load data first." << std::endl;
        return false;
    }
    if (custom_labels.n_elem != contact_matrix.n_rows) {
        Rcpp::Rcerr << "Error: Size of custom label vector (" << custom_labels.n_elem
                  << ") does not match the number of loci (" << contact_matrix.n_rows
                  << ")." << std::endl;
        return false;
    }
    Rcpp::Rcout << "Assigning custom cluster labels..." << std::endl;
    this->cluster_labels = custom_labels;
    cluster_adjacencies.clear();
    backbone_contact_matrix.reset();
    return true;
}

void Chromosome::calculate_backbone_contact_matrix(const std::vector<arma::uvec>& indices_by_cluster) {
    arma::uword num_clusters = indices_by_cluster.size();
    
    // Determine the contact source: prior contact matrix, calculated from prior positions, or current contact matrix
    arma::mat contact_source = contact_matrix;
    if (!prior_contact_matrix.is_empty()) {
        Rcpp::Rcout << "Using prior contact matrix for backbone calculation..." << std::endl;
        // logger->info("Using prior contact matrix ({} x {}) for backbone", prior_contact_matrix.n_rows, prior_contact_matrix.n_cols);
        contact_source = prior_contact_matrix;
    } 
    // Calculate backbone matrix using geometric mean
    backbone_contact_matrix.set_size(num_clusters, num_clusters);
    for (arma::uword c1 = 0; c1 < num_clusters; ++c1) {
        for (arma::uword c2 = c1; c2 < num_clusters; ++c2) {
            arma::mat sub_contact_matrix = contact_source.submat(indices_by_cluster[c1], indices_by_cluster[c2]);
            
            // Calculate geometric mean with +1/-1 for zeros
            arma::mat sub_contact_matrix_plus_one = sub_contact_matrix + 1.0;
            double log_sum = arma::accu(arma::log(sub_contact_matrix_plus_one));
            double geo_mean = std::exp(log_sum / sub_contact_matrix.n_elem) - 1.0;
            
            backbone_contact_matrix(c1, c2) = geo_mean;
            backbone_contact_matrix(c2, c1) = geo_mean;
        }
    }
}

void Chromosome::build_cluster_relationships(double k) {
    if (cluster_labels.is_empty()) {
        Rcpp::Rcerr << "Error: Cannot build relationships. Assign clusters first." << std::endl;
        return;
    }
    if (k < 0.0 || k > 1.0) {
        Rcpp::Rcerr << "Warning: k should be between 0.0 and 1.0." << std::endl;
    }
    arma::uvec unique_labels = arma::unique(cluster_labels);
    arma::uword num_clusters = unique_labels.n_elem;
    Rcpp::Rcout << "Building relationships for " << num_clusters << " clusters using geometric mean..." << std::endl;
    std::vector<arma::uvec> indices_by_cluster(num_clusters);
    for (arma::uword i = 0; i < cluster_labels.n_elem; ++i) {
        indices_by_cluster[cluster_labels(i)].insert_rows(0, arma::uvec{i});
    }
    
    // Calculate backbone contact matrix
    calculate_backbone_contact_matrix(indices_by_cluster);
    
    cluster_adjacencies.resize(num_clusters);
    arma::uword num_neighbors = static_cast<arma::uword>(std::round(k * (num_clusters > 1 ? num_clusters - 1 : 0)));
    if (num_clusters > 2 && num_neighbors < 2) {
        num_neighbors = 2;
    }

    for (arma::uword c = 0; c < num_clusters; ++c) {
        cluster_adjacencies[c].self_indices = indices_by_cluster[c];

        std::vector<std::pair<double, arma::uword>> neighbor_candidates;
        for (arma::uword other_c = 0; other_c < num_clusters; ++other_c) {
            if (c != other_c) {
                neighbor_candidates.push_back({backbone_contact_matrix(c, other_c), other_c});
            }
        }
        
        std::sort(neighbor_candidates.rbegin(), neighbor_candidates.rend());

        std::vector<std::pair<double, arma::uword>> top_neighbors;
        if (num_neighbors > 0 && !neighbor_candidates.empty()) {
            top_neighbors.assign(neighbor_candidates.begin(), neighbor_candidates.begin() + std::min((size_t)num_neighbors, neighbor_candidates.size()));
        }

        auto is_neighbor = [&](arma::uword cluster_idx) {
            for(const auto& p : top_neighbors) if(p.second == cluster_idx) return true;
            return false;
        };

        if (c > 0 && !is_neighbor(c - 1)) {
            if (!top_neighbors.empty()) top_neighbors.pop_back();
            top_neighbors.push_back({backbone_contact_matrix(c, c - 1), c - 1});
        }
        if (c < num_clusters - 1 && !is_neighbor(c + 1)) {
            if (!top_neighbors.empty()) {
                 // Sort again to find the lowest to drop
                std::sort(top_neighbors.begin(), top_neighbors.end()); 
                // Ensure we don't drop the one we just added if it's the lowest
                if (top_neighbors[0].second == c-1 && top_neighbors.size() > 1) {
                    top_neighbors.erase(top_neighbors.begin()+1);
                } else {
                    top_neighbors.erase(top_neighbors.begin());
                }
            }
            top_neighbors.push_back({backbone_contact_matrix(c, c + 1), c + 1});
        }

        std::set<arma::uword> final_neighbor_labels;
        for(const auto& p : top_neighbors) {
            final_neighbor_labels.insert(p.second);
        }

        std::vector<arma::uword> temp_neighbor_indices, temp_stranger_indices;
        for (arma::uword other_c = 0; other_c < num_clusters; ++other_c) {
            if (c == other_c) continue;
            const arma::uvec& other_indices = indices_by_cluster[other_c];
            if (final_neighbor_labels.count(other_c)) {
                temp_neighbor_indices.insert(temp_neighbor_indices.end(), other_indices.begin(), other_indices.end());
            } else {
                temp_stranger_indices.insert(temp_stranger_indices.end(), other_indices.begin(), other_indices.end());
            }
        }
        cluster_adjacencies[c].neighbor_indices = arma::conv_to<arma::uvec>::from(temp_neighbor_indices);
        cluster_adjacencies[c].stranger_indices = arma::conv_to<arma::uvec>::from(temp_stranger_indices);
    }

    // --- Log the relationships ---
    // logger->info("--- Cluster Adjacency Report (Statistics-based) ---");
    for (arma::uword c = 0; c < num_clusters; ++c) {
        const auto& adj = cluster_adjacencies[c];
        std::string neighbors_str = "None";
        if (!adj.neighbor_indices.is_empty()) {
            arma::uvec neighbor_labels = arma::unique(cluster_labels(adj.neighbor_indices));
            std::stringstream ss;
            neighbor_labels.t().print(ss);
            neighbors_str = ss.str();
            neighbors_str.erase(std::remove(neighbors_str.begin(), neighbors_str.end(), '\n'), neighbors_str.end());
        }
        std::string strangers_str = "None";
        if (!adj.stranger_indices.is_empty()) {
            arma::uvec stranger_labels = arma::unique(cluster_labels(adj.stranger_indices));
            std::stringstream ss;
            stranger_labels.t().print(ss);
            strangers_str = ss.str();
            strangers_str.erase(std::remove(strangers_str.begin(), strangers_str.end(), '\n'), strangers_str.end());
        }
        // logger->info("Cluster {}: Neighbors [{}], Strangers [{}]", c, neighbors_str, strangers_str);
    }
}

void Chromosome::build_cluster_relationships_by_distance(int window_size) {
    if (cluster_labels.is_empty()) {
        Rcpp::Rcerr << "Error: Cannot build relationships. Assign clusters first." << std::endl;
        return;
    }
    arma::uvec unique_labels = arma::unique(cluster_labels);
    arma::uword num_clusters = unique_labels.n_elem;
    Rcpp::Rcout << "Building relationships for " << num_clusters << " clusters using distance window..." << std::endl;

    std::vector<arma::uvec> indices_by_cluster(num_clusters);
    for (arma::uword i = 0; i < cluster_labels.n_elem; ++i) {
        indices_by_cluster[cluster_labels(i)].insert_rows(0, arma::uvec{i});
    }

    // Calculate backbone contact matrix
    calculate_backbone_contact_matrix(indices_by_cluster);

    cluster_adjacencies.resize(num_clusters);
    for (int c = 0; c < num_clusters; ++c) {
        cluster_adjacencies[c].self_indices = indices_by_cluster[c];

        std::vector<arma::uword> temp_neighbor_indices, temp_stranger_indices;
        for (int other_c = 0; other_c < num_clusters; ++other_c) {
            if (c == other_c) continue;

            if (std::abs(c - other_c) <= window_size) {
                const arma::uvec& other_indices = indices_by_cluster[other_c];
                temp_neighbor_indices.insert(temp_neighbor_indices.end(), other_indices.begin(), other_indices.end());
            } else {
                const arma::uvec& other_indices = indices_by_cluster[other_c];
                temp_stranger_indices.insert(temp_stranger_indices.end(), other_indices.begin(), other_indices.end());
            }
        }
        cluster_adjacencies[c].neighbor_indices = arma::conv_to<arma::uvec>::from(temp_neighbor_indices);
        cluster_adjacencies[c].stranger_indices = arma::conv_to<arma::uvec>::from(temp_stranger_indices);
    }

    // --- Log the relationships ---
    // logger->info("--- Cluster Adjacency Report (Distance-based, window={}) ---", window_size);
    for (arma::uword c = 0; c < num_clusters; ++c) {
        const auto& adj = cluster_adjacencies[c];
        std::string neighbors_str = "None";
        if (!adj.neighbor_indices.is_empty()) {
            arma::uvec neighbor_labels = arma::unique(cluster_labels(adj.neighbor_indices));
            std::stringstream ss;
            neighbor_labels.t().print(ss);
            neighbors_str = ss.str();
            neighbors_str.erase(std::remove(neighbors_str.begin(), neighbors_str.end(), '\n'), neighbors_str.end());
        }
        std::string strangers_str = "None";
        if (!adj.stranger_indices.is_empty()) {
            arma::uvec stranger_labels = arma::unique(cluster_labels(adj.stranger_indices));
            std::stringstream ss;
            stranger_labels.t().print(ss);
            strangers_str = ss.str();
            strangers_str.erase(std::remove(strangers_str.begin(), strangers_str.end(), '\n'), strangers_str.end());
        }
        // logger->info("Cluster {}: Neighbors [{}], Strangers [{}]", c, neighbors_str, strangers_str);
    }
}


void Chromosome::initialize_positions() {
    if (contact_matrix.is_empty() || cluster_labels.is_empty()) {
        Rcpp::Rcerr << "Error: Load data and assign clusters before initializing positions." << std::endl;
        return;
    }
    arma::uword n_loci = contact_matrix.n_rows;
    Rcpp::Rcout << "Initializing " << n_loci << " loci positions along the x-axis..." << std::endl;

    // Create an n_loci x 3 matrix filled with zeros
    position_matrix.zeros(n_loci, 3);

    // Set the x-coordinates (first column) to be evenly spaced by 1.0/500.0
    position_matrix.col(0) = arma::regspace<arma::vec>(0, n_loci - 1) * (1.0 / 500.0);
    
    // y and z coordinates (columns 1 and 2) are already zero.

    arma::uvec unique_labels = arma::unique(cluster_labels);
    arma::uword num_clusters = unique_labels.n_elem;

    cluster_center_position_matrix.set_size(num_clusters, 3);
    Rcpp::Rcout << "Calculating " << num_clusters << " cluster center positions..." << std::endl;
    for (arma::uword c = 0; c < num_clusters; ++c) {
        arma::uvec indices = arma::find(cluster_labels == c);
        cluster_center_position_matrix.row(c) = arma::mean(position_matrix.rows(indices), 0);
    }
    
    Rcpp::Rcout << "Calculating initial hybrid pairwise distance matrix..." << std::endl;
    pairwise_distance_matrix = calculate_pairwise_distances(position_matrix, position_matrix);

    for (arma::uword c1 = 0; c1 < num_clusters; ++c1) {
        const arma::uvec& stranger_idx_c1 = cluster_adjacencies[c1].stranger_indices;
        if (stranger_idx_c1.is_empty()) continue;
        
        arma::uvec unique_stranger_labels = arma::unique(cluster_labels(stranger_idx_c1));

        for (arma::uword c2_label : unique_stranger_labels) {
            if (c1 >= c2_label) continue;
            
            const arma::uvec& c1_indices = cluster_adjacencies[c1].self_indices;
            const arma::uvec& c2_indices = cluster_adjacencies[c2_label].self_indices;

            double dist_c1_c2 = arma::norm(cluster_center_position_matrix.row(c1) - cluster_center_position_matrix.row(c2_label));
            
            pairwise_distance_matrix.submat(c1_indices, c2_indices).fill(dist_c1_c2);
            pairwise_distance_matrix.submat(c2_indices, c1_indices).fill(dist_c1_c2);
        }
    }

    // Initialize best state
    max_log_likelihood = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    best_position_matrix = position_matrix;
    position_matrix.save(chromosome_name + "/initial_positions.txt", arma::raw_ascii);
    pairwise_distance_matrix.save(chromosome_name + "/initial_distances.txt", arma::raw_ascii);
    Rcpp::Rcout << "Initial log-likelihood: " << max_log_likelihood << std::endl;
}


void Chromosome::initialize_structure_from_clusters(int sub_iterations, int sub_burn_in, double initial_sd) {
    if (cluster_adjacencies.empty()) {
        Rcpp::Rcerr << "Error: Build cluster relationships before initializing from clusters." << std::endl;
        return;
    }

    arma::uword num_clusters = cluster_adjacencies.size();
    initial_cluster_structures.resize(num_clusters);
    initial_beta0s.resize(num_clusters);
    initial_beta1s.resize(num_clusters);

    // Create initialization directory
    std::string init_dir = chromosome_name + "/initialization";
    std::filesystem::create_directories(init_dir);

    Rcpp::Rcout << "\n--- Initializing Structures for Each Cluster and Backbone (in parallel) ---" << std::endl;
    auto init_start_time = std::chrono::steady_clock::now();

    #pragma omp parallel for
    for (arma::uword c = 0; c <= num_clusters; ++c) {
        if (c < num_clusters) {
            // --- Process Individual Clusters ---
            const arma::uvec& self_idx = cluster_adjacencies[c].self_indices;
            arma::mat sub_contacts = contact_matrix.submat(self_idx, self_idx);
            
            std::string cluster_name = init_dir + "/cluster_" + std::to_string(c);
            
            #pragma omp critical
            std::cout << "--- Processing " << cluster_name << " (" << self_idx.n_elem << " loci) on thread " << omp_get_thread_num() << " ---" << std::endl;

            Chromosome cluster_chr(cluster_name, sub_contacts);
            cluster_chr.assign_clusters(1);
            cluster_chr.build_cluster_relationships_by_distance(0); 
            cluster_chr.initialize_positions();
            cluster_chr.set_beta0(3.0);
            cluster_chr.set_beta1(-2.0);
            cluster_chr.run_mcmc(sub_iterations, sub_burn_in, initial_sd);

            initial_cluster_structures[c] = cluster_chr.get_best_position_matrix();
            initial_beta0s[c] = cluster_chr.get_beta0();
            initial_beta1s[c] = cluster_chr.get_beta1();
            
            arma::mat best_sub_pos = cluster_chr.get_best_position_matrix();
            best_sub_pos.save(cluster_name + "/final_positions.txt", arma::raw_ascii);
        } else {
            // --- Process Backbone Structure ---
            std::string backbone_name = init_dir + "/backbone";
            
            #pragma omp critical
            std::cout << "--- Processing " << backbone_name << " (" << this->backbone_contact_matrix.n_rows << " clusters) on thread " << omp_get_thread_num() << " ---" << std::endl;
            
            arma::mat backbone_contacts_double = this->backbone_contact_matrix;
            Chromosome backbone_chr(backbone_name, backbone_contacts_double);
            backbone_chr.build_cluster_relationships_by_distance(0);
            backbone_chr.initialize_positions();
            backbone_chr.set_beta0(1.0);
            backbone_chr.set_beta1(-1.0);
            backbone_chr.run_mcmc(sub_iterations, sub_burn_in, initial_sd);

            #pragma omp critical
            {
                this->backbone_structure = backbone_chr.get_best_position_matrix();
                this->backbone_beta0 = backbone_chr.get_beta0();
                this->backbone_beta1 = backbone_chr.get_beta1();
                this->backbone_structure.save(backbone_name + "/final_positions.txt", arma::raw_ascii);
            }
        }
    }
    std::cout << "--- Finished Initializing All Cluster and Backbone Structures ---" << std::endl;
    auto init_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> init_duration = init_end_time - init_start_time;
    // logger->info("Initialization completed in {:.2f} seconds.", init_duration.count());
}

arma::mat Chromosome::rotation_solver(arma::vec sv, arma::vec dv) const {
    arma::vec startvec = arma::normalise(sv);
    arma::vec destvec = arma::normalise(dv);
    arma::vec v = arma::cross(startvec, destvec);
    double s = arma::norm(v);
    double c = arma::dot(startvec, destvec);

    // If s is close to zero, vectors are collinear, no rotation needed.
    if (s < 1e-9) {
        return arma::eye<arma::mat>(3, 3);
    }

    arma::mat vx = {
        {0, -v(2), v(1)},
        {v(2), 0, -v(0)},
        {-v(1), v(0), 0}
    }; 
    arma::mat vxvx = vx * vx;
    arma::mat identity(3, 3, arma::fill::eye);
    double scaler = (1 - c) / (s * s);
    arma::mat r = identity + vx + vxvx * scaler;
    return r;
}

arma::mat Chromosome::rotate_positions(const arma::mat& positions, const arma::mat& rotation) const {
    // Validate input
    if (positions.is_empty()) {
        // logger->warn("rotate_positions called with empty positions matrix.");
        return positions;
    }
    if (rotation.is_empty()) {
        // logger->warn("rotate_positions called with empty rotation matrix. Returning original positions.");
        return positions;
    }

    // arma::uword n = positions.n_rows;
    arma::uword d = positions.n_cols;

    if (rotation.n_rows != d || rotation.n_cols != d) {
        // logger->warn("rotate_positions: rotation matrix dimension ({}, {}) does not match positions dimensionality ({}). Returning original positions.", rotation.n_rows, rotation.n_cols, d);
        return positions;
    }

    // Compute centroid (row vector)
    arma::rowvec centroid = arma::mean(positions, 0);

    // Convert to relative positions (n x d), apply rotation, and add centroid back
    arma::mat relative = positions.each_row() - centroid;
    arma::mat rotated_relative = relative * rotation.t();
    arma::mat result = rotated_relative;
    result.each_row() += centroid;

    return result;
}


void Chromosome::scale_structure() {
    arma::mat beta_mat(initial_beta0s.size(), 2);
    beta_mat.col(0) = arma::conv_to<arma::vec>::from(initial_beta0s);
    beta_mat.col(1) = arma::conv_to<arma::vec>::from(initial_beta1s);
    
    double median_beta0 = arma::median(beta_mat.col(0));
    double median_beta1 = arma::median(beta_mat.col(1));

    this->beta0 = median_beta0;
    this->beta1 = median_beta1;

    // logger->info("Median beta0 for assembly: {:.4f}", median_beta0);
    // logger->info("Median beta1 for assembly: {:.4f}", median_beta1);


    // now scale backbone structure
    arma::mat backbone_distances = calculate_pairwise_distances(this->backbone_structure, this->backbone_structure);
    std::vector<double> backbone_ratios_vec;
    backbone_ratios_vec.reserve(backbone_contact_matrix.n_elem);
    for(arma::uword i = 0; i < backbone_contact_matrix.n_rows; ++i) {
        for(arma::uword j = i + 1; j < backbone_contact_matrix.n_cols; ++j) {
            if (backbone_contact_matrix(i, j) > 0 && backbone_distances(i, j) > 0) {
                double expected_dist = exp((log(backbone_contact_matrix(i, j)) - median_beta0) / median_beta1);
                backbone_ratios_vec.push_back(expected_dist / backbone_distances(i, j));
            }
        }
    }
    arma::vec backbone_ratios(backbone_ratios_vec);
    double backbone_scale = backbone_ratios.is_empty() ? 1.0 : arma::median(backbone_ratios);
    backbone_scale = std::max(0.5, std::min(10.0, backbone_scale));
    // logger->info("Scaling backbone structure by a factor of {:.4f}", backbone_scale);
    arma::rowvec backbone_center = arma::mean(backbone_structure, 0);
    arma::mat relative_backbone_pos = backbone_structure.each_row() - backbone_center;
    this->backbone_structure = relative_backbone_pos * backbone_scale;
    this->backbone_structure.each_row() += backbone_center;
    

    arma::uword n_loci = contact_matrix.n_rows;
    position_matrix.set_size(n_loci, 3);
    arma::uword num_clusters = initial_cluster_structures.size();

    for (arma::uword c = 0; c < num_clusters; ++c) {
        arma::mat cluster_pos = initial_cluster_structures[c];
        const arma::uvec& cluster_indices = cluster_adjacencies[c].self_indices;
        const arma::mat& cluster_contacts = contact_matrix.submat(cluster_indices, cluster_indices);
    
        // --- Scaling ---
        arma::mat current_distances = calculate_pairwise_distances(cluster_pos, cluster_pos);
        std::vector<double> ratios_vec;
        ratios_vec.reserve(cluster_contacts.n_elem);
        for(arma::uword i = 0; i < cluster_contacts.n_rows; ++i) {
            for(arma::uword j = i + 1; j < cluster_contacts.n_cols; ++j) {
                if (cluster_contacts(i, j) > 0 && current_distances(i, j) > 0) {
                    double expected_dist = exp((log(cluster_contacts(i, j)) - median_beta0) / median_beta1);
                    ratios_vec.push_back(expected_dist / current_distances(i, j));
                }
            }
        }
        arma::vec ratios(ratios_vec);
        double scale = ratios.is_empty() ? 1.0 : arma::median(ratios);
        scale = std::max(0.5, std::min(2.0, scale));
        // logger->info("Scaling cluster {} by a factor of {:.4f}", c, scale);

        // calculate cluster center position
        arma::rowvec cluster_center = arma::mean(cluster_pos, 0);
        arma::mat relative_pos = cluster_pos.each_row() - cluster_center;
        arma::mat final_cluster_pos = relative_pos * scale;
        final_cluster_pos.each_row() += backbone_structure.row(c);
        position_matrix.rows(cluster_indices) = final_cluster_pos;
    }
    position_matrix.save(chromosome_name + "/post_scale_positions.txt", arma::raw_ascii);
}


void Chromosome::assemble_global_structure() {
    if (initial_cluster_structures.empty() || backbone_structure.is_empty()) {
        Rcpp::Rcerr << "Error: Must run initialize_structure_from_clusters() before assembling." << std::endl;
        return;
    }

    Rcpp::Rcout << "\n--- Assembling Global Structure ---" << std::endl;
    auto assembly_start_time = std::chrono::steady_clock::now();
    scale_structure();
    // arma::uword n_loci = contact_matrix.n_rows;
    arma::uword num_clusters = initial_cluster_structures.size();

    for (arma::uword c = 1; c < num_clusters; ++c) {
        // --- Rotation and Alignment (for c > 0) ---
        const arma::uvec& cluster_indices = cluster_adjacencies[c].self_indices;
        arma::rowvec prev_center = backbone_structure.row(c - 1);
        arma::rowvec current_center = backbone_structure.row(c);
        arma::uvec prev_cluster_indices = cluster_adjacencies[c - 1].self_indices;
        arma::uword prev_last_locus_idx = prev_cluster_indices.max();
        arma::uword curr_first_locus_idx = cluster_indices.min();
        arma::mat final_cluster_pos = position_matrix.rows(cluster_indices);

        arma::rowvec prev_last_locus_pos = position_matrix.row(prev_last_locus_idx);
        arma::rowvec current_first_locus_pos = final_cluster_pos.row(0);
        Rcpp::Rcout << "Aligning cluster " << c << " to previous cluster " << (c - 1) << "..." << std::endl;
        Rcpp::Rcout << "Previous last idx: " << prev_last_locus_idx << "prev cluster size: " << prev_cluster_indices.n_elem << std::endl;
        current_first_locus_pos.print("Current first locus position: ");
        prev_last_locus_pos.print("Previous last locus position: ");

        arma::vec start_vec = arma::trans(final_cluster_pos.row(0) - current_center); // Vector from origin to first point
        arma::vec dest_vec = arma::trans(prev_last_locus_pos - current_center); // Vector from origin to last point of previous cluster

        arma::mat R = rotation_solver(start_vec, dest_vec);
        final_cluster_pos = rotate_positions(final_cluster_pos, R);
        position_matrix.rows(cluster_indices) = final_cluster_pos;
        current_first_locus_pos = position_matrix.row(curr_first_locus_idx);
        prev_last_locus_pos = position_matrix.row(prev_last_locus_idx);

        arma::vec current_vec = current_first_locus_pos.t() - prev_last_locus_pos.t();
        double current_dist = arma::norm(current_vec);
        current_vec.print("Current vector: ");
        double contact = contact_matrix(prev_last_locus_idx, curr_first_locus_idx);
        double expected_dist = current_dist;
        if (contact > 0) {
            expected_dist = exp((log(contact) - this->beta0) / this->beta1);
        }

        double translation_distance = expected_dist - current_dist;
        Rcpp::Rcout << "Aligning cluster " << c << ": current_dist = " << current_dist
                  << ", expected_dist = " << expected_dist
                  << ", translation_distance = " << translation_distance << std::endl;
        for (arma::uword c1 = c; c1 < num_clusters; ++c1) {
            const arma::uvec& c1_indices = cluster_adjacencies[c1].self_indices;
            arma::mat c1_pos = position_matrix.rows(c1_indices);
            // --- Final Translation to match expected distance ---
            arma::rowvec new_cluster_center = backbone_structure.row(c1) - translation_distance * normalise(dest_vec).t(); // FLAG1
            arma::rowvec current_cluster_center = backbone_structure.row(c1);
            c1_pos.each_row() += (new_cluster_center - current_cluster_center);
            backbone_structure.row(c1) = new_cluster_center;
            position_matrix.rows(c1_indices) = c1_pos;

        }
        Rcpp::Rcout << "New previous locus pos " << position_matrix.row(prev_last_locus_idx) << " new current first pos " << position_matrix.row(curr_first_locus_idx) << std::endl;
        double new_current_dist = arma::norm(position_matrix.row(curr_first_locus_idx).t() - position_matrix.row(prev_last_locus_idx).t());
        Rcpp::Rcout << "Post-alignment distance between clusters " << (c - 1) << " and " << c << ": " << new_current_dist << std::endl;
    }

    Rcpp::Rcout << "--- Global Assembly Complete ---" << std::endl;  
     
    auto assembly_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> assembly_duration = assembly_end_time - assembly_start_time;
    // logger->info("Assembly completed in {:.2f} seconds.", assembly_duration.count());

    arma::uvec unique_labels = arma::unique(cluster_labels);

    cluster_center_position_matrix.set_size(num_clusters, 3);
    Rcpp::Rcout << "Calculating " << num_clusters << " cluster center positions..." << std::endl;
    for (arma::uword c = 0; c < num_clusters; ++c) {
        arma::uvec indices = arma::find(cluster_labels == c);
        cluster_center_position_matrix.row(c) = arma::mean(position_matrix.rows(indices), 0);
    }
    Rcpp::Rcout << "Calculating initial hybrid pairwise distance matrix..." << std::endl;
    pairwise_distance_matrix = calculate_pairwise_distances(position_matrix, position_matrix);

    for (arma::uword c1 = 0; c1 < num_clusters; ++c1) {
        const arma::uvec& stranger_idx_c1 = cluster_adjacencies[c1].stranger_indices;
        if (stranger_idx_c1.is_empty()) continue;
        
        arma::uvec unique_stranger_labels = arma::unique(cluster_labels(stranger_idx_c1));

        for (arma::uword c2_label : unique_stranger_labels) {
            if (c1 >= c2_label) continue;
            
            const arma::uvec& c1_indices = cluster_adjacencies[c1].self_indices;
            const arma::uvec& c2_indices = cluster_adjacencies[c2_label].self_indices;

            double dist_c1_c2 = arma::norm(cluster_center_position_matrix.row(c1) - cluster_center_position_matrix.row(c2_label));
            
            pairwise_distance_matrix.submat(c1_indices, c2_indices).fill(dist_c1_c2);
            pairwise_distance_matrix.submat(c2_indices, c1_indices).fill(dist_c1_c2);
        }
    }

    // Initialize best state
    max_log_likelihood = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    best_position_matrix = position_matrix;
    position_matrix.save(chromosome_name + "/initial_positions.txt", arma::raw_ascii);
    pairwise_distance_matrix.save(chromosome_name + "/initial_distances.txt", arma::raw_ascii);
    Rcpp::Rcout << "Initial log-likelihood: " << max_log_likelihood << std::endl;
    Rcpp::Rcout << "Position matrix has shape: " << position_matrix.n_rows << " x " << position_matrix.n_cols << std::endl;
    Rcpp::Rcout << "DEBUG: End of assemble_global_structure(), this=" << this 
              << ", position_matrix.is_empty()=" << position_matrix.is_empty() 
              << ", position_matrix.memptr()=" << (void*)position_matrix.memptr() << std::endl;
}


arma::mat Chromosome::calculate_pairwise_distances(const arma::mat& pos1, const arma::mat& pos2) const {
    arma::mat p1_2 = arma::sum(arma::square(pos1), 1);
    arma::mat p2_2 = arma::sum(arma::square(pos2), 1);
    arma::mat p1p2 = -2 * pos1 * pos2.t();
    p1p2.each_col() += p1_2;
    p1p2.each_row() += p2_2.t();
    return arma::sqrt(p1p2.clamp(0, p1p2.max()));
}


bool Chromosome::fit_distance_priors_from_file(const std::string& positions_file, const std::vector<int>& distances) {
    arma::mat positions;
    bool loaded = positions.load(positions_file, arma::raw_ascii);
    
    if (!loaded) {
        Rcpp::Rcerr << "Error: Could not load position matrix from " << positions_file << std::endl;
        // logger->error("Failed to load position matrix from {}", positions_file);
        return false;
    }
    
    return fit_distance_priors(positions, distances);
}


bool Chromosome::fit_distance_priors(const arma::mat& positions, const std::vector<int>& distances) {
    if (positions.n_rows < 2) {
        Rcpp::Rcerr << "Error: Need at least 2 positions to fit prior." << std::endl;
        // logger->error("Cannot fit prior: position matrix has fewer than 2 rows");
        return false;
    }
    
    // Initialize distance priors structure
    distance_priors = DistancePriors(distances);
    prior_position_matrix = positions;  // Store the positions used to fit the priors
    
    Rcpp::Rcout << "\n=== Fitting Distance Priors ===" << std::endl;
    // logger->info("Fitting gamma priors for {} distance(s)", distances.size());
    
    // Fit each distance separately
    for (int dist : distances) {
        if (dist <= 0) {
            Rcpp::Rcerr << "Warning: Skipping invalid distance " << dist << " (must be positive)" << std::endl;
            continue;
        }
        
        if (static_cast<arma::uword>(dist) >= positions.n_rows) {
            Rcpp::Rcerr << "Warning: Distance " << dist << " exceeds position matrix size, skipping" << std::endl;
            continue;
        }
        
        // Collect all distances at this separation
        arma::uword n_pairs = positions.n_rows - dist;
        arma::vec dist_values(n_pairs);
        
        for (arma::uword i = 0; i < n_pairs; ++i) {
            dist_values(i) = arma::norm(positions.row(i + dist) - positions.row(i), 2);
        }
        
        // Check for invalid distances
        if (arma::any(dist_values <= 0)) {
            Rcpp::Rcerr << "Warning: Found non-positive distances for separation " << dist << ", skipping" << std::endl;
            // logger->warn("Skipping distance {}: found non-positive values", dist);
            continue;
        }
        
        // Method of Moments (MOM) for Gamma distribution
        double mean = arma::mean(dist_values);
        double variance = arma::var(dist_values, 1);  // normalized by n-1
        
        if (variance <= 0) {
            Rcpp::Rcerr << "Warning: Non-positive variance for distance " << dist << ", skipping" << std::endl;
            // logger->warn("Skipping distance {}: variance = {}", dist, variance);
            continue;
        }
        
        double shape_mom = (mean * mean) / variance;
        double rate_mom = mean / variance;
        
        // Store the fitted prior
        distance_priors.set_prior(dist, GammaPrior(shape_mom, rate_mom));
        
        Rcpp::Rcout << "Distance " << dist << " (n=" << n_pairs << " pairs):" << std::endl;
        Rcpp::Rcout << "  Shape = " << shape_mom << ", Rate = " << rate_mom << std::endl;
        Rcpp::Rcout << "  Mean = " << mean << ", SD = " << std::sqrt(variance) << std::endl;
        
        // logger->info("Distance {}: shape={:.4f}, rate={:.4f}, mean={:.4f}, sd={:.4f} ({} pairs)", dist, shape_mom, rate_mom, mean, std::sqrt(variance), n_pairs);
    }
    
    Rcpp::Rcout << "=== Distance Prior Fitting Complete ===" << std::endl;
    
    return true;
}


// Legacy methods for backward compatibility
bool Chromosome::fit_adjacent_distance_prior_from_file(const std::string& positions_file) {
    return fit_distance_priors_from_file(positions_file, {1});
}


bool Chromosome::fit_adjacent_distance_prior(const arma::mat& positions) {
    return fit_distance_priors(positions, {1});
}


// ========================================================================
// MCMC Sampling Helper Methods
// ========================================================================

void Chromosome::sample_beta_parameters(double& current_ll, double sd_b0, double sd_b1, 
                                       int& accepted_b0, int& accepted_b1) {
    double delta_ll = 0.0;
    
    // --- Sample beta0 ---
    double proposed_beta0 = this->beta0 + arma::randn() * sd_b0;
    double old_beta0 = this->beta0;
    this->beta0 = proposed_beta0;
    delta_ll = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix) - current_ll;
    if (arma::randu() < std::exp(delta_ll)) {
        accepted_b0++;
        current_ll += delta_ll;
        
        // Update theta matrix if it exists (for convoluted sampling)
        if (!theta_matrix.is_empty()) {
            theta_matrix = arma::exp(beta0 + beta1 * arma::log(pairwise_distance_matrix + 1e-10));
            theta_matrix.diag().zeros();
        }
    } else {
        this->beta0 = old_beta0;
    }
    
    // --- Sample beta1 ---
    double proposed_beta1;
    do {
        proposed_beta1 = this->beta1 + arma::randn() * sd_b1;
    } while (proposed_beta1 >= 0);
    
    double old_beta1 = this->beta1;
    this->beta1 = proposed_beta1;
    delta_ll = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix) - current_ll;
    if (arma::randu() < std::exp(delta_ll)) {
        accepted_b1++;
        current_ll += delta_ll;
        
        // Update theta matrix if it exists (for convoluted sampling)
        if (!theta_matrix.is_empty()) {
            theta_matrix = arma::exp(beta0 + beta1 * arma::log(pairwise_distance_matrix + 1e-10));
            theta_matrix.diag().zeros();
        }
    } else {
        this->beta1 = old_beta1;
    }
}

void Chromosome::sample_cluster_centers(double& current_ll, double sd_center, int& accepted_center) {
    arma::uword num_clusters = cluster_center_position_matrix.n_rows;
    
    if (num_clusters <= 1) {
        return;  // No cluster centers to sample
    }
    
    for (arma::uword c = 0; c < num_clusters; ++c) {
        const auto& adj = cluster_adjacencies[c];
        const arma::uvec& self_idx = adj.self_indices;
        const arma::uvec& neighbor_idx = adj.neighbor_indices;
        
        arma::rowvec delta = arma::randn<arma::rowvec>(3) * sd_center;
        arma::mat proposed_pos_c = position_matrix.rows(self_idx);
        proposed_pos_c.each_row() += delta;
        arma::rowvec proposed_center_c = cluster_center_position_matrix.row(c) + delta;

        double delta_ll = 0.0;
        
        // Calculate likelihood change for neighbors
        arma::mat proposed_dists_neighbors;
        if (!neighbor_idx.is_empty()) {
            arma::mat current_dists_neighbors = pairwise_distance_matrix.submat(self_idx, neighbor_idx);
            proposed_dists_neighbors = calculate_pairwise_distances(proposed_pos_c, position_matrix.rows(neighbor_idx));
            delta_ll += 2 * calculate_log_likelihood(proposed_dists_neighbors, contact_matrix.submat(self_idx, neighbor_idx))
                      - 2 * calculate_log_likelihood(current_dists_neighbors, contact_matrix.submat(self_idx, neighbor_idx));
        }

        // Calculate likelihood change for stranger clusters
        arma::uvec unique_stranger_labels;
        if (!adj.stranger_indices.is_empty()) {
            unique_stranger_labels = arma::unique(cluster_labels(adj.stranger_indices));
        }
        
        for (arma::uword s_label : unique_stranger_labels) {
            arma::uvec s_idx = cluster_adjacencies[s_label].self_indices;
            double current_dist_cs = arma::norm(cluster_center_position_matrix.row(c) - cluster_center_position_matrix.row(s_label));
            double prop_dist_cs = arma::norm(proposed_center_c - cluster_center_position_matrix.row(s_label));
            arma::mat current_dist_mat(self_idx.n_elem, s_idx.n_elem, arma::fill::value(current_dist_cs));
            arma::mat prop_dist_mat(self_idx.n_elem, s_idx.n_elem, arma::fill::value(prop_dist_cs));
            delta_ll += 2 * calculate_log_likelihood(prop_dist_mat, contact_matrix.submat(self_idx, s_idx))
                      - 2 * calculate_log_likelihood(current_dist_mat, contact_matrix.submat(self_idx, s_idx));
        }

        // Accept or reject
        if (arma::randu() < std::exp(delta_ll)) {
            position_matrix.rows(self_idx) = proposed_pos_c;
            cluster_center_position_matrix.row(c) = proposed_center_c;
            
            if (!neighbor_idx.is_empty()) {
                pairwise_distance_matrix.submat(self_idx, neighbor_idx) = proposed_dists_neighbors;
                pairwise_distance_matrix.submat(neighbor_idx, self_idx) = proposed_dists_neighbors.t();
            }
            
            for (arma::uword s_label : unique_stranger_labels) {
                arma::uvec s_idx = cluster_adjacencies[s_label].self_indices;
                double new_dist_cs = arma::norm(proposed_center_c - cluster_center_position_matrix.row(s_label));
                pairwise_distance_matrix.submat(self_idx, s_idx).fill(new_dist_cs);
                pairwise_distance_matrix.submat(s_idx, self_idx).fill(new_dist_cs);
            }
            
            accepted_center++;
            current_ll += delta_ll;
        }
    }
}

void Chromosome::sample_locus_positions(double& current_ll, double sd_locus, int& accepted_locus) {
    arma::uword n_loci = position_matrix.n_rows;
    
    for (arma::uword j = 0; j < n_loci; ++j) {
        arma::uword c_j = cluster_labels(j);
        const auto& adj = cluster_adjacencies[c_j];
        
        // Get indices of neighbors (excluding j itself)
        arma::uvec self_idx_minus_j = adj.self_indices;
        self_idx_minus_j.shed_row(arma::as_scalar(arma::find(self_idx_minus_j == j)));
        arma::uvec combined_idx = arma::join_cols(self_idx_minus_j, adj.neighbor_indices);
        
        // Skip if this locus has zero contacts with all other loci (if option enabled)
        if (skip_zero_contact_loci) {
            arma::rowvec contacts_j_row = contact_matrix.row(j);
            if (arma::accu(contacts_j_row) == 0) {
                continue;  // No contact data for this locus
            }
        }
        
        arma::rowvec pos_j = position_matrix.row(j);
        arma::rowvec proposed_pos_j;
        
        // Sample new position: either from prior positions or random perturbation
        if (sample_from_prior && !prior_position_matrix.is_empty() && j < prior_position_matrix.n_rows) {
            // Sample from prior position for this locus with added noise
            proposed_pos_j = prior_position_matrix.row(j) + arma::randn<arma::rowvec>(3) * sd_locus;
        } else {
            // Standard random perturbation around current position
            proposed_pos_j = pos_j + arma::randn<arma::rowvec>(3) * sd_locus;
        }

        double delta_ll = 0.0;
        arma::mat proposed_dists;
        
        if (!combined_idx.is_empty()) {
            arma::mat current_dists = pairwise_distance_matrix.submat(arma::uvec{j}, combined_idx);
            arma::mat current_contacts = contact_matrix.submat(arma::uvec{j}, combined_idx);
            proposed_dists = calculate_pairwise_distances(proposed_pos_j, position_matrix.rows(combined_idx));
            delta_ll = 2 * (calculate_log_likelihood(proposed_dists, current_contacts) 
                          - calculate_log_likelihood(current_dists, current_contacts));
        }
        
        // Add distance prior contribution if informative priors exist
        if (distance_priors.has_informative_priors()) {
            double current_log_prior = distance_priors.log_prior_for_locus(position_matrix, pos_j, j);
            double proposed_log_prior = distance_priors.log_prior_for_locus(position_matrix, proposed_pos_j, j);
            delta_ll += (proposed_log_prior - current_log_prior);
        }

        // Accept or reject
        if (arma::randu() < std::exp(delta_ll)) {
            position_matrix.row(j) = proposed_pos_j;
            
            if (!combined_idx.is_empty()) {
                pairwise_distance_matrix.submat(arma::uvec{j}, combined_idx) = proposed_dists;
                pairwise_distance_matrix.submat(combined_idx, arma::uvec{j}) = proposed_dists.t();
            }
            
            accepted_locus++;
            current_ll += delta_ll;
        }
    }
}


// ========================================================================
// Sample Locus Positions Using Convoluted Matrices
// ========================================================================

void Chromosome::sample_locus_positions_convoluted(double& current_ll, double sd_locus, int& accepted_locus, int k) {
    arma::uword n_loci = position_matrix.n_rows;
    int half_k = k / 2;
    
    // Initialize theta matrix if empty
    if (theta_matrix.is_empty()) {
        theta_matrix = arma::exp(beta0 + beta1 * arma::log(pairwise_distance_matrix + 1e-10));
        theta_matrix.diag().zeros();  // Zero diagonal
        // logger->info("Initialized theta matrix for convoluted sampling");
    }
    
    for (arma::uword j = 0; j < n_loci; ++j) {
        arma::uword c_j = cluster_labels(j);
        const auto& adj = cluster_adjacencies[c_j];
        
        // Get indices of neighbors (excluding j itself)
        arma::uvec self_idx_minus_j = adj.self_indices;
        self_idx_minus_j.shed_row(arma::as_scalar(arma::find(self_idx_minus_j == j)));
        arma::uvec combined_idx = arma::join_cols(self_idx_minus_j, adj.neighbor_indices);
        
        // Skip if this locus has zero contacts with all other loci (if option enabled)
        if (skip_zero_contact_loci) {
            arma::rowvec contacts_j_row = contact_matrix.row(j);
            if (arma::accu(contacts_j_row) == 0) {
                continue;  // No contact data for this locus
            }
        }
        
        arma::rowvec pos_j = position_matrix.row(j);
        arma::rowvec proposed_pos_j;
        
        // Sample new position: either from prior positions or random perturbation
        if (sample_from_prior && !prior_position_matrix.is_empty() && j < prior_position_matrix.n_rows) {
            proposed_pos_j = prior_position_matrix.row(j) + arma::randn<arma::rowvec>(3) * sd_locus;
        } else {
            proposed_pos_j = pos_j + arma::randn<arma::rowvec>(3) * sd_locus;
        }

        // Compute new distances and thetas only for neighbor loci
        arma::mat proposed_dists;
        arma::mat proposed_thetas;
        if (!combined_idx.is_empty()) {
            proposed_dists = calculate_pairwise_distances(proposed_pos_j, position_matrix.rows(combined_idx));
            proposed_thetas = arma::exp(beta0 + beta1 * arma::log(proposed_dists + 1e-10));
        }
        
        double delta_ll = 0.0;
        
        // For each neighbor pair (i, m), check if their convoluted theta is affected
        // Only consider neighbors to reduce computation
        for (arma::uword idx1 = 0; idx1 < combined_idx.n_elem; ++idx1) {
            arma::uword i = combined_idx(idx1);
            
            for (arma::uword idx2 = idx1 + 1; idx2 < combined_idx.n_elem; ++idx2) {
                arma::uword m = combined_idx(idx2);
                
                // Define window boundaries for convoluted cell (i, m)
                arma::uword i_start = (i >= static_cast<arma::uword>(half_k)) ? i - half_k : 0;
                arma::uword i_end = std::min(n_loci, i + half_k + 1);
                arma::uword m_start = (m >= static_cast<arma::uword>(half_k)) ? m - half_k : 0;
                arma::uword m_end = std::min(n_loci, m + half_k + 1);
                
                // Check if window intersects with locus j (row j or column j in the crucifix)
                bool window_contains_j = (j >= i_start && j < i_end) || (j >= m_start && j < m_end);
                
                if (!window_contains_j) {
                    continue;  // This convoluted cell is unaffected
                }
                
                // Compute current convoluted theta (sum over window, upper triangular only)
                double current_conv_theta = 0.0;
                for (arma::uword r = i_start; r < i_end; ++r) {
                    for (arma::uword c = m_start; c < m_end; ++c) {
                        if (r < c) {  // upper triangular only
                            current_conv_theta += theta_matrix(r, c);
                        }
                    }
                }
                
                // Compute proposed convoluted theta (substitute affected thetas involving j)
                double proposed_conv_theta = current_conv_theta;
                for (arma::uword r = i_start; r < i_end; ++r) {
                    for (arma::uword c = m_start; c < m_end; ++c) {
                        if (r < c) {  // upper triangular only
                            // Check if this element involves locus j
                            if (r == j || c == j) {
                                // Subtract old theta
                                proposed_conv_theta -= theta_matrix(r, c);
                                
                                // Add new theta (lookup from proposed_thetas)
                                arma::uword other_locus = (r == j) ? c : r;
                                auto it = std::find(combined_idx.begin(), combined_idx.end(), other_locus);
                                if (it != combined_idx.end()) {
                                    arma::uword other_idx = std::distance(combined_idx.begin(), it);
                                    proposed_conv_theta += proposed_thetas(0, other_idx);
                                }
                            }
                        }
                    }
                }
                
                // Compute likelihood change for this convoluted cell
                double c_im = convoluted_contact_matrix(i, m);
                if (current_conv_theta > 1e-10) {
                    delta_ll -= c_im * std::log(current_conv_theta) - current_conv_theta;
                }
                if (proposed_conv_theta > 1e-10) {
                    delta_ll += c_im * std::log(proposed_conv_theta) - proposed_conv_theta;
                }
            }
        }
        
        // Add distance prior contribution if informative priors exist
        if (distance_priors.has_informative_priors()) {
            double current_log_prior = distance_priors.log_prior_for_locus(position_matrix, pos_j, j);
            double proposed_log_prior = distance_priors.log_prior_for_locus(position_matrix, proposed_pos_j, j);
            delta_ll += (proposed_log_prior - current_log_prior);
        }

        // Accept or reject
        if (arma::randu() < std::exp(delta_ll)) {
            // Update position matrix
            position_matrix.row(j) = proposed_pos_j;
            
            // Update pairwise distance matrix (only for neighbors)
            if (!combined_idx.is_empty()) {
                pairwise_distance_matrix.submat(arma::uvec{j}, combined_idx) = proposed_dists;
                pairwise_distance_matrix.submat(combined_idx, arma::uvec{j}) = proposed_dists.t();
            }
            pairwise_distance_matrix(j, j) = 0;
            
            // Update theta matrix (only for neighbors)
            if (!combined_idx.is_empty()) {
                theta_matrix.submat(arma::uvec{j}, combined_idx) = proposed_thetas;
                theta_matrix.submat(combined_idx, arma::uvec{j}) = proposed_thetas.t();
            }
            theta_matrix(j, j) = 0;
            
            accepted_locus++;
            current_ll += delta_ll;
        }
    }
}


// ========================================================================
// Main MCMC Function
// ========================================================================


void Chromosome::run_mcmc(int iterations, int burn_in, double initial_sd, double sd_floor, double sd_ceiling, bool save_samples, int sample_interval) {
    if (pairwise_distance_matrix.is_empty()) {
        Rcpp::Rcerr << "Error: Positions must be initialized before running MCMC." << std::endl;
        return;
    }
    double current_ll = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    Rcpp::Rcout << "\n--- Starting MCMC Sampling ---" << std::endl;
    // logger->info("--- Starting MCMC Sampling ---");
    // logger->info("Initial beta0: {:.4f}, Initial beta1: {:.4f}, Initial loglikelihood: {:.4f}", this->beta0, this->beta1, current_ll);
    // logger->info("Iterations: {}, Burn-in: {}", iterations, burn_in);
    
    // Set up sample saving if requested
    std::string samples_file;
    std::ofstream samples_stream;
    int samples_saved = 0;
    
    if (save_samples) {
        // Create samples directory for MCMC trace
        std::string samples_dir = chromosome_name + "/samples";
        std::filesystem::create_directories(samples_dir);
        
        // Create timestamped filename
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        
        std::stringstream ss;
        ss << samples_dir << "/mcmc_samples_" << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S") 
           << "_" << std::setfill('0') << std::setw(3) << ms.count() << ".txt";
        samples_file = ss.str();
        
        samples_stream.open(samples_file);
        if (!samples_stream.is_open()) {
            Rcpp::Rcerr << "Warning: Could not open samples file " << samples_file << ". Continuing without saving samples." << std::endl;
            save_samples = false;
        } else {
            // logger->info("Saving MCMC samples to: {}", samples_file);
            // logger->info("Sample interval: every {} iterations after burn-in", sample_interval);
            
            // Write header
            samples_stream << "# MCMC Samples - HiSpa Chromosome Analysis\n";
            samples_stream << "# Generated: " << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S") << "\n";
            samples_stream << "# Iterations: " << iterations << ", Burn-in: " << burn_in << ", Sample interval: " << sample_interval << "\n";
            samples_stream << "# Matrix dimensions: " << position_matrix.n_rows << " loci x " << position_matrix.n_cols << " coordinates (X,Y,Z)\n";
            samples_stream << "# Columns: iteration log_likelihood beta0 beta1 [position_matrix: locus1_x locus1_y locus1_z locus2_x locus2_y locus2_z ...]\n";
            samples_stream << std::fixed << std::setprecision(6);
        }
    }

    // Create intermediate results directory for best structures
    std::string intermediate_dir = chromosome_name + "/intermediate_results";
    std::filesystem::create_directories(intermediate_dir);
    // logger->info("Saving intermediate best results to: {}", intermediate_dir);

    // --- Adaptive MCMC setup ---
    double sd_b0 = initial_sd, sd_b1 = initial_sd, sd_center = initial_sd, sd_locus = initial_sd;
    int accepted_b0 = 0, accepted_b1 = 0, accepted_center = 0, accepted_locus = 0;
    int total_accepted_b0 = 0, total_accepted_b1 = 0, total_accepted_center = 0, total_accepted_locus = 0;
    
    mcmc_trace_beta0.set_size(iterations - burn_in);
    mcmc_trace_beta1.set_size(iterations - burn_in);
    mcmc_trace_log_likelihood.set_size(iterations - burn_in);
    mcmc_trace_cluster_centers.reserve(iterations - burn_in);
    mcmc_trace_block_durations.reserve((iterations / 50) + 1);

    arma::uword num_clusters = cluster_center_position_matrix.n_rows;
    arma::uword n_loci = position_matrix.n_rows;

    auto block_start_time = std::chrono::steady_clock::now();
    auto total_start_time = std::chrono::steady_clock::now();
    
    for (int i = 0; i < iterations; ++i) {
        // Sample beta parameters
        int iter_accepted_b0 = 0, iter_accepted_b1 = 0;
        sample_beta_parameters(current_ll, sd_b0, sd_b1, iter_accepted_b0, iter_accepted_b1);
        accepted_b0 += iter_accepted_b0;
        accepted_b1 += iter_accepted_b1;
        total_accepted_b0 += iter_accepted_b0;
        total_accepted_b1 += iter_accepted_b1;
        
        // Sample cluster center positions
        int iter_accepted_center = 0;
        sample_cluster_centers(current_ll, sd_center, iter_accepted_center);
        accepted_center += iter_accepted_center;
        total_accepted_center += iter_accepted_center;
        
        // Sample individual locus positions
        int iter_accepted_locus = 0;
        sample_locus_positions(current_ll, sd_locus, iter_accepted_locus);
        accepted_locus += iter_accepted_locus;
        total_accepted_locus += iter_accepted_locus;

        // --- Track best state and store samples ---
        if (current_ll > max_log_likelihood) {
            max_log_likelihood = current_ll;
            best_position_matrix = position_matrix;
        }

        // Save samples after burn-in period
        if (save_samples && i >= burn_in && (i - burn_in) % sample_interval == 0) {
            samples_stream << (i + 1) << " " << current_ll << " " << this->beta0 << " " << this->beta1;
            
            // Save position matrix (flattened row-wise: x1 y1 z1 x2 y2 z2 ...)
            for (arma::uword row = 0; row < position_matrix.n_rows; ++row) {
                for (arma::uword col = 0; col < position_matrix.n_cols; ++col) {
                    samples_stream << " " << position_matrix(row, col);
                }
            }
            
            samples_stream << "\n";
            samples_saved++;
            
            // Periodic flush to ensure data is written
            if (samples_saved % 100 == 0) {
                samples_stream.flush();
            }
        }

        // --- Adapt proposal SDs and log block progress ---
        if ((i + 1) % 50 == 0) {
            auto block_end_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> block_duration = block_end_time - block_start_time;
            mcmc_trace_block_durations.push_back(block_duration.count());
            
            double rate_b0 = (double)accepted_b0 / 50.0;
            double rate_b1 = (double)accepted_b1 / 50.0;
            double rate_center = (num_clusters > 1) ? (double)accepted_center / (50.0 * num_clusters) : 0.0;
            double rate_locus = (double)accepted_locus / (50.0 * n_loci);

            // logger->info("-------------------- Block Summary (Iter {}) --------------------", i + 1);
            // logger->info("Time for last 50 iterations: {:.2f} seconds.", block_duration.count());
            // logger->info("Acceptance Rates -> beta0: {:.2f}%, beta1: {:.2f}%, center: {:.2f}%, locus: {:.2f}%", rate_b0 * 100.0, rate_b1 * 100.0, rate_center * 100.0, rate_locus * 100.0);
            // logger->info("Proposal SDs -> beta0: {:.4f}, beta1: {:.4f}, center: {:.4f}, locus: {:.4f}", sd_b0, sd_b1, sd_center, sd_locus);
            // logger->info("Current State -> max LL: {:.4f}, beta0: {:.4f}, beta1: {:.4f}", max_log_likelihood, this->beta0, this->beta1);

            // Save intermediate best results
            std::string iter_str = std::to_string(i + 1);
            best_position_matrix.save(intermediate_dir + "/positions_iter" + iter_str + ".txt", arma::raw_ascii);
            
            // Save parameters
            std::ofstream param_file(intermediate_dir + "/parameters_iter" + iter_str + ".txt");
            param_file << "iteration " << (i + 1) << std::endl;
            param_file << "max_log_likelihood " << max_log_likelihood << std::endl;
            param_file << "current_log_likelihood " << current_ll << std::endl;
            param_file << "beta0 " << this->beta0 << std::endl;
            param_file << "beta1 " << this->beta1 << std::endl;
            param_file.close();

            if (rate_b0 > 0.3) sd_b0 = std::min(sd_b0*2.0, sd_ceiling); else if (rate_b0 < 0.2) sd_b0 = std::max(sd_b0/2.0, sd_floor);
            if (rate_b1 > 0.3) sd_b1 = std::min(sd_b1*2.0, sd_ceiling); else if (rate_b1 < 0.2) sd_b1 = std::max(sd_b1/2.0, sd_floor);
            if (num_clusters > 1) {
                if (rate_center > 0.3) sd_center = std::min(sd_center*2.0, sd_ceiling); else if (rate_center < 0.2) sd_center = std::max(sd_center/2.0, sd_floor);
            }
            if (rate_locus > 0.3) sd_locus = std::min(sd_locus*2.0, sd_ceiling); else if (rate_locus < 0.2) sd_locus = std::max(sd_locus/2.0, sd_floor);
            
            accepted_b0=0; accepted_b1=0; accepted_center=0; accepted_locus=0;
            block_start_time = std::chrono::steady_clock::now();
        }

        if ((i + 1) % 1000 == 0) {
            Rcpp::Rcout << "Iter " << i+1 << "/" << iterations << ". LL: " << max_log_likelihood << ". Best LL: " << max_log_likelihood << std::endl;
        }
    }

    // Close samples file if it was opened
    if (save_samples && samples_stream.is_open()) {
        samples_stream.close();
        Rcpp::Rcout << "Saved " << samples_saved << " samples to: " << samples_file << std::endl;
        // logger->info("Saved {} samples to: {}", samples_saved, samples_file);
    }

    Rcpp::Rcout << "--- MCMC Finished ---" << std::endl;
    Rcpp::Rcout << "Final max log-likelihood found: " << max_log_likelihood << std::endl;
    Rcpp::Rcout << "Beta0 acceptance rate: " << (double)total_accepted_b0 / iterations * 100.0 << "%" << std::endl;
    Rcpp::Rcout << "Beta1 acceptance rate: " << (double)total_accepted_b1 / iterations * 100.0 << "%" << std::endl;
    if (num_clusters > 1) {
        Rcpp::Rcout << "Center pos acceptance rate: " << (double)total_accepted_center / (iterations * num_clusters) * 100.0 << "%" << std::endl;
    }
    Rcpp::Rcout << "Locus pos acceptance rate: " << (double)total_accepted_locus / (iterations * n_loci) * 100.0 << "%" << std::endl;
    auto total_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_duration = total_end_time - total_start_time;
    // logger->info("Total MCMC duration: {:.2f} seconds.", total_duration.count());
    best_position_matrix.save(chromosome_name + "/final_positions.txt", arma::raw_ascii);
    // Save final parameters
    std::ofstream final_param_file(chromosome_name + "/final_parameters.txt");
    final_param_file << "max_log_likelihood " << max_log_likelihood << std::endl;
    final_param_file << "beta0 " << this->beta0 << std::endl;
    final_param_file << "beta1 " << this->beta1 << std::endl;
    final_param_file.close();
    log_memory_usage();
}

void Chromosome::log_memory_usage() {
    #if defined(__APPLE__) && defined(__MACH__)
        mach_task_basic_info_data_t info;
        mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) == KERN_SUCCESS) {
            double memory_usage_mb = static_cast<double>(info.resident_size) / (1024 * 1024);
            // logger->info("Peak memory usage (RSS): {:.2f} MB", memory_usage_mb);
        }
    #elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
        std::ifstream status_file("/proc/self/status");
        std::string line;
        while (std::getline(status_file, line)) {
            if (line.rfind("VmRSS:", 0) == 0) {
                std::istringstream iss(line);
                std::string key;
                long memory_kb;
                iss >> key >> memory_kb;
                double memory_mb = static_cast<double>(memory_kb) / 1024.0;
                // logger->info("Peak memory usage (RSS): {:.2f} MB", memory_mb);
                break;
            }
        }
    #else
        // logger->info("Memory usage reporting not supported on this platform.");
    #endif
}


// ========================================================================
// Run MCMC with Convoluted Sampling
// ========================================================================

void Chromosome::run_mcmc_convoluted(int iterations, int burn_in, double initial_sd, double sd_floor, double sd_ceiling, bool save_samples, int sample_interval, int k) {
    if (pairwise_distance_matrix.is_empty()) {
        Rcpp::Rcerr << "Error: Positions must be initialized before running MCMC." << std::endl;
        return;
    }
    
    // Compute convoluted contact matrix
    if (convoluted_contact_matrix.is_empty()) {
        Rcpp::Rcout << "Computing convoluted contact matrix (half_k=" << k << ")..." << std::endl;
        // logger->info("Computing convoluted contact matrix with half_k={}", k);
        compute_convoluted_contact_matrix(k);
    }
    
    // Initialize theta matrix for convoluted sampling
    theta_matrix = arma::exp(beta0 + beta1 * arma::log(pairwise_distance_matrix + 1e-10));
    theta_matrix.diag().zeros();
    // logger->info("Initialized theta matrix for convoluted sampling");
    
    // Compute initial log-likelihood using convoluted matrices
    // For convoluted likelihood: sum over (i,j) of [c_conv(i,j) * log(theta_conv(i,j)) - theta_conv(i,j)]
    // where theta_conv(i,j) = sum of thetas in window around (i,j)
    arma::mat convoluted_theta = convolute_contacts(theta_matrix, k);
    double current_ll = 0.0;
    arma::uword n = convoluted_contact_matrix.n_rows;
    for (arma::uword i = 0; i < n; ++i) {
        for (arma::uword j = i + 1; j < n; ++j) {
            double c_ij = convoluted_contact_matrix(i, j);
            double theta_ij = convoluted_theta(i, j);
            if (theta_ij > 1e-10) {
                current_ll += c_ij * std::log(theta_ij) - theta_ij;
            }
        }
    }
    
    Rcpp::Rcout << "\n--- Starting MCMC Sampling (Convoluted) ---" << std::endl;
    // logger->info("--- Starting MCMC Sampling (Convoluted) ---");
    // logger->info("Convolution half_window_size: {}", k);
    // logger->info("Initial beta0: {:.4f}, Initial beta1: {:.4f}, Initial loglikelihood: {:.4f}", this->beta0, this->beta1, current_ll);
    // logger->info("Iterations: {}, Burn-in: {}", iterations, burn_in);
    
    // Set up sample saving if requested
    std::string samples_file;
    std::ofstream samples_stream;
    int samples_saved = 0;
    
    if (save_samples) {
        std::string samples_dir = chromosome_name + "/samples";
        std::filesystem::create_directories(samples_dir);
        
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        
        std::stringstream ss;
        ss << samples_dir << "/mcmc_samples_convoluted_" << std::put_time(std::localtime(&time_t), "%Y%m%d_%H%M%S") 
           << "_" << std::setfill('0') << std::setw(3) << ms.count() << ".txt";
        samples_file = ss.str();
        
        samples_stream.open(samples_file);
        if (!samples_stream.is_open()) {
            Rcpp::Rcerr << "Warning: Could not open samples file " << samples_file << ". Continuing without saving samples." << std::endl;
            save_samples = false;
        } else {
            // logger->info("Saving MCMC samples to: {}", samples_file);
            // logger->info("Sample interval: every {} iterations after burn-in", sample_interval);
            
            samples_stream << "# MCMC Samples (Convoluted) - HiSpa Chromosome Analysis\n";
            samples_stream << "# Generated: " << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S") << "\n";
            samples_stream << "# Convolution half_k: " << k << "\n";
            samples_stream << "# Iterations: " << iterations << ", Burn-in: " << burn_in << ", Sample interval: " << sample_interval << "\n";
            samples_stream << "# Matrix dimensions: " << position_matrix.n_rows << " loci x " << position_matrix.n_cols << " coordinates (X,Y,Z)\n";
            samples_stream << "# Columns: iteration log_likelihood beta0 beta1 [position_matrix: locus1_x locus1_y locus1_z locus2_x locus2_y locus2_z ...]\n";
            samples_stream << std::fixed << std::setprecision(6);
        }
    }

    std::string intermediate_dir = chromosome_name + "/intermediate_results";
    std::filesystem::create_directories(intermediate_dir);
    // logger->info("Saving intermediate best results to: {}", intermediate_dir);

    // --- Adaptive MCMC setup ---
    double sd_b0 = initial_sd, sd_b1 = initial_sd, sd_center = initial_sd, sd_locus = initial_sd;
    int accepted_b0 = 0, accepted_b1 = 0, accepted_center = 0, accepted_locus = 0;
    int total_accepted_b0 = 0, total_accepted_b1 = 0, total_accepted_center = 0, total_accepted_locus = 0;
    
    mcmc_trace_beta0.set_size(iterations - burn_in);
    mcmc_trace_beta1.set_size(iterations - burn_in);
    mcmc_trace_log_likelihood.set_size(iterations - burn_in);
    mcmc_trace_cluster_centers.reserve(iterations - burn_in);
    mcmc_trace_block_durations.reserve((iterations / 50) + 1);

    arma::uword num_clusters = cluster_center_position_matrix.n_rows;
    arma::uword n_loci = position_matrix.n_rows;

    auto block_start_time = std::chrono::steady_clock::now();
    auto total_start_time = std::chrono::steady_clock::now();
    
    for (int i = 0; i < iterations; ++i) {
        // Sample beta parameters (updates theta_matrix internally when accepted)
        int iter_accepted_b0 = 0, iter_accepted_b1 = 0;
        sample_beta_parameters(current_ll, sd_b0, sd_b1, iter_accepted_b0, iter_accepted_b1);
        accepted_b0 += iter_accepted_b0;
        accepted_b1 += iter_accepted_b1;
        total_accepted_b0 += iter_accepted_b0;
        total_accepted_b1 += iter_accepted_b1;
        
        // Sample cluster center positions (not yet implemented for convoluted version)
        // TODO: Implement sample_cluster_centers_convoluted
        // int iter_accepted_center = 0;
        // For now, skip cluster center sampling in convoluted mode
        // sample_cluster_centers(current_ll, sd_center, iter_accepted_center);
        // accepted_center += iter_accepted_center;
        // total_accepted_center += iter_accepted_center;
        
        // Sample individual locus positions using convoluted matrices
        int iter_accepted_locus = 0;
        sample_locus_positions_convoluted(current_ll, sd_locus, iter_accepted_locus, k);
        accepted_locus += iter_accepted_locus;
        total_accepted_locus += iter_accepted_locus;

        // Track best state
        if (current_ll > max_log_likelihood) {
            max_log_likelihood = current_ll;
            best_position_matrix = position_matrix;
        }

        // Save samples after burn-in
        if (save_samples && i >= burn_in && (i - burn_in) % sample_interval == 0) {
            samples_stream << (i + 1) << " " << current_ll << " " << this->beta0 << " " << this->beta1;
            
            for (arma::uword row = 0; row < position_matrix.n_rows; ++row) {
                for (arma::uword col = 0; col < position_matrix.n_cols; ++col) {
                    samples_stream << " " << position_matrix(row, col);
                }
            }
            
            samples_stream << "\n";
            samples_saved++;
            
            if (samples_saved % 100 == 0) {
                samples_stream.flush();
            }
        }

        // Adapt proposal SDs and log progress
        if ((i + 1) % 50 == 0) {
            auto block_end_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> block_duration = block_end_time - block_start_time;
            mcmc_trace_block_durations.push_back(block_duration.count());
            
            double rate_b0 = (double)accepted_b0 / 50.0;
            double rate_b1 = (double)accepted_b1 / 50.0;
            double rate_center = (num_clusters > 1) ? (double)accepted_center / (50.0 * num_clusters) : 0.0;
            double rate_locus = (double)accepted_locus / (50.0 * n_loci);

            // logger->info("-------------------- Block Summary (Iter {}) --------------------", i + 1);
            // logger->info("Time for last 50 iterations: {:.2f} seconds.", block_duration.count());
            // logger->info("Acceptance Rates -> beta0: {:.2f}%, beta1: {:.2f}%, center: {:.2f}%, locus: {:.2f}%", rate_b0 * 100.0, rate_b1 * 100.0, rate_center * 100.0, rate_locus * 100.0);
            // logger->info("Proposal SDs -> beta0: {:.4f}, beta1: {:.4f}, center: {:.4f}, locus: {:.4f}", sd_b0, sd_b1, sd_center, sd_locus);
            // logger->info("Current State -> max LL: {:.4f}, beta0: {:.4f}, beta1: {:.4f}", max_log_likelihood, this->beta0, this->beta1);

            // Save intermediate results
            std::string iter_str = std::to_string(i + 1);
            best_position_matrix.save(intermediate_dir + "/positions_iter" + iter_str + ".txt", arma::raw_ascii);
            
            std::ofstream param_file(intermediate_dir + "/parameters_iter" + iter_str + ".txt");
            param_file << "iteration " << (i + 1) << std::endl;
            param_file << "max_log_likelihood " << max_log_likelihood << std::endl;
            param_file << "current_log_likelihood " << current_ll << std::endl;
            param_file << "beta0 " << this->beta0 << std::endl;
            param_file << "beta1 " << this->beta1 << std::endl;
            param_file << "convolution_half_k " << k << std::endl;
            param_file.close();

            if (rate_b0 > 0.3) sd_b0 = std::min(sd_b0*2.0, sd_ceiling); else if (rate_b0 < 0.2) sd_b0 = std::max(sd_b0/2.0, sd_floor);
            if (rate_b1 > 0.3) sd_b1 = std::min(sd_b1*2.0, sd_ceiling); else if (rate_b1 < 0.2) sd_b1 = std::max(sd_b1/2.0, sd_floor);
            if (num_clusters > 1) {
                if (rate_center > 0.3) sd_center = std::min(sd_center*2.0, sd_ceiling); else if (rate_center < 0.2) sd_center = std::max(sd_center/2.0, sd_floor);
            }
            if (rate_locus > 0.3) sd_locus = std::min(sd_locus*2.0, sd_ceiling); else if (rate_locus < 0.2) sd_locus = std::max(sd_locus/2.0, sd_floor);
            
            accepted_b0=0; accepted_b1=0; accepted_center=0; accepted_locus=0;
            block_start_time = std::chrono::steady_clock::now();
        }

        if ((i + 1) % 1000 == 0) {
            Rcpp::Rcout << "Iter " << i+1 << "/" << iterations << ". LL: " << current_ll << ". Best LL: " << max_log_likelihood << std::endl;
        }
    }

    if (save_samples && samples_stream.is_open()) {
        samples_stream.close();
        Rcpp::Rcout << "Saved " << samples_saved << " samples to: " << samples_file << std::endl;
        // logger->info("Saved {} samples to: {}", samples_saved, samples_file);
    }

    Rcpp::Rcout << "--- MCMC Finished (Convoluted) ---" << std::endl;
    Rcpp::Rcout << "Final max log-likelihood found: " << max_log_likelihood << std::endl;
    Rcpp::Rcout << "Beta0 acceptance rate: " << (double)total_accepted_b0 / iterations * 100.0 << "%" << std::endl;
    Rcpp::Rcout << "Beta1 acceptance rate: " << (double)total_accepted_b1 / iterations * 100.0 << "%" << std::endl;
    if (num_clusters > 1) {
        Rcpp::Rcout << "Center pos acceptance rate: " << (double)total_accepted_center / (iterations * num_clusters) * 100.0 << "%" << std::endl;
    }
    Rcpp::Rcout << "Locus pos acceptance rate: " << (double)total_accepted_locus / (iterations * n_loci) * 100.0 << "%" << std::endl;
    auto total_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_duration = total_end_time - total_start_time;
    // logger->info("Total MCMC duration: {:.2f} seconds.", total_duration.count());
    log_memory_usage();
}

// --- GETTERS ---
const std::string& Chromosome::get_name() const { return chromosome_name; }
const arma::mat& Chromosome::get_contact_matrix() const { return contact_matrix; }
const arma::uvec& Chromosome::get_cluster_labels() const { return cluster_labels; }
const std::vector<ClusterAdjacency>& Chromosome::get_cluster_adjacencies() const { return cluster_adjacencies; }
const arma::mat& Chromosome::get_backbone_matrix() const { return backbone_contact_matrix; }
const arma::mat& Chromosome::get_position_matrix() const { return position_matrix; }
const arma::mat& Chromosome::get_cluster_center_position_matrix() const { return cluster_center_position_matrix; }
const arma::mat& Chromosome::get_pairwise_distance_matrix() const { return pairwise_distance_matrix; }
double Chromosome::get_beta0() const { return beta0; }
double Chromosome::get_beta1() const { return beta1; }
const arma::vec& Chromosome::get_mcmc_trace_beta0() const { return mcmc_trace_beta0; }
const arma::vec& Chromosome::get_mcmc_trace_beta1() const { return mcmc_trace_beta1; }
const std::vector<arma::mat>& Chromosome::get_mcmc_trace_cluster_centers() const { return mcmc_trace_cluster_centers; }
const arma::mat& Chromosome::get_best_position_matrix() const { return best_position_matrix; }
const arma::vec& Chromosome::get_mcmc_trace_log_likelihood() const { Rcpp::Rcout << "mcmc trace length: " << mcmc_trace_log_likelihood.n_elem << std::endl; return mcmc_trace_log_likelihood; }
const std::vector<double>& Chromosome::get_mcmc_trace_block_durations() const { return mcmc_trace_block_durations; }
const std::vector<arma::mat>& Chromosome::get_initial_cluster_structures() const { return initial_cluster_structures; }
const arma::mat& Chromosome::get_backbone_structure() const { return backbone_structure; }
const std::vector<double>& Chromosome::get_initial_beta0s() const { return initial_beta0s; }
const std::vector<double>& Chromosome::get_initial_beta1s() const { return initial_beta1s; }
double Chromosome::get_backbone_beta0() const { return backbone_beta0; }
double Chromosome::get_backbone_beta1() const { return backbone_beta1; }
const DistancePriors& Chromosome::get_distance_priors() const { return distance_priors; }
const GammaPrior& Chromosome::get_adjacent_distance_prior() const { 
    // Legacy accessor - returns adjacent distance prior (distance=1)
    return distance_priors.get_prior(1); 
}
const arma::mat& Chromosome::get_prior_position_matrix() const { return prior_position_matrix; }
const arma::mat& Chromosome::get_prior_contact_matrix() const { return prior_contact_matrix; }
bool Chromosome::get_skip_zero_contact_loci() const { return skip_zero_contact_loci; }


// --- SETTERS ---
void Chromosome::set_beta0(double val) { this->beta0 = val; }
void Chromosome::set_beta1(double val) { this->beta1 = val; }

void Chromosome::set_distance_priors(const DistancePriors& priors) {
    distance_priors = priors;
    // logger->info("Set distance priors with {} active distance(s)", priors.active_distances.size());
}

void Chromosome::set_distance_prior(int distance, double shape, double rate) {
    distance_priors.set_prior(distance, GammaPrior(shape, rate));
    // logger->info("Set distance prior for d={}: shape={:.4f}, rate={:.4f}", distance, shape, rate);
}

void Chromosome::set_distance_prior(int distance, const GammaPrior& prior) {
    distance_priors.set_prior(distance, prior);
    // logger->info("Set distance prior for d={}: shape={:.4f}, rate={:.4f}", distance, prior.shape, prior.rate);
}

// Legacy setters for backward compatibility
void Chromosome::set_adjacent_distance_prior(double shape, double rate) {
    set_distance_prior(1, shape, rate);
}

void Chromosome::set_adjacent_distance_prior(const GammaPrior& prior) {
    set_distance_prior(1, prior);
}

void Chromosome::set_position_matrix(const arma::mat& positions) {
    if (positions.n_rows != contact_matrix.n_rows) {
        Rcpp::Rcerr << "Error: Position matrix has " << positions.n_rows 
                  << " rows, but contact matrix has " << contact_matrix.n_rows << " rows" << std::endl;
        return;
    }
    if (positions.n_cols != 3) {
        Rcpp::Rcerr << "Error: Position matrix must have 3 columns (x, y, z), but has " 
                  << positions.n_cols << " columns" << std::endl;
        return;
    }
    
    position_matrix = positions;
    
    // Update pairwise distance matrix
    pairwise_distance_matrix = calculate_pairwise_distances(position_matrix, position_matrix);
    
    // logger->info("Set position matrix ({} x {}) and updated pairwise distances",  positions.n_rows, positions.n_cols);
}

void Chromosome::set_prior_contact_matrix(const arma::mat& contacts) {
    if (contacts.n_rows != contacts.n_cols) {
        Rcpp::Rcerr << "Error: Contact matrix must be square, but has dimensions " 
                  << contacts.n_rows << " x " << contacts.n_cols << std::endl;
        // logger->error("Failed to set prior contact matrix: not square ({} x {})", contacts.n_rows, contacts.n_cols);
        return;
    }
    
    prior_contact_matrix = contacts;
    
    // logger->info("Set prior contact matrix ({} x {})", contacts.n_rows, contacts.n_cols);
    Rcpp::Rcout << "Prior contact matrix stored (" << contacts.n_rows << " x " 
              << contacts.n_cols << ")" << std::endl;
}

bool Chromosome::load_prior_contact_matrix_from_file(const std::string& filename) {
    arma::mat contacts;
    bool loaded = contacts.load(filename, arma::raw_ascii);
    
    if (!loaded) {
        Rcpp::Rcerr << "Error: Could not load prior contact matrix from " << filename << std::endl;
        // logger->error("Failed to load prior contact matrix from {}", filename);
        return false;
    }
    
    if (contacts.n_rows != contacts.n_cols) {
        Rcpp::Rcerr << "Error: Prior contact matrix must be square, but has dimensions " 
                  << contacts.n_rows << " x " << contacts.n_cols << std::endl;
        // logger->error("Prior contact matrix from {} is not square", filename);
        return false;
    }
    
    set_prior_contact_matrix(contacts);
    Rcpp::Rcout << "Successfully loaded prior contact matrix from: " << filename << std::endl;
    
    return true;
}

void Chromosome::set_skip_zero_contact_loci(bool skip) {
    skip_zero_contact_loci = skip;
    // logger->info("Set skip_zero_contact_loci to: {}", skip);
    Rcpp::Rcout << "Skip zero-contact loci: " << (skip ? "enabled" : "disabled") << std::endl;
}

void Chromosome::set_sample_from_prior(bool sample) {
    sample_from_prior = sample;
    // logger->info("Set sample_from_prior to: {}", sample);
    Rcpp::Rcout << "Sample from prior positions: " << (sample ? "enabled" : "disabled") << std::endl;
}

bool Chromosome::get_sample_from_prior() const {
    return sample_from_prior;
}

// =========================================================================
// Getters and Setters for Convoluted Matrices
// =========================================================================

const arma::mat& Chromosome::get_convoluted_contact_matrix() const {
    return convoluted_contact_matrix;
}

void Chromosome::set_convoluted_contact_matrix(const arma::mat& convoluted_contacts) {
    convoluted_contact_matrix = convoluted_contacts;
    // logger->info("Set convoluted contact matrix ({} x {})", convoluted_contacts.n_rows, convoluted_contacts.n_cols);
}

void Chromosome::compute_convoluted_contact_matrix(int k) {
    if (contact_matrix.is_empty()) {
        Rcpp::Rcerr << "Error: Cannot compute convoluted contact matrix - contact matrix is empty" << std::endl;
        // logger->error("Failed to compute convoluted contact matrix: contact matrix is empty");
        return;
    }
    
    convoluted_contact_matrix = convolute_contacts(contact_matrix, k);
    // logger->info("Computed convoluted contact matrix with half_window_size={}", k);
    Rcpp::Rcout << "Convoluted contact matrix computed (half_window_size=" << k << ")" << std::endl;
}

