# HiSpaR - R Package for Hi-C Spatial Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HiSpaR is an R package that provides R bindings for the HiSpa C++ library, enabling hierarchical Bayesian inference of 3D chromatin structures from Hi-C contact matrices.

## Features

- 🧬 **Hierarchical Bayesian modeling** of 3D chromatin structure
- 🎯 **MCMC sampling** with adaptive proposals for efficient inference
- 📊 **Cluster-based analysis** for large-scale Hi-C data
- 🔄 **Prior information integration** from reference datasets
- ⚡ **High-performance C++** backend via Rcpp and RcppArmadillo
- 🔧 **OpenMP parallelization** for multi-core processing

## Installation

### Prerequisites

Before installing HiSpaR, ensure you have the following:

#### System Requirements
- **R** (≥ 4.0.0)
- **C++17 compatible compiler**
  - macOS: Clang 12.0+ (comes with Xcode Command Line Tools)
  - Linux: GCC 7.0+ or Clang 5.0+
  - Windows: Rtools 4.0+

#### Required Libraries

**Armadillo C++ Library** (≥ 9.0):
```bash
# macOS (Homebrew)
brew install armadillo

# Ubuntu/Debian
sudo apt-get update
sudo apt-get install libarmadillo-dev

# Fedora/RHEL
sudo dnf install armadillo-devel
```

**OpenMP**:
```bash
# macOS (Homebrew)
brew install libomp

# Usually included with GCC on Linux
# For Ubuntu/Debian: sudo apt-get install libomp-dev
```

#### R Package Dependencies
```r
install.packages(c("Rcpp", "RcppArmadillo", "devtools"))
```

### Building and Installing

#### Option 1: Install from Source (Recommended)

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Navigate to the package directory
setwd("/path/to/HiSpa/R-package")

# Install the package
devtools::install()
```

#### Option 2: Using R CMD

```bash
# In terminal, from HiSpa/R-package directory
cd /path/to/HiSpa/R-package

# Build the package
R CMD build .

# Install the package
R CMD INSTALL HiSpaR_0.99.0.tar.gz
```

### Platform-Specific Notes

#### macOS
If you encounter OpenMP errors, create or edit `~/.R/Makevars`:
```bash
PKG_CXXFLAGS = -Xclang -fopenmp
PKG_LIBS = -lomp
```

#### Linux
Usually works out of the box if OpenMP is installed with GCC.

#### Windows
1. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
2. Install Armadillo (more complex on Windows, may require pre-compiled libraries)

## Quick Start

```r
library(HiSpaR)

# Prepare your Hi-C contact matrix file (space/tab delimited text file)
# The file should contain an n x n matrix of contact counts

# Run HiSpa analysis with file path
# All results are automatically saved to the output directory
hispa_analyze(
  input_file = "contact_matrix.txt",
  output_dir = "output",
  mcmc_iterations = 6000,
  mcmc_burn_in = 1000,
  num_clusters = 0,  # Auto-detect
  verbose = TRUE
)

# Results are saved in output directory:
# - output/final_positions.txt (final 3D structure)
# - output/initial_positions.txt (starting positions)
# - output/log_likelihood_trace.txt (MCMC convergence)
# - output/block_timings.txt (performance metrics)
# - output/mcmc_log.txt (detailed log)

# Read and analyze results
final_positions <- read.table("output/final_positions.txt")
ll_trace <- scan("output/log_likelihood_trace.txt")

# Plot convergence
plot(ll_trace, type = "l",
     xlab = "Iteration", ylab = "Log-Likelihood")
```

## Main Functions

### `hispa_analyze()`
Run complete HiSpa MCMC analysis on a Hi-C contact matrix following the exact workflow from the C++ implementation.

**Key Parameters:**
- `input_file`: Path to contact matrix file (space/tab delimited text format)
- `output_dir`: Output directory path where all results will be saved
- `mcmc_iterations`: Number of MCMC iterations (default: 6000)
- `num_clusters`: Number of clusters (default: 0 for auto-detect as sqrt(n))
- `mcmc_burn_in`: Burn-in period (default: 0)
- `mcmc_initial_sd`: Initial proposal SD (default: 0.1)
- `mcmc_sd_floor`: Minimum SD (default: 0.0001)
- `mcmc_sd_ceil`: Maximum SD (default: 0.3)
- `use_cluster_init`: Use cluster-based initialization (default: FALSE)
- `cluster_init_iterations`: Iterations for cluster MCMC (default: 1000)
- `cluster_initial_sd`: SD for cluster initialization (default: 0.1)
- `save_samples`: Save MCMC trace (default: FALSE)
- `sample_interval`: Save interval (default: 50)
- `verbose`: Enable detailed output (default: TRUE)

**Output Files:** All results are automatically saved to `output_dir`:
- `final_positions.txt` - Final 3D coordinates (n × 3 matrix)
- `initial_positions.txt` - Initial positions before MCMC
- `log_likelihood_trace.txt` - MCMC log-likelihood values
- `block_timings.txt` - Computation time per MCMC block
- `mcmc_log.txt` - Detailed analysis log

**Workflow:**
1. Load contact matrix from input file
2. Assign loci to clusters (k-means)
3. Build cluster relationships
4. Initialize structure (random or cluster-based)
5. Assemble global structure
6. Run main MCMC algorithm
7. Save all results to output directory

## Examples

See the `examples/` directory for complete workflows:
- **basic_analysis.R**: Simple Hi-C analysis with random initialization
- **with_priors.R**: Using cluster-based initialization
- **convolution_example.R**: Comparing different parameter settings

### Example: Basic Analysis

```r
library(HiSpaR)

# Generate synthetic data and save to file
n <- 50
contact_mat <- matrix(rpois(n*n, lambda = 10), n, n)
contact_mat <- (contact_mat + t(contact_mat)) / 2
write.table(contact_mat, "contact_matrix.txt", 
            row.names = FALSE, col.names = FALSE)

# Run analysis - results saved automatically
hispa_analyze(
  input_file = "contact_matrix.txt",
  output_dir = "example_output",
  mcmc_iterations = 4000,
  mcmc_burn_in = 1000
)

# Read results from output directory
final_positions <- as.matrix(read.table("example_output/final_positions.txt"))
ll_trace <- scan("example_output/log_likelihood_trace.txt")

# Plot convergence
plot(ll_trace, type = "l",
     xlab = "Iteration", ylab = "Log-Likelihood")

# 3D visualization (requires rgl)
if (require("rgl")) {
  plot3d(final_positions, col = rainbow(n), size = 5)
}

```

### Example: Cluster Initialization

```r
# Use cluster-based initialization for better convergence
hispa_analyze(
  input_file = "contact_matrix.txt",
  output_dir = "output",
  mcmc_iterations = 6000,
  mcmc_burn_in = 1000,
  use_cluster_init = TRUE,
  cluster_init_iterations = 1000,
  num_clusters = 5
)

# Results are saved in output/
```

## Documentation

After installation, access help documentation:
```r
?hispa_analyze
```

## Performance Tips

1. **Start with fewer iterations** for testing (e.g., 1000-2000)
2. **Use clustering** for large matrices (>200 loci)
3. **Enable convolution** for noisy data: `use_convoluted_sampling = TRUE`
4. **Parallel processing**: OpenMP automatically uses multiple cores
5. **Save samples sparingly**: Only when needed for diagnostics

## Troubleshooting

### "Armadillo not found" error
Ensure Armadillo is installed and can be found by pkg-config:
```bash
pkg-config --modversion armadillo
```

### OpenMP linking errors on macOS
Add to `~/.R/Makevars`:
```
PKG_CXXFLAGS = -Xclang -fopenmp
PKG_LIBS = -lomp
LDFLAGS = -L/usr/local/opt/libomp/lib
CPPFLAGS = -I/usr/local/opt/libomp/include
```

### Compilation errors
1. Ensure C++17 support: `clang++ --version` or `g++ --version`
2. Update R packages: `update.packages()`
3. Reinstall Rcpp: `install.packages("Rcpp", type = "source")`

## Citation

If you use HiSpaR in your research, please cite:

```
[Citation information to be added]
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

- **Issues**: https://github.com/masterStormtrooper/HiSpa/issues
- **Email**: [your contact]

## Acknowledgments

HiSpaR wraps the HiSpa C++ library developed by [authors].
