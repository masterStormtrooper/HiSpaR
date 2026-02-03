# HiSpaR - R Package for Hi-C Spatial Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HiSpaR is an R package that provides R bindings for the HiSpa C++ library, enabling hierarchical Bayesian inference of 3D chromatin structures from Hi-C contact matrices.

## Features

- 🧬 **Hierarchical Bayesian modeling** of 3D chromatin structure
- 🎯 **MCMC sampling** with adaptive proposals for efficient inference
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


**Workflow:**
1. Load contact matrix from input file / matrix
2. Assign loci to clusters (k-means)
3. Build cluster relationships
4. Initialize structure (random or cluster-based)
5. Assemble global structure
6. Run main MCMC algorithm
7. Save all results to output directory

### Example: Basic Analysis

```r
library(HiSpaR)

# loads a hic matrix, from Su et al., 2020
data(su1_contact_mat)

# Run analysis - results saved automatically
results <- hispa_analyze(
  su1_contact_mat,
  output_dir = "su1",
  mcmc_iterations = 4000,
  mcmc_burn_in = 1000,
  use_cluster_init = TRUE # strongly recommended
)

# Read results from output directory
final_pos <- results$positions

# 3D visualization (requires plotly)
library(plotly)
plot_ly(x = final_pos[,1], y = final_pos[,2], z = final_pos[,3], type = 'scatter3d', mode = 'markers+lines')
```

## Documentation

After installation, access help documentation:
```r
?hispa_analyze
```

## Performance Tips

1. **Start with fewer iterations** for testing (e.g., 1000-2000)
2. **Use clustering** for large matrices (>200 loci)
3. **Parallel processing**: OpenMP automatically uses multiple cores
4. **Save samples sparingly**: Only when needed for diagnostics

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



## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

- **Issues**: https://github.com/masterStormtrooper/HiSpa/issues
- **Email**: lyc22@mails.tsinghua.edu.cn

## Acknowledgments

HiSpaR wraps the HiSpa C++ library developed by Yingcheng Luo.
