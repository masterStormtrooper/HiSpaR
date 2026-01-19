# Known Issues

## Terminal Abort After Successful Completion (macOS)

### Description
On some systems (particularly macOS with certain R configurations), the `hispa_analyze()` function may cause the R session to abort with exit code 134 after successful completion. This happens during R's cleanup phase, AFTER all output files have been correctly generated.

### Impact
- **Data is NOT affected**: All output files (`final_positions.txt`, `initial_positions.txt`, etc.) are correctly created before the abort
- **Results are valid**: The MCMC algorithm completes successfully
- **Only affects**: Interactive R sessions when called directly

### Workaround
The issue does not occur when:
1. Running in non-interactive mode (R CMD BATCH, Rscript)
2. Building vignettes with knitr
3. Running automated tests

### Root Cause
This is a known interaction between the spdlog logging library and R's memory management system. The standalone C++ version does not have this issue.

### Status
We are investigating a permanent fix. The abort does not affect the scientific validity of the results.
