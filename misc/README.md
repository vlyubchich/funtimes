# Misc Directory

This directory contains experimental, deprecated, or archived code that is not part of the main package.

## Contents

### Archived/Removed Functions

- **TopoCBN.R** - Topological clustering function using Betti numbers. This function was removed from the package in version 8.2 due to removal of the `TDA` package dependency. The implementation is preserved here for reference.

### Experimental/Incomplete Code

- **CRAD.R** and **CRAD.py** - Experimental clustering using robust anomaly detection. This appears to be incomplete/under development and is not currently integrated into the package.
- **CRAD.Rd** - Draft documentation for CRAD function.

### Old Versions

- **causality_predVAR_v0.R** - Previous version of the `causality_predVAR()` function, retained for reference during development.

## Recommendations

1. **CRAD files** - Consider either:
   - Completing the implementation and integrating into the package, or
   - Removing if the approach was abandoned

2. **TopoCBN.R** - Keep for archival purposes since it was a documented feature removed for dependency reasons. Users may still want to access it.

3. **causality_predVAR_v0.R** - Can be removed once the new version is stable and well-tested.
