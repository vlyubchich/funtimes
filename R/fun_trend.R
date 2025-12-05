# Internal Helper Functions for Trend Testing and Estimation
#
# This file contains internal (non-exported) helper functions for trend testing
# and estimation methods. These functions are used by notrend_test(), wavk_test(),
# sync_test(), and related trend analysis methods.
#
# Functions:
#   - mann_kendall_tau(): Computes Kendall's tau for Mann-Kendall trend test
#
# @keywords internal
# @noRd


##### Mann-Kendall Trend Test Statistic #####
# Internal implementation to replace Kendall::MannKendall() dependency
# Computes Kendall's tau correlation between time series and time index
#
# @param x numeric vector (time series)
# @return numeric scalar, Kendall's tau (from -1 to 1)
mann_kendall_tau <- function(x) {
    if (length(x) < 2) stop("The time series must have at least 2 observations")
    
    # Use base R's C-optimized cor() function
    # Comparing 'x' against 'seq_along(x)' (time) is equivalent to Mann-Kendall
    stats::cor(x, seq_along(x), method = "kendall")
}
