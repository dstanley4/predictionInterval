# Confidence limits for d, found by inverting the noncentral t distribution:
# each limit is the value of the noncentrality parameter at which the observed
# t sits at the required percentile (Steiger & Fouladi, 1997; Cumming & Finch,
# 2001). The noncentral probabilities come from pt() in the stats package.
#
# pt() is accurate only while the noncentrality parameter is no larger than
# 37.62 in absolute value, and it returns 1 when the observed t is beyond
# about 38.6 with a very large df. Whenever the observed t or one of its
# limits reaches that region, the same noncentral t probability is obtained
# by numerical integration of its definition instead (prob_t_below_integral).

# pt() is trusted only while the observed t and both limits are inside this
nct_trusted_range <- 37.5


check_conf_level <- function(conf_level) {
     if (length(conf_level) != 1 || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
          stop("The confidence level must be a single number between 0 and 1.", call. = FALSE)
     }
}


# Limits for the noncentrality parameter of an observed t.
get_nct_limits <- function(t_value, df, conf_level) {
     check_conf_level(conf_level)
     if (length(t_value) != 1 || !is.finite(t_value)) {
          stop("The t-value must be a single finite number.", call. = FALSE)
     }
     if (length(df) != 1 || !is.finite(df) || df <= 0) {
          stop("The degrees of freedom must be a positive, finite number.", call. = FALSE)
     }
     t_value <- as.numeric(t_value)
     df <- as.numeric(df)

     # pt() is used when the observed t and both limits are inside the range
     # where it can be trusted. Its warnings are not passed on because the
     # limits are recalculated by integration whenever that range is left.
     limits <- tryCatch({
          pt_trusted <- abs(t_value) <= nct_trusted_range
          if (pt_trusted) {
               limits <- solve_nct_limits(t_value, df, conf_level, prob_t_below = function(t_value, df, delta) {
                    suppressWarnings(stats::pt(t_value, df, ncp = delta))
               })
               pt_trusted <- max(abs(c(limits$LL, limits$UL))) <= nct_trusted_range
          }
          if (!pt_trusted) {
               limits <- solve_nct_limits(t_value, df, conf_level, prob_t_below = prob_t_below_integral)
          }
          limits
     }, error = function(e) {
          stop("A confidence interval could not be calculated for this t-value.", call. = FALSE)
     })

     return(limits)
}


# Probability of a t at or below t_value when the noncentrality is delta,
# by numerical integration. A noncentral t is (Z + delta) / S, where Z is
# standard normal and S is the square root of an independent chi-square
# divided by its df, so the probability is pnorm(t * S - delta) averaged over
# the distribution of S. Used where pt() cannot be trusted.
prob_t_below_integral <- function(t_value, df, delta) {
     if (t_value < 0) {
          return(1 - prob_t_below_integral(-t_value, df, -delta))
     }
     if (t_value == 0) {
          return(stats::pnorm(-delta))
     }

     S_lowest  <- sqrt(stats::qchisq(1e-18, df) / df)
     S_highest <- sqrt(stats::qchisq(1e-18, df, lower.tail = FALSE) / df)

     # pnorm(t * S - delta) is 0 for S below this range and 1 for S above it
     lower <- max(S_lowest,  (delta - 8.3) / t_value)
     upper <- min(S_highest, (delta + 8.3) / t_value)
     if (upper <= S_lowest) {
          return(1)
     }
     if (lower >= S_highest) {
          return(0)
     }

     prob_S_above <- stats::pchisq(df * upper^2, df, lower.tail = FALSE)
     integrand <- function(s) {
          stats::pnorm(t_value * s - delta) * 2 * df * s * stats::dchisq(df * s^2, df)
     }
     prob_S_within <- stats::integrate(integrand, lower = lower, upper = upper,
                                       rel.tol = 1e-10, subdivisions = 1000L)$value
     min(1, max(0, prob_S_above + prob_S_within))
}


# Finds both limits, given a function for the noncentral t probability
solve_nct_limits <- function(t_value, df, conf_level, prob_t_below) {
     tail_prob <- (1 - conf_level) / 2

     # Probability of a t at or below the observed one; falls as delta rises
     prob_below <- function(delta) {
          prob_t_below(t_value, df, delta)
     }

     find_limit <- function(target) {
          step <- max(1, abs(t_value))
          lower <- t_value - step
          upper <- t_value + step
          for (attempt in 1:200) {
               prob_lower <- prob_below(lower)
               prob_upper <- prob_below(upper)
               if (is.na(prob_lower) || is.na(prob_upper)) {
                    stop("A confidence limit could not be calculated for this t-value.", call. = FALSE)
               }
               if (prob_lower >= target && prob_upper <= target) {
                    break
               }
               if (prob_lower < target) {
                    lower <- lower - step
               }
               if (prob_upper > target) {
                    upper <- upper + step
               }
          }
          stats::uniroot(function(delta) prob_below(delta) - target,
                         lower = lower, upper = upper,
                         tol = 1e-10, maxiter = 1000)$root
     }

     output <- list()
     output$LL <- find_limit(1 - tail_prob)
     output$UL <- find_limit(tail_prob)
     return(output)
}


# Confidence interval for d (standardized mean difference, two independent groups)
get_d_ci <- function(d, n1, n2, conf_level = .95) {
     if (length(d) != 1 || !is.finite(d)) {
          stop("A confidence interval for d could not be calculated because d is not a finite number.", call. = FALSE)
     }
     n1 <- as.numeric(n1)
     n2 <- as.numeric(n2)
     d_to_t <- sqrt((n1 * n2) / (n1 + n2))
     t_value <- d * d_to_t
     t_limits <- get_nct_limits(t_value = t_value, df = n1 + n2 - 2, conf_level = conf_level)

     output <- list()
     output$LL <- t_limits$LL / d_to_t
     output$UL <- t_limits$UL / d_to_t
     return(output)
}
