# Reference values in this file were calculated independently of the package
# code. For pi.d the confidence limits came from the noncentral t distribution
# (obtained by integrating over its normal component, without the noncentral
# version of pt()) and were combined with Zou (2007), Equation 15. For pi.r
# and pi.m the formulas were written out again.

expect_close <- function(value, reference, tolerance=1e-6) {
     expect_lt(abs(value - reference), tolerance)
}

expect_intervals <- function(x, ci_LL, ci_UL, pi_LL, pi_UL) {
     expect_close(x$lower_confidence_interval, ci_LL)
     expect_close(x$upper_confidence_interval, ci_UL)
     expect_close(x$lower_prediction_interval, pi_LL)
     expect_close(x$upper_prediction_interval, pi_UL)
}

test_that("pi.d matches independent values", {
     x <- pi.d(d=.65, n1=50, n2=50, rep.n1=100, rep.n2=100)
     expect_intervals(x, 0.2460344, 1.0508146, 0.1563735, 1.1419550)

     x <- pi.d(d=-.8, n1=30, n2=40, rep.n1=60, rep.n2=80, prob.level=.90)
     expect_intervals(x, -1.2101160, -0.3842633, -1.3041634, -0.2928479)

     x <- pi.d(d=0, n1=10, n2=10, rep.n1=10, rep.n2=10)
     expect_intervals(x, -0.8765225, 0.8765225, -1.2395901, 1.2395901)

     x <- pi.d(d=2, n1=25, n2=25, rep.n1=40, rep.n2=40, prob.level=.99)
     expect_intervals(x, 1.0987782, 2.8939308, 0.8553551, 3.1417719)
})

test_that("pi.d is correct when the replication is large enough to be beyond the range of pt()", {
     #replication t = 40
     x <- pi.d(d=.8, n1=100, n2=100, rep.n1=5000, rep.n2=5000)
     expect_intervals(x, 0.5109134, 1.0872110, 0.5080598, 1.0900882)
})

test_that("pi.d uses the original cell sizes when replication cell sizes are not given", {
     x <- pi.d(d=.5, n1=20, n2=20)
     expect_intervals(x, -0.1330181, 1.1266613, -0.3907392, 1.3907392)
     expect_equal(x$replication_N1, 20)
     expect_equal(x$replication_N2, 20)
})

test_that("pi.d uses each original cell size for its own replication cell", {
     x <- pi.d(d=.5, n1=20, n2=80)
     expect_intervals(x, 0.0037974, 0.9937033, -0.1999714, 1.1999714)
     expect_equal(x$replication_N1, 20)
     expect_equal(x$replication_N2, 80)
})

test_that("pi.d matches independent values when the limits are far from d", {
     x <- pi.d(d=0, n1=10, n2=10, prob.level=.99)
     expect_intervals(x, -1.1519459, 1.1519459, -1.6290975, 1.6290975)

     x <- pi.d(d=.1, n1=20, n2=20, prob.level=.99)
     expect_intervals(x, -0.7157349, 0.9144269, -1.0526988, 1.2526988)
})

test_that("pi.d reports the probability level that was requested", {
     x <- pi.d(d=.5, n1=20, n2=20, prob.level=.90)
     expect_equal(x$prob_level, .90)
     expect_intervals(x, -0.0317965, 1.0253598, -0.2475362, 1.2475362)
     expect_output(print(x), "90% CI[-0.03, 1.03]", fixed=TRUE)
     expect_output(print(x), "90% PI\\[-0\\.25, ?1\\.25\\]")
})

test_that("pi.d gives no warnings", {
     expect_warning(pi.d(d=2, n1=25, n2=25, rep.n1=40, rep.n2=40, prob.level=.99), NA)
     expect_warning(pi.d(d=.8, n1=5000, n2=5000), NA)
     expect_warning(pi.d(d=3.75, n1=200, n2=200), NA)
})

test_that("pi.d returns and prints the expected output", {
     x <- pi.d(d=.65, n1=50, n2=50, rep.n1=100, rep.n2=100)
     expect_s3_class(x, "d_prediction_interval")
     expect_equal(x$original_d, .65)
     expect_equal(x$prob_level, .95)
     expect_output(print(x), "Original study: d = 0.65, N1 = 50, N2 = 50, 95% CI[0.25, 1.05]", fixed=TRUE)
     expect_output(print(x), "Prediction interval: 95% PI[0.16, 1.14]", fixed=TRUE)
})

test_that("The interpretation text reads as intended", {
     x <- pi.d(d=.65, n1=50, n2=50, rep.n1=100, rep.n2=100)
     expect_equal(x$method_text, paste(
          "The original d-value is 0.65, with a prediction interval of 95% PI[0.16, 1.14]",
          "based on replication cell sizes of N1 = 100 and N2 = 100.",
          "If the replication d-value differs from the original d-value only due to sampling error,",
          "there is a 95% chance the replication result will fall in this interval.",
          "If the replication d-value falls outside of this interval,",
          "factors beyond sampling error are likely also responsible for the difference."))
     expect_equal(x$ri_text, "The prediction interval is 95% PI[0.16, 1.14].")

     x <- pi.r(r=.35, n=100, rep.n=200)
     expect_equal(x$method_text, paste(
          "The original correlation is 0.35, with a prediction interval of 95% PI[0.13, 0.56]",
          "based on a replication sample size of N = 200.",
          "If the replication correlation differs from the original correlation only due to sampling error,",
          "there is a 95% chance the replication result will fall in this interval.",
          "If the replication correlation falls outside of this interval,",
          "factors beyond sampling error are likely also responsible for the difference."))

     x <- pi.m(M=2.53, SD=1.02, n=40, rep.n=80)
     expect_equal(x$method_text, paste(
          "The original mean is 2.53, with a prediction interval of 95% PI[2.13, 2.93]",
          "based on a replication sample size of N = 80.",
          "If the replication mean differs from the original mean only due to sampling error,",
          "there is a 95% chance the replication result will fall in this interval.",
          "If the replication mean falls outside of this interval,",
          "factors beyond sampling error are likely also responsible for the difference."))
})

test_that("The printed output has the intended layout and ends with a line break", {
     x <- pi.r(r=.35, n=100, rep.n=200)
     output_file <- tempfile()
     capture.output(print(x), file=output_file)
     printed <- readChar(output_file, file.info(output_file)$size)
     unlink(output_file)
     #a file written on Windows has a carriage return before each line break
     printed <- gsub("\r\n", "\n", printed, fixed=TRUE)
     expect_equal(printed, paste0(
          "\nOriginal study: r = 0.35, N = 100, 95% CI[0.16, 0.51]",
          "\nReplication study: N = 200",
          "\nPrediction interval: 95% PI[0.13, 0.56]\n",
          "\nInterpretation:\n", x$method_text, "\n"))
})

test_that("Every accepted probability level can be printed", {
     for (level in c(.07, .14, .28, .29, .55, .56, .57, .58, .90, .95, .99)) {
          percent <- sprintf("%d%% CI[", round(level*100))
          expect_output(print(pi.d(d=.5, n1=20, n2=20, prob.level=level)), percent, fixed=TRUE)
          expect_output(print(pi.r(r=.35, n=100, prob.level=level)), percent, fixed=TRUE)
          expect_output(print(pi.m(M=2.53, SD=1.02, n=40, prob.level=level)), percent, fixed=TRUE)
     }
     set.seed(1)
     expect_output(print(pi.d.demo(number.trials=5, prob.level=.57)), "57% Prediction interval capture percentage", fixed=TRUE)
     expect_output(print(pi.r.demo(number.trials=5, prob.level=.57)), "57% Prediction interval capture percentage", fixed=TRUE)
     expect_output(print(pi.m.demo(number.trials=5, prob.level=.57)), "57% Prediction interval capture percentage", fixed=TRUE)
})

test_that("Sample sizes must be whole numbers", {
     message_expected <- "Sample sizes must be whole numbers"
     expect_error(pi.d(d=.5, n1=20.5, n2=20), message_expected)
     expect_error(pi.d(d=.5, n1=20, n2=20, rep.n1=30.5, rep.n2=30), message_expected)
     expect_error(pi.d(d=.5, n1=NA, n2=20), message_expected)
     expect_error(pi.r(r=.35, n=100.5), message_expected)
     expect_error(pi.r(r=.35, n=100, rep.n=150.5), message_expected)
     expect_error(pi.m(M=2.53, SD=1.02, n=40.5), message_expected)
     expect_error(pi.m(M=2.53, SD=1.02, n=40, rep.n=80.5), message_expected)
     expect_error(pi.d.demo(n1=20.5, n2=20, number.trials=5), message_expected)
     expect_error(pi.r.demo(n=100.5, number.trials=5), message_expected)
     expect_error(pi.m.demo(n=40.5, number.trials=5), message_expected)

     #whole numbers are accepted whether stored as integers or not
     expect_equal(pi.d(d=.5, n1=20L, n2=20L)$lower_prediction_interval, pi.d(d=.5, n1=20, n2=20)$lower_prediction_interval)
     expect_equal(pi.r(r=.35, n=100L)$lower_prediction_interval, pi.r(r=.35, n=100)$lower_prediction_interval)
})

test_that("pi.d, pi.r and pi.m describe the interval as a prediction interval (PI)", {
     intervals <- list(pi.d(d=.65, n1=50, n2=50, rep.n1=100, rep.n2=100),
                       pi.r(r=.35, n=100, rep.n=200),
                       pi.m(M=2.53, SD=1.02, n=40, rep.n=80))
     for (x in intervals) {
          for (txt in c(x$method_text, x$ri_text, paste(capture.output(print(x)), collapse=" "))) {
               expect_true(grepl("prediction interval", txt, ignore.case=TRUE))
               expect_true(grepl("% PI[", txt, fixed=TRUE))
               expect_false(grepl("replication interval", txt, ignore.case=TRUE))
               expect_false(grepl("RI[", txt, fixed=TRUE))
          }
     }
     expect_true(grepl("95% PI[0.16, 1.14]", intervals[[1]]$method_text, fixed=TRUE))
})

test_that("pi.d does not change the warn option", {
     warn_before <- getOption("warn")
     x <- pi.d(d=.65, n1=50, n2=50, rep.n1=100, rep.n2=100)
     x <- pi.d(d=.8, n1=5000, n2=5000)
     expect_equal(getOption("warn"), warn_before)
})

test_that("prob.level must be a value with at most two decimals from .01 to .99", {
     accepted <- c(.01, .50, .57, .90, .95, .99, 1 - .05)
     for (level in accepted) {
          expect_equal(pi.d(d=.5, n1=20, n2=20, prob.level=level)$prob_level, level)
          expect_equal(pi.r(r=.35, n=100, prob.level=level)$prob_level, level)
          expect_equal(pi.m(M=2.53, SD=1.02, n=40, prob.level=level)$prob_level, level)
     }

     rejected <- list(.999, .995, .975, .683, .9999, 1 - 1e-12, 1, 0, .001, 95, -.95, NA, Inf, "a", c(.90, .95), NULL)
     for (level in rejected) {
          expect_error(pi.d(d=.5, n1=20, n2=20, prob.level=level), "prob.level must be")
          expect_error(pi.r(r=.35, n=100, prob.level=level), "prob.level must be")
          expect_error(pi.m(M=2.53, SD=1.02, n=40, prob.level=level), "prob.level must be")
     }

     #the demos stop before any trials are run
     expect_error(pi.d.demo(number.trials=5, prob.level=.999), "prob.level must be")
     expect_error(pi.r.demo(number.trials=5, prob.level=.999), "prob.level must be")
     expect_error(pi.m.demo(number.trials=5, prob.level=.999), "prob.level must be")
})

test_that("pi.d stops when the interval cannot be calculated", {
     expect_error(pi.d(d=Inf, n1=50, n2=50))
     expect_error(pi.d(d=.65, n1=50, n2=50, prob.level=95))
     expect_error(pi.d(d=.65, n1=1, n2=1))
})

test_that("pi.r matches independent values", {
     x <- pi.r(r=.35, n=100, rep.n=200)
     expect_intervals(x, 0.1649195, 0.5112702, 0.1315151, 0.5558678)
     expect_s3_class(x, "r_prediction_interval")
})

test_that("pi.m matches independent values", {
     x <- pi.m(M=2.53, SD=1.02, n=40, rep.n=80)
     expect_intervals(x, 2.2037882, 2.8562118, 2.1304737, 2.9295263)
     expect_s3_class(x, "M_prediction_interval")

     x_var <- pi.m(M=2.53, VAR=1.02^2, n=40, rep.n=80)
     expect_equal(x_var$lower_prediction_interval, x$lower_prediction_interval)
     expect_equal(x_var$upper_prediction_interval, x$upper_prediction_interval)
})

test_that("pi.d.demo runs and reports each trial", {
     set.seed(1)
     demo <- pi.d.demo(n1=50, n2=50, rep.n1=100, rep.n2=100, pop.d=.50, number.trials=20)
     expect_s3_class(demo, "replication_demo_d")
     trials <- demo$results_each_trial
     expect_equal(nrow(trials), 20)
     expect_true(all(trials$ci.LL < trials$d & trials$d < trials$ci.UL))
     expect_true(all(trials$pi.LL < trials$ci.LL & trials$ci.UL < trials$pi.UL))
     expect_equal(demo$in_prediction_interval_count, sum(trials$rep.d.in.pi))
     expect_equal(demo$percent_in_ri, 100 * sum(trials$rep.d.in.pi) / 20)
})

# The illustrative trials table as printed, read back into a data frame
printed_trials <- function(demo) {
     old_options <- options(width=200)
     on.exit(options(old_options))
     printed <- capture.output(print(demo))
     header <- grep("ci\\.LL +ci\\.UL", printed)[1]
     read.table(text=printed[header:(header+5)], header=TRUE)
}

test_that("The demos print each trial's own confidence and prediction limits", {
     set.seed(1)
     demos <- list(pi.m.demo(n=150, mu=0, sigma=1, number.trials=20),
                   pi.r.demo(n=100, rho=.50, number.trials=20),
                   pi.d.demo(n1=50, n2=50, rep.n1=100, rep.n2=100, pop.d=.50, number.trials=20))
     #every demo prints the limits to two decimals
     for (demo in demos) {
          shown  <- printed_trials(demo)
          stored <- demo$results_each_trial[1:5,]
          expect_equal(nrow(shown), 5)
          #the upper limit column once repeated the lower limit
          expect_true(all(shown$ci.LL < shown$ci.UL))
          expect_equal(shown$ci.LL, round(stored$ci.LL, 2))
          expect_equal(shown$ci.UL, round(stored$ci.UL, 2))
          expect_equal(shown$pi.LL, round(stored$pi.LL, 2))
          expect_equal(shown$pi.UL, round(stored$pi.UL, 2))
     }
})

test_that("The notes under the illustrative trials label the columns correctly", {
     set.seed(1)
     expect_output(print(pi.r.demo(n=100, rho=.50, number.trials=20)),
                   "pi.LL = lower-limit prediction interval, pi.UL = upper-limit prediction interval", fixed=TRUE)
     expect_output(print(pi.d.demo(n1=50, n2=50, pop.d=.50, number.trials=20)),
                   "pi.LL = lower-limit prediction interval, pi.UL = upper-limit prediction interval", fixed=TRUE)
     expect_output(print(pi.m.demo(n=150, mu=0, sigma=1, number.trials=20)),
                   "pi.LL = lower-limit prediction interval, pi.UL = upper-limit prediction interval, rep.M = replication mean.", fixed=TRUE)
})

test_that("pi.r.demo and pi.m.demo run", {
     set.seed(1)
     demo <- pi.r.demo(n=100, rho=.50, number.trials=20)
     expect_s3_class(demo, "replication_demo_r")
     expect_equal(nrow(demo$results_each_trial), 20)

     demo <- pi.m.demo(n=150, mu=0, sigma=1, number.trials=20)
     expect_s3_class(demo, "replication_demo_M")
     expect_equal(nrow(demo$results_each_trial), 20)
})
