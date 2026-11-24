# Reference values in this file were calculated independently of the package
# code and without the noncentral version of pt(): the noncentral t
# probability was obtained by integrating over its normal component.

expect_close <- function(value, reference, tolerance=1e-6) {
     expect_lt(abs(value - reference), tolerance)
}

test_that("d confidence interval matches independent values", {
     ci <- get_d_ci(d=0.5, n1=20, n2=20)
     expect_close(ci$LL, -0.1330181)
     expect_close(ci$UL, 1.1266613)

     ci <- get_d_ci(d=1.2, n1=15, n2=30)
     expect_close(ci$LL, 0.5246049)
     expect_close(ci$UL, 1.8633919)

     ci <- get_d_ci(d=0, n1=10, n2=10)
     expect_close(ci$LL, -0.8765225)
     expect_close(ci$UL, 0.8765225)

     ci <- get_d_ci(d=2.5, n1=8, n2=8)
     expect_close(ci$LL, 1.1343061)
     expect_close(ci$UL, 3.8183402)

     ci <- get_d_ci(d=0.65, n1=50, n2=50)
     expect_close(ci$LL, 0.2460344)
     expect_close(ci$UL, 1.0508146)
})

test_that("d confidence interval matches independent values at other confidence levels", {
     ci <- get_d_ci(d=0.5, n1=20, n2=20, conf_level=.90)
     expect_close(ci$LL, -0.0317965)
     expect_close(ci$UL, 1.0253598)

     ci <- get_d_ci(d=0.5, n1=20, n2=20, conf_level=.99)
     expect_close(ci$LL, -0.3308299)
     expect_close(ci$UL, 1.3246694)

     ci <- get_d_ci(d=-1.2, n1=15, n2=30, conf_level=.80)
     expect_close(ci$LL, -1.6309545)
     expect_close(ci$UL, -0.7555657)

     ci <- get_d_ci(d=0.3, n1=5, n2=7, conf_level=.50)
     expect_close(ci$LL, -0.1048697)
     expect_close(ci$UL, 0.6900437)
})

test_that("d confidence interval is correct when the noncentrality is beyond the range of pt()", {
     #t = 40
     ci <- get_d_ci(d=0.8, n1=5000, n2=5000)
     expect_close(ci$LL, 0.7592440)
     expect_close(ci$UL, 0.8407188)

     ci <- get_d_ci(d=0.8, n1=5000, n2=5000, conf_level=.99)
     expect_close(ci$LL, 0.7464448)
     expect_close(ci$UL, 0.8535207)

     #t = 38.7
     ci <- get_d_ci(d=5, n1=120, n2=120)
     expect_close(ci$LL, 4.4832367)
     expect_close(ci$UL, 5.5138340)

     #t = 80
     ci <- get_d_ci(d=1.6, n1=5000, n2=5000)
     expect_close(ci$LL, 1.5549317)
     expect_close(ci$UL, 1.6450067)

     #t = 37.5, upper limit beyond 37.62
     ci <- get_d_ci(d=3.75, n1=200, n2=200)
     expect_close(ci$LL, 3.4231396)
     expect_close(ci$UL, 4.0750036)

     #a negative d mirrors the limits
     ci <- get_d_ci(d=-3.75, n1=200, n2=200)
     expect_close(ci$LL, -4.0750036)
     expect_close(ci$UL, -3.4231396)

     #t = 39.1 with a very large df, where pt() returns 1
     ci <- get_d_ci(d=0.35, n1=25000, n2=25000)
     expect_close(ci$LL, 0.3323341)
     expect_close(ci$UL, 0.3676624)

     ci <- get_d_ci(d=10, n1=50, n2=50)
     expect_close(ci$LL, 8.5454194)
     expect_close(ci$UL, 11.4487060)

     #small groups, t = 20, upper limit beyond 37.62
     ci <- get_d_ci(d=20, n1=2, n2=2)
     expect_close(ci$LL, 3.0291006)
     expect_close(ci$UL, 38.4958055)
})

test_that("pt() and numerical integration agree where both can be used", {
     pt_prob <- function(t_value, df, delta) {
          suppressWarnings(stats::pt(t_value, df, ncp=delta))
     }
     limits_pt       <- solve_nct_limits(t_value=20, df=98, conf_level=.95, prob_t_below=pt_prob)
     limits_integral <- solve_nct_limits(t_value=20, df=98, conf_level=.95, prob_t_below=prob_t_below_integral)
     expect_close(limits_pt$LL, limits_integral$LL)
     expect_close(limits_pt$UL, limits_integral$UL)

     limits_pt       <- solve_nct_limits(t_value=2.2, df=6, conf_level=.95, prob_t_below=pt_prob)
     limits_integral <- solve_nct_limits(t_value=2.2, df=6, conf_level=.95, prob_t_below=prob_t_below_integral)
     expect_close(limits_pt$LL, limits_integral$LL)
     expect_close(limits_pt$UL, limits_integral$UL)
})

test_that("Numerical integration handles a t of zero and a negative t", {
     expect_equal(prob_t_below_integral(t_value=0, df=10, delta=0.7), stats::pnorm(-0.7))

     #independent values (integral over the normal component)
     expect_close(prob_t_below_integral(t_value=2, df=10, delta=1), 0.8076115625, tolerance=1e-8)
     expect_close(prob_t_below_integral(t_value=-2, df=10, delta=-1), 0.1923884375, tolerance=1e-8)
     expect_close(prob_t_below_integral(t_value=45, df=198, delta=44), 0.6465719464, tolerance=1e-8)
     expect_close(prob_t_below_integral(t_value=-45, df=198, delta=-44), 0.3534280536, tolerance=1e-8)

     #limits for a negative t mirror those for a positive t
     limits_positive <- get_nct_limits(t_value=40, df=9998, conf_level=.95)
     limits_negative <- get_nct_limits(t_value=-40, df=9998, conf_level=.95)
     expect_close(limits_negative$LL, -limits_positive$UL, tolerance=1e-8)
     expect_close(limits_negative$UL, -limits_positive$LL, tolerance=1e-8)
})

test_that("Large group sizes stored as integers do not overflow", {
     ci_integer <- get_d_ci(d=0.3, n1=50000L, n2=50000L)
     ci_double  <- get_d_ci(d=0.3, n1=50000, n2=50000)
     expect_equal(ci_integer$LL, ci_double$LL)
     expect_equal(ci_integer$UL, ci_double$UL)
     expect_true(ci_integer$LL < 0.3)
     expect_true(ci_integer$UL > 0.3)
})

test_that("Invalid input gives an error", {
     expect_error(get_d_ci(d=Inf, n1=5, n2=5))
     expect_error(get_d_ci(d=NaN, n1=5, n2=5))
     expect_error(get_d_ci(d=NA, n1=5, n2=5))
     expect_error(get_d_ci(d=0.5, n1=20, n2=20, conf_level=95))
     expect_error(get_d_ci(d=0.5, n1=20, n2=20, conf_level=1))
     expect_error(get_d_ci(d=0.5, n1=20, n2=20, conf_level=0))
     expect_error(get_d_ci(d=0.5, n1=1, n2=1))
     expect_error(get_nct_limits(t_value=2, df=0, conf_level=.95))
     expect_error(get_nct_limits(t_value=Inf, df=10, conf_level=.95))
})
