
# Interval Calculation ----------------------------------------------------

#' \emph{d}-value (i.e., standardized mean difference) Prediction Interval
#' @param d Original study: Sample \emph{d}-value (standardized mean difference) created with pooled variance denominator. See formulas 4.18 and 4.19 (p.26) in Borenstein, Hedges, Higgins, & Rothstein (2009).
#' @param n1 Original study: Sample size for group 1
#' @param n2 Original study: Sample size for group 2
#' @param rep.n1 (optional) Replication study: Sample size for group 1. If not specified, n1 is used.
#' @param rep.n2 (optional) Replication study: Sample size for group 2. If not specified, n2 is used.
#' @param prob.level (optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @return The prediction interval and related statistics in list format.
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2009). \emph{Introduction to meta-analysis}. John Wiley & Sons.\cr\cr
#'Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are based on central and noncentral distributions. \emph{Educational and Psychological Measurement, 61(4)}, 532-574.
#' @examples
#' pi.d(d=.65,n1=50,n2=50,rep.n1=100,rep.n2=100)
#' @export
pi.d <- function (d,n1,n2,rep.n1=NA, rep.n2=NA,prob.level=.95) {
     prob_level <- prob.level
     if (is.na(rep.n1)) {
          rep.n1=n1
     }
     if (is.na(rep.n2)) {
          rep.n2=n2
     }

     original_d <-d
     original_N1 <- n1
     original_N2 <- n2
     replication_N1 <- rep.n1
     replication_N2 <- rep.n2

     original_CI <- ci.smd(smd=original_d,n.1=original_N1,n.2=original_N2,conf.level=prob_level)
     l1 <- original_CI$Lower.Conf.Limit.smd
     u1 <- original_CI$Upper.Conf.Limit.smd

     mod_CI <- ci.smd(smd=original_d,n.1=replication_N1,n.2=replication_N2,conf.level=prob_level)
     l2 <- mod_CI$Lower.Conf.Limit.smd
     u2 <- mod_CI$Upper.Conf.Limit.smd

     prediction_interval_values <- mas_interval(original_r=original_d,l1=l1,u1=u1,l2=l2,u2=u2)

     LL <- prediction_interval_values$LL
     UL <- prediction_interval_values$UL
     prediction_interval_metrics <- list()
     prediction_interval_metrics$original_d <- original_d
     prediction_interval_metrics$original_N1 <- original_N1
     prediction_interval_metrics$original_N2 <- original_N2
     prediction_interval_metrics$replication_N1 <- replication_N1
     prediction_interval_metrics$replication_N2 <- replication_N2
     prediction_interval_metrics$lower_prediction_interval <- LL
     prediction_interval_metrics$upper_prediction_interval <- UL
     prediction_interval_metrics$lower_confidence_interval <- l1
     prediction_interval_metrics$upper_confidence_interval <- u1
     prediction_interval_metrics$prob_level <- prob_level

     percent_level <- as.integer(round(prob_level*100))
     method_text <- get_method_text_d(original_d,LL,UL,replication_N1,replication_N2,percent_level,"d-value")
     prediction_interval_metrics$method_text <- method_text$txt_combined
     prediction_interval_metrics$ri_text <- method_text$txt_ri

     class(prediction_interval_metrics) <- "d_prediction_interval"

     return(prediction_interval_metrics)
}




#' @export
print.d_prediction_interval <- function(x,...) {
     conf_per <- x$prob_level * 100

     cat(sprintf("\nOriginal study: d = %1.2f, N1 = %d, N2 = %d, %d%% CI[%1.2f, %1.2f]",x$original_d,x$original_N1,x$original_N2,conf_per,x$lower_confidence_interval,x$upper_confidence_interval))
     cat(sprintf("\nReplication study: N1 = %d, N2 = %d",x$replication_N1,x$replication_N2))
     cat(sprintf("\nPrediction interval: %d%% PI[%1.2f,%1.2f].\n\n",conf_per,x$lower_prediction_interval,x$upper_prediction_interval))
     cat("\nInterpretation:\n")
     cat(x$method_text)

}




# Simulation --------------------------------------------------------------





#' Simulation to demonstrate the meaning of the d-value prediction interval
#' @param n1 Original study: Cell size 1
#' @param n2 Original study: Cell size 2
#' @param rep.n1 (optional) Replication study: Cell size 1. If not specified, n is used.
#' @param rep.n2 (optional) Replication study: Cell size 2. If not specified, n is used.
#' @param pop.d All samples are drawn from a common population. This specifies the population correlation.
#' @param number.trials Indicate the number of pairs of sample (original, replication) that should be used. 10,000 or higher suggested for stable results.
#' @param prob.level (optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @param bias.correction Apply bias correction formula to d-values.
#' @return The prediction interval capture percentage and related statistics in list format.
#' @examples
#' pi.d.demo(n1=50,n2=50,rep.n1=100,rep.n2=100,pop.d=.50,number.trials=10)
#' @export
pi.d.demo <- function(n1=50,n2=50,rep.n1=NA,rep.n2=NA,pop.d=.50,number.trials=10000,prob.level=.95,bias.correction = FALSE) {
     bias_correction <- bias.correction
     number_trials <- number.trials
     original_N1 <- n1
     original_N2 <- n2
     prob_level <- prob.level
     if (missing(rep.n1)) {rep.n1 <- n1}
     if (missing(rep.n2)) {rep.n2 <- n2}
     replication_N1 <- rep.n1
     replication_N2 <- rep.n2
     pop_d <- pop.d

     output <-pbapply::pbreplicate(number_trials,get_orig_rep_d(original_N1,original_N2,replication_N1,replication_N2,pop_d,prob_level,bias_correction))

     #output <- c(original_N1,original_N2,original_d,ci_LL,ci_UL,replication_N1,replication_N2,replication_d,ri_LL,ri_UL,is.in.ci,is.in.pi)

     n1 <-output[,1,]
     n2 <-output[,2,]
     d <-output[,3,]
     ci.LL <-output[,4,]
     ci.UL <-output[,5,]
     rep.n1 <-output[,6,]
     rep.n2 <-output[,7,]
     rep.d <-output[,8,]
     pi.LL <-output[,9,]
     pi.UL <-output[,10,]
     rep.d.in.ci <-as.logical(output[,11,])
     rep.d.in.pi <-as.logical(output[,12,])
     output_df <- data.frame(n1,n2,d,ci.LL,ci.UL,rep.n1,rep.n2,pi.LL,pi.UL,rep.d,rep.d.in.ci,rep.d.in.pi)


     in_prediction_interval_count <- sum(rep.d.in.pi)
     in_confidence_interval_count <- sum(rep.d.in.ci)
     percent_in_ri <- (in_prediction_interval_count/(number_trials))*100
     percent_in_ci <- (in_confidence_interval_count/(number_trials))*100
     replication_demo_output <- list()
     replication_demo_output$percent_in_ri <- percent_in_ri
     replication_demo_output$percent_in_ci <- percent_in_ci
     replication_demo_output$in_prediction_interval_count<-in_prediction_interval_count
     replication_demo_output$in_confidence_interval_count<-in_confidence_interval_count
     replication_demo_output$results_each_trial <- output_df
     replication_demo_output$pop_d <- pop_d
     replication_demo_output$original_N1 <- original_N1
     replication_demo_output$original_N2 <- original_N2
     replication_demo_output$replication_N1 <- replication_N1
     replication_demo_output$replication_N2 <- replication_N2
     replication_demo_output$prob_level <- prob_level

     class(replication_demo_output) <- "replication_demo_d"

     return(replication_demo_output)
}

#' @export
print.replication_demo_d <- function(x,...) {
     num_trials <- dim(x$results_each_trial)[1]
     cat(sprintf("\nPopulation d-value: %1.2f\n",x$pop_d))
     cat(sprintf("\nOriginal cell sizes: %d %d\nReplication cell sizes: %d %d",x$original_N1,x$original_N2,x$replication_N1,x$replication_N2))

     percent_level <- round(x$prob_level*100)
     cat(sprintf("\n\n%d%% Prediction interval capture percentage: %2.1f%% (%d of %d trials)",percent_level,x$percent_in_ri,x$in_prediction_interval_count,num_trials))
     cat(sprintf("\n%d%% Confidence interval capture percentage: %2.1f%% (%d of %d trials)",percent_level,x$percent_in_ci,x$in_confidence_interval_count,num_trials))

     table_out <- x$results_each_trial
     cat("\n\nIllustrative Trials:\n\n")
     print(table_out[1:5,],row.names=FALSE,digits = 2)
     cat("\n")


     cat("Note: n1 = original cell 1 size, n2 = original cell 2 size, d = original d-value,")
     cat("\nci.LL = lower-limit confidence interval, ci.UL = upper-limit confidence interval,")
     cat("\nrep.n1 = replication cell 1 size, rep.n2 = replication cell 2 size,")
     cat("\npi.LL = lower-limit prediction interval, pi.UL = upper-limit prediction interval,")
     cat("\nrep.d = replication d-value.\n")

}

get_orig_rep_d <- function(original_N1,original_N2,replication_N1,replication_N2,pop_d,prob_level,bias_correction) {

     #original sample
     original_cell1 <- rnorm(original_N1) + pop_d
     original_cell2 <- rnorm(original_N2)
     original_d <- get_d_value(original_cell1,original_cell2,bias_correction)

     #confidence and prediction interval
     prediction_interval_metrics <- pi.d(d=original_d,n1=original_N1,n2=original_N2,rep.n1=replication_N1,rep.n2=replication_N2,prob.level=prob_level)
     confidence_interval  <- c(prediction_interval_metrics$lower_confidence_interval,prediction_interval_metrics$upper_confidence_interval)
     prediction_interval <- c(prediction_interval_metrics$lower_prediction_interval,prediction_interval_metrics$upper_prediction_interval)

     #replication sample
     replication_cell1 <- rnorm(replication_N1) + pop_d
     replication_cell2 <- rnorm(replication_N2)
     replication_d <- get_d_value(replication_cell1,replication_cell2,bias_correction)

     #check if replication is in interval
     cur_result_in_prediction_interval <- FALSE
     is.in.ci <- is_value_in_interval(replication_d, confidence_interval)
     is.in.pi <- is_value_in_interval(replication_d, prediction_interval)

     ci_LL <- confidence_interval[1]
     ci_UL <- confidence_interval[2]
     ri_LL <- prediction_interval[1]
     ri_UL <- prediction_interval[2]

     output <- c(original_N1,original_N2,original_d,ci_LL,ci_UL,replication_N1,replication_N2,replication_d,ri_LL,ri_UL,is.in.ci,is.in.pi)
     output <- t(output)
     return(output)
}

