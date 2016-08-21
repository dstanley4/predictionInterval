
# Interval Calculation ----------------------------------------------------

#' Correlation prediction interval
#' @param r Original study: Correlation
#' @param n Original study: Sample size
#' @param rep.n (optional) Replication study: Sample size. If not specified, n is used.
#' @param prob.level (optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @return The prediction interval and related statistics in list format.
#' @examples
#' pi.r(r=.35,n=100,rep.n=200)
#' @export
pi.r <- function (r,n,rep.n=NA,prob.level = .95) {
     original_r <- r
     original_N <- n
     prob_level <- prob.level
     if (is.na(rep.n)) {rep.n<-original_N}
     replication_N <- rep.n

     original_CI <- ci_r(original_r, original_N,conf.level = prob_level)
     l1 <- original_CI$lower_conf_limit_r
     u1 <- original_CI$upper_conf_limit_r


     mod_CI <- ci_r(original_r, replication_N,conf.level = prob_level)
     l2 <- mod_CI$lower_conf_limit_r
     u2 <- mod_CI$upper_conf_limit_r

     prediction_interval_values <- mas_interval(original_r=r,l1=l1,u1=u1,l2=l2,u2=u2)


     LL <- prediction_interval_values$LL
     UL <- prediction_interval_values$UL
     prediction_interval_metrics <- list()
     prediction_interval_metrics$original_r <- original_r
     prediction_interval_metrics$original_N <- original_N
     prediction_interval_metrics$replication_N <- replication_N
     prediction_interval_metrics$lower_prediction_interval <- LL
     prediction_interval_metrics$upper_prediction_interval <- UL
     prediction_interval_metrics$lower_confidence_interval <- l1
     prediction_interval_metrics$upper_confidence_interval <- u1
     prediction_interval_metrics$prob_level <- prob_level

     percent_level <- as.integer(round(prob_level*100))
     method_text <- get_method_text(original_r,LL,UL,replication_N,percent_level,"correlation")
     prediction_interval_metrics$method_text <- method_text$txt_combined
     prediction_interval_metrics$ri_text <- method_text$txt_ri

     class(prediction_interval_metrics) <- "r_prediction_interval"

     return(prediction_interval_metrics)
}




#' @export
print.r_prediction_interval <- function(x,...) {
     conf_per <- x$prob_level * 100

     cat(sprintf("\nOriginal study: r = %1.2f, N = %d, %d%% CI[%1.2f, %1.2f]",x$original_r,x$original_N,conf_per,x$lower_confidence_interval,x$upper_confidence_interval))
     cat(sprintf("\nReplication study: N = %d",x$replication_N))
     cat(sprintf("\nPrediction interval: %d%% PI[%1.2f,%1.2f].\n\n",conf_per,x$lower_prediction_interval,x$upper_prediction_interval))
     cat("\nInterpretation:\n")
     cat(x$method_text)

}




# Simulation --------------------------------------------------------------


#' Simulation to demonstrate the meaning of the correlation prediction interval
#' @param n Original study: Sample size
#' @param rep.n (optional) Replication study: Sample size. If not specified, n is used.
#' @param rho All samples are drawn from a common population. This specifies the population correlation.
#' @param number.trials Indicate the number of pairs of sample (original, replication) that should be used. 10,000 or higher suggested for stable results.
#' @param prob.level (optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @param bias.correction Apply bias correction formula to d-values.
#' @return The prediction interval capture percentage and related statistics in list format.
#' @examples
#' pi.r.demo(n=100,rho=.50,number.trials=10)
#' @export
pi.r.demo <- function(n=100,rep.n=NA,rho=.50,number.trials=10000,prob.level=.95, bias.correction=FALSE) {
     bias_correction <- bias.correction
     number_trials <- number.trials
     original_N <- n
     prob_level <- prob.level
     if (missing(rep.n)) {rep.n <- n}
     replication_N <- rep.n

     output <-pbapply::pbreplicate(number_trials,get_orig_rep_r(original_N=original_N,replication_N=replication_N,rho=rho,prob_level=prob_level,bias_correction=bias_correction))


     n <-output[,1,]
     r <-output[,2,]
     ci.LL <-output[,3,]
     ci.UL <-output[,4,]
     rep.n <-output[,5,]
     pi.LL <-output[,6,]
     pi.UL <-output[,7,]
     rep.r <-output[,8,]
     rep.r.in.ci <- as.logical(output[,9,])
     rep.r.in.pi <- as.logical(output[,10,])
     output_df <- data.frame(n,r,ci.LL,ci.UL,rep.n,pi.LL,pi.UL,rep.r,rep.r.in.ci,rep.r.in.pi)


     in_prediction_interval_count <- sum(rep.r.in.pi)
     in_confidence_interval_count <- sum(rep.r.in.ci)
     percent_in_pi <- (in_prediction_interval_count/(number_trials))*100
     percent_in_ci <- (in_confidence_interval_count/(number_trials))*100
     replication_demo_output <- list()
     replication_demo_output$percent_in_pi <- percent_in_pi
     replication_demo_output$percent_in_ci <- percent_in_ci
     replication_demo_output$in_prediction_interval_count<-in_prediction_interval_count
     replication_demo_output$in_confidence_interval_count<-in_confidence_interval_count
     replication_demo_output$results_each_trial <- output_df
     replication_demo_output$rho <- rho
     replication_demo_output$original_N <- original_N
     replication_demo_output$replication_N <- replication_N
     replication_demo_output$prob_level <- prob_level

     class(replication_demo_output) <- "replication_demo_r"

     return(replication_demo_output)
}


#' @export
print.replication_demo_r <- function(x,...) {
     num_trials <- dim(x$results_each_trial)[1]
     cat(sprintf("\nPopulation correlation: %1.2f\n",x$rho))
     cat(sprintf("\nOriginal sample size: %d\nReplication sample size: %d",x$original_N,x$replication_N))

     percent_level <- round(x$prob_level*100)
     cat(sprintf("\n\n%d%% Prediction interval capture percentage: %2.1f%% (%d of %d trials)",percent_level,x$percent_in_pi,x$in_prediction_interval_count,num_trials))
     cat(sprintf("\n%d%% Confidence interval capture percentage: %2.1f%% (%d of %d trials)",percent_level,x$percent_in_ci,x$in_confidence_interval_count,num_trials))

     table_out <- x$results_each_trial[1:5,]
     table_out$r <- round(table_out$r,2)
     table_out$ci.LL <- round(table_out$ci.LL,2)
     table_out$ci.UL <- round(table_out$ci.LL,2)
     table_out$pi.LL <- round(table_out$pi.LL,2)
     table_out$pi.UL <- round(table_out$pi.UL,2)
     table_out$rep.r <- round(table_out$rep.r,2)

     cat("\n\nIllustrative Trials:\n\n")
     print(table_out,row.names=FALSE,digits = 2)
     cat("\n")

     cat("Note: n = original sample size, r = original correlation,")
     cat("\nci.LL = lower-limit confidence interval, ci.UL = upper-limit confidence interval, rep.n = replication sample size,")
     cat("\npi.UL = lower-limit prediction interval, pi.UL = upper-limit prediction interval,")
     cat("\nrep.r = replication correlation.\n")
}




get_orig_rep_r <- function(original_N,replication_N,rho,prob_level,bias_correction) {
     sigma <- diag(2)
     sigma[2,1] <- rho
     sigma[1,2] <- rho

     #original sample
     original_sample <- MASS::mvrnorm(n=original_N,mu=rep(0,2),Sigma=sigma,empirical=FALSE)
     original_r <- get_cor(original_sample,bias_correction=bias_correction)

     #confidence and prediction interval
     prediction_interval_metrics <- pi.r(r=original_r,n=original_N,rep.n=replication_N,prob.level=prob_level)
     confidence_interval  <- c(prediction_interval_metrics$lower_confidence_interval,prediction_interval_metrics$upper_confidence_interval)
     prediction_interval <- c(prediction_interval_metrics$lower_prediction_interval,prediction_interval_metrics$upper_prediction_interval)

     #replication sample
     replication_sample <- MASS::mvrnorm(n=replication_N,mu=rep(0,2),Sigma=sigma,empirical=FALSE)
     replication_r <- get_cor(replication_sample,bias_correction=bias_correction)

     #check if replication is in interval
     is.in.ci <- is_value_in_interval(replication_r, confidence_interval)
     is.in.pi <- is_value_in_interval(replication_r, prediction_interval)

     ci_LL <- confidence_interval[1]
     ci_UL <- confidence_interval[2]
     pi_LL <- prediction_interval[1]
     pi_UL <- prediction_interval[2]
     output <- c(original_N,original_r,ci_LL,ci_UL,replication_N,pi_LL,pi_UL,replication_r,is.in.ci,is.in.pi)
     output <- t(output)
     return(output)
}


