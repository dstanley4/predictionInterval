# Interval Calculation ----------------------------------------------------

#' Prediction interval for the mean
#' @param M Original study: Mean
#' @param SD Original study: Standard deviation. Provide this or variance - not both.
#' @param VAR Original study: Variance. Provide this or standard deviation - not both.
#' @param n Original study: Sample size
#' @param rep.n (optional) Replication study: Sample size. If not specified, n is used.
#' @param prob.level (optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @return The prediction interval and related statistics in list format.
#' @examples
#' pi.m(M=2.53,SD=1.02,n=40,rep.n=80)
#' @export
pi.m <- function (M,SD=NA,VAR=NA,n,rep.n=NA,prob.level = .95) {
     original_M <- M
     original_N <- n
     prob_level <- prob.level
     if (is.na(rep.n)) {rep.n<-original_N}
     replication_N <- rep.n

     if (is.na(VAR)) {
          original_VAR <- SD*SD
     } else {
          original_VAR <- VAR
     }


     #Confidence Interval for Mean
     one_tail_prob <- 1-((1-prob_level)/2)
     mean_SE <- sqrt(original_VAR/original_N)
     original_df <- original_N-1
     ci_LL <- original_M - qt(one_tail_prob,original_df) * mean_SE
     ci_UL <- original_M + qt(one_tail_prob,original_df) * mean_SE

     #Prediction Interval for Mean
     difference_SE <- sqrt(original_VAR/original_N + original_VAR/replication_N)
     #replication_df <- original_N + replication_N - 2
     replication_df <- original_N - 1
     ri_LL <- original_M - qt(one_tail_prob,replication_df) * difference_SE
     ri_UL <- original_M + qt(one_tail_prob,replication_df) * difference_SE

     prediction_interval_metrics <- list()
     prediction_interval_metrics$original_M <- original_M
     prediction_interval_metrics$original_VAR <- original_VAR
     prediction_interval_metrics$original_SD <- sqrt(original_VAR)
     prediction_interval_metrics$original_N <- original_N
     prediction_interval_metrics$replication_N <- replication_N
     prediction_interval_metrics$lower_prediction_interval <- ri_LL
     prediction_interval_metrics$upper_prediction_interval <- ri_UL
     prediction_interval_metrics$lower_confidence_interval <- ci_LL
     prediction_interval_metrics$upper_confidence_interval <- ci_UL
     prediction_interval_metrics$prob_level <- prob_level

     percent_level <- as.integer(round(prob_level*100))
     method_text <- get_method_text(original_M,ri_LL,ri_UL,replication_N,percent_level,"mean")
     prediction_interval_metrics$method_text <- method_text$txt_combined
     prediction_interval_metrics$ri_text <- method_text$txt_ri

     class(prediction_interval_metrics) <- "M_prediction_interval"

     return(prediction_interval_metrics)
}

#' @export
print.M_prediction_interval <- function(x,...) {
     conf_per <- x$prob_level * 100
     cat(sprintf("\nOriginal study: M = %1.2f, SD = %1.2f, N = %d, %d%% CI[%1.2f, %1.2f]",x$original_M,x$original_SD,x$original_N,conf_per,x$lower_confidence_interval,x$upper_confidence_interval))
     cat(sprintf("\nReplication study: N = %d",x$replication_N))
     cat(sprintf("\nPrediction interval: %d%% PI[%1.2f,%1.2f].\n\n",conf_per,x$lower_prediction_interval,x$upper_prediction_interval))
     cat("\nInterpretation:\n")
     cat(x$method_text)

}




# Simulation --------------------------------------------------------------





#' Simulation to demonstrate the meaning of the prediction interval for the mean
#' @param n Original study: Sample size
#' @param rep.n (optional) Replication study: Sample size. If not specified, n is used.
#' @param mu All samples are drawn from a common population. This specifies the population correlation.
#' @param sigma All samples are drawn from a common population. This specifies the population standard deviation.
#' @param number.trials Indicate the number of pairs of sample (original, replication) that should be used. 10,000 or higher suggested for stable results.
#' @param prob.level (optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.
#' @param show.all.trials Show original correlation, prediction interval, replication correlation, and whether replication effect is in the interval.
#' @return The prediction interval capture percentage and related statistics in list format.
#' @examples
#' pi.m.demo(n=150,mu=0,sigma=1,number.trials=10)
#' @export
pi.m.demo <- function(n=10,rep.n=NA,mu=0,sigma=1,number.trials=10000,prob.level=.95,show.all.trials=FALSE) {
     number_trials <- number.trials
     original_N <- n
     prob_level <- prob.level
     if (missing(rep.n)) {rep.n <- n}
     replication_N <- rep.n

     output <-pbapply::pbreplicate(number_trials,get_orig_rep_m(original_N=original_N,replication_N=replication_N,mu=mu,sigma=sigma,prob_level=prob_level))

     #output <- c(original_N,original_M,original_VAR,ci_LL,ci_UL,replication_N,ri_LL,ri_UL,replication_M,is.in.ci,is.in.pi)

     n <- output[,1,]
     M <- output[,2,]
     SD <- sqrt(output[,3,])
     ci.LL <- output[,4,]
     ci.UL <- output[,5,]
     rep.n <- output[,6,]
     pi.LL <- output[,7,]
     pi.UL <- output[,8,]
     rep.M <-output[,9,]
     rep.M.in.ci <- as.logical(output[,10,])
     rep.M.in.pi <- as.logical(output[,11,])
     output_df <- data.frame(n,M,SD,ci.LL,ci.UL,rep.n,pi.LL,pi.UL,rep.M,rep.M.in.ci,rep.M.in.pi)


     in_prediction_interval_count <- sum(rep.M.in.pi)
     in_confidence_interval_count <- sum(rep.M.in.ci)
     percent_in_ri <- (in_prediction_interval_count/(number_trials))*100
     percent_in_ci <- (in_confidence_interval_count/(number_trials))*100
     replication_demo_output <- list()
     replication_demo_output$percent_in_ri <- percent_in_ri
     replication_demo_output$percent_in_ci <- percent_in_ci
     replication_demo_output$in_prediction_interval_count <- in_prediction_interval_count
     replication_demo_output$in_confidence_interval_count <- in_confidence_interval_count
     replication_demo_output$results_each_trial <- output_df
     replication_demo_output$mu <- mu
     replication_demo_output$sigma <- sigma
     replication_demo_output$original_N <- original_N
     replication_demo_output$replication_N <- replication_N
     replication_demo_output$prob_level <- prob_level

     class(replication_demo_output) <- "replication_demo_M"

     return(replication_demo_output)
}




#' @export
print.replication_demo_M <- function(x,...) {
     num_trials <- dim(x$results_each_trial)[1]
     cat(sprintf("\nPopulation mean: %1.2f\nPopulation standard deviation: %1.2f\n",x$mu, x$sigma))
     cat(sprintf("\nOriginal sample size: %d\nReplication sample size: %d",x$original_N,x$replication_N))
     percent_level <- round(x$prob_level*100)
     cat(sprintf("\n\n%d%% Prediction interval capture percentage: %2.1f%% (%d of %d trials)",percent_level,x$percent_in_ri,x$in_prediction_interval_count,num_trials))
     cat(sprintf("\n%d%% Confidence interval capture percentage: %2.1f%% (%d of %d trials)",percent_level,x$percent_in_ci,x$in_confidence_interval_count,num_trials))

     table_out <- x$results_each_trial[1:5,]
     cat("\n\nIllustrative Trials:\n\n")
     table_out$M <- round(table_out$M,2)
     table_out$SD <- round(table_out$SD,2)
     table_out$ci.LL <- round(table_out$ci.LL,2)
     table_out$ci.UL <- round(table_out$ci.LL,2)
     table_out$pi.LL <- round(table_out$pi.LL,2)
     table_out$pi.UL <- round(table_out$pi.UL,2)
     table_out$rep.M <- round(table_out$rep.M,2)
     print(table_out,row.names=FALSE)
     cat("\n")

     cat("\nNote: n = original sample size, M = original mean, SD = original standard deviation,")
     cat("\nci.LL = lower-limit confidence interval, ci.UL = upper-limit confidence interval, rep.n = replication sample size,")
     cat("\npi.LL = lower-limit prediction interval, pi.UL = upper-limit prediction interval, rep.M = replication mean\n")
}


get_orig_rep_m <- function(original_N,replication_N,mu,sigma,prob_level) {

     #original sample
     original_sample <- rnorm(n=original_N,mean=mu,sd=sigma)
     original_M <- mean(original_sample)
     original_VAR <- var(original_sample)


     #confidence and prediction intervals
     prediction_interval_metrics <- pi.m(M=original_M,VAR=original_VAR,n=original_N,rep.n=replication_N,prob.level=prob_level)
     confidence_interval  <- c(prediction_interval_metrics$lower_confidence_interval,prediction_interval_metrics$upper_confidence_interval)
     prediction_interval <- c(prediction_interval_metrics$lower_prediction_interval,prediction_interval_metrics$upper_prediction_interval)

     #replication sample
     replication_sample <- rnorm(n=replication_N,mean=mu,sd=sigma)
     replication_M <- mean(replication_sample)

     is.in.ci <- is_value_in_interval(replication_M, confidence_interval)
     is.in.pi <- is_value_in_interval(replication_M, prediction_interval)

     #check if replication is in interval
     ci_LL <- confidence_interval[1]
     ci_UL <- confidence_interval[2]
     ri_LL <- prediction_interval[1]
     ri_UL <- prediction_interval[2]
     output <- c(original_N,original_M,original_VAR,ci_LL,ci_UL,replication_N,ri_LL,ri_UL,replication_M,is.in.ci,is.in.pi)
     output <- t(output)
     return(output)
}


