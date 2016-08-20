# Using Modified Asymetrical Method (MAS) 
# Zou (2007), Equation 15
mas_interval <- function (original_r,l1,u1,l2,u2) {
     r1 <- original_r
     r2 <- r1
     LL <- r1 - sqrt((r1-l1)^2 + (u2-r2)^2)
     UL <- r1 + sqrt((u1-r1)^2 + (r2-l2)^2)
     
     replication_interval_out <- list()
     replication_interval_out$LL <- LL
     replication_interval_out$UL <- UL
     return(replication_interval_out)
}

# 
# quantile.r <- function (probability.in.tail,ncr,n,is.lower.tail) {
#      ncz <- atanh(ncr)
#      ncz_se <- 1 / sqrt(n-3)
#      z_tail <- stats::qnorm(p=probability.in.tail, mean=ncz, sd=ncz_se,lower.tail=is.lower.tail)
#      r_tail <- tanh(z_tail)
#      return(r_tail)
# }

# ci_m <- function(M,SD=NA,VAR=NA,n,conf_level) {
#      if (is.na(VAR)) {
#           VAR = SD * SD
#      }
#      one_tail_prob <- 1-((1-conf_level)/2) #turn 95 to 97.5 
#      mean_SE <- sqrt(VAR/n)
#      original_df <- n-1
#      l <- M - qt(one_tail_prob,original_df) * mean_SE
#      u <- M + qt(one_tail_prob,original_df) * mean_SE
#      
#      conf_interval_output <- list()
#      conf_interval_output$M <- M
#      conf_interval_output$SD <- SD
#      conf_interval_output$N <- n
#      conf_interval_output$lower_conf_limit_M <- l
#      conf_interval_output$upper_conf_limit_M <- u
#      return(conf_interval_output)
# }
# 

ci_r <- function(r,n,conf.level) {
     probability_in_tail <- (1- conf.level)/2
     
     obs_z <- atanh(r)
     obs_z_se <- 1/sqrt(n-3)
     
     high_z <- stats::qnorm(p=probability_in_tail, mean=obs_z,sd=obs_z_se,lower.tail=FALSE)
     high_r <- tanh(high_z)
     
     low_z <- stats::qnorm(p=probability_in_tail, mean=obs_z,sd=obs_z_se,lower.tail=TRUE)
     low_r <- tanh(low_z)
     
     
     conf_interval_output <- list()

     conf_interval_output$r <- r
     conf_interval_output$lower_conf_limit_r <- low_r
     conf_interval_output$upper_conf_limit_r <- high_r
     return(conf_interval_output)
}


# get_cor_data <- function(r, n) {
#      x <- scale(matrix(rnorm(n*2),n,2),center=T,scale=T)
#      svd_out <- svd(x)
#      u <- svd_out$u
#      
#      y=u[,1]
#      e=u[,2]
#      
#      if (round(r,3)!=0) {
#           r2 <- r*r
#           e_weight<-sqrt((1-r2)/r2)
#           x = y + e_weight*e
#           
#           vec_out=matrix(NA,nrow=n,ncol=2)
#           vec_out[,1]=scale(x,center=T,scale=T)
#           vec_out[,2]=scale(y,center=T,scale=T)
#      } else {
#           vec_out=matrix(NA,nrow=n,ncol=2)
#           vec_out[,1]=scale(y,center=T,scale=T)
#           vec_out[,2]=scale(e,center=T,scale=T)
#      }
#      df_out <- data.frame(vec_out)
#      names(df_out) <- c("x","y")
#      return(df_out)
# }

# get_d_value <- function(group1.data,group2.data) {
#      #confirm this is the unbiased one
#      d <- MBESS::smd(Group.1 = group1.data,Group.2 = group2.data,Unbiased = TRUE)
#      return(d)
# }


get_d_value <- function(group1.data,group2.data, bias_correction) {
     #confirm this is the unbiased one
     m1 <- mean(group1.data)
     m2 <- mean(group2.data)
     v1 <- var(group1.data)
     v2 <- var(group2.data)
     n1 <- length(group1.data)
     n2 <- length(group2.data)
     sp <- sqrt((v1*(n1-1) + v2*(n2-1))/(n1+n2-2))
     d_out <- (m1-m2)/sp
     
     if (bias_correction==TRUE) {
          m <- length(group1.data) + length(group2.data) - 2
          cm <- 1 - 3/(4*m-1)
          d_out <- d_out * cm
     }     
     
     return(d_out)
}

get_method_text <- function(orig_stat,LL,UL,replication_N,percent_level,stat_string) {
     method_text0 <- sprintf("The original %s is %1.2f -",stat_string,orig_stat)
     method_text1 <- sprintf("with a prediction interval %d%% PI[%1.2f, %1.2f] based a replication sample size of N = %d.",percent_level,LL,UL,replication_N)
     method_text2 <- sprintf("If the replication %s differs from the original %s only due to sampling error, there is a %d%% chance the replication result will fall in this interval.",stat_string,stat_string,percent_level)
     method_text3 <- sprintf("If the replication %s falls outside of this range, factors beyond sampling error are likely also responsible for the difference.",stat_string)
     
     output <- list()
     output$txt0 <- method_text0
     output$txt1 <- method_text1
     output$txt2 <- method_text2
     output$txt3 <- method_text3
     output$txt_combined <- combine_method_text(output)
     output$txt_ri <- sprintf("The prediction interval is %d%% PI[%1.2f, %1.2f].",percent_level,LL,UL)
     
     return(output)
}


combine_method_text <- function(x) {
     xout <- paste(x$txt0,x$txt1,x$txt2,x$txt3)
     return(xout)
}

get_method_text_d <- function(orig_stat,LL,UL,replication_N1,replication_N2,percent_level,stat_string) {
     output <- get_method_text(orig_stat,LL,UL,replication_N1,percent_level,stat_string) 
     output$txt1 <-  method_text1 <- sprintf("with a replication interval %d%% RI[%1.2f, %1.2f] based a replication cell sizes N1 = %d and N2 = %d.",percent_level,LL,UL,replication_N1,replication_N2)
     output$txt_combined <- combine_method_text(output)
     return(output)
}
     
     
     

get_cor <- function(x,bias_correction=FALSE) {
     n <- dim(x)[1]
     r <- cor(x[,1],x[,2])
     if (bias_correction==TRUE) {
          rho_hat <- r*(1+(1-r^2)/(2*n))
          r <- rho_hat
     }
     return(r)
}

is_value_in_interval <- function(value,interval) {
     is_in_interval <- FALSE
     check_interval<-findInterval(value,interval,rightmost.closed = TRUE)
     if (check_interval==1) {
          is_in_interval <- TRUE
     }
     return(is_in_interval)
}

