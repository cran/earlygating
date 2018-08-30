
avg_oc_wr_ne_rct <- function(N_c, N_e, delta, delta_power, confidence, e_a=0.5, e_b=0.5, c_a=0.5, c_b=0.5,
                             h_a=0.5, h_b=0.5, N_h = NULL, RR_h = NULL, w=NULL, alpha_c, beta_c,
                             trues = seq(0,1,0.01), plot = T, coresnum = NULL, legend = T,
                             legend.pos = "topleft"){

  req_resp_rct <- function(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b){
    resp <- rep(NA, N_c+1)
    for(i in 0:N_c){

      suc_c <- i
      c_a_new <- c_a + suc_c
      c_b_new <- c_b + N_c - suc_c

      crit <- 0
      e_a_temp <- e_a
      e_b_temp <- e_b + N_e

      while(integrate(function(y) { dbeta(y, e_a_temp, e_b_temp)*
          sapply(y, function(y) {
            pbeta(y-delta, c_a_new, c_b_new)
          })
      }, delta, 1)$value <= confidence){
        e_a_temp <- e_a + crit
        e_b_temp <- e_b - crit + N_e
        crit <- crit + 1
        if(crit>N_e){
          crit <- N_e + 1
          break
        }
      }

      resp[i+1] <- crit
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }

  oc_rct <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, trues = seq(0,1,0.01), plot = T, legend=T, legend.pos){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b)
    resp <- resp_cl
    resp[is.na(resp)] <- N_e+1
    for(j in 1:length(trues)){
      alphas <- pbinom(resp-1, N_e, trues[j], lower.tail = FALSE)
      weighted_alphas <- alphas * dbinom(0:N_c, N_c, trues[j])

      powers <- suppressWarnings(pbinom(resp-1, N_e, trues[j]+delta_power, lower.tail = FALSE))
      powers[is.nan(powers)] <- 0
      weighted_powers <- powers * dbinom(0:N_c, N_c, trues[j])

      alpha <- sum(weighted_alphas)
      power <- sum(weighted_powers)

      ocs[j,1] <- alpha
      ocs[j,2] <- power
    }

    if(plot){

      plot(0:N_c,resp_cl, xlab = "Successes in Control Group", ylab = "Successes in Experimental Group",
           main = "Minimum required number of successes in Exp Group to make \n  GO decision w/r to Successes in Con Group", type = "s")

      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      lines(trues, ocs[,2], type = "l")
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
    }

    return(ocs)
  }


  req_resp_rct_dyn <- function(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w){
    resp <- rep(NA, N_c+1)
    for(i in 0:N_c){

      suc_c <- i

      w1 <- w*beta(suc_c + N_h*RR_h + h_a, N_h + N_c - suc_c - N_h*RR_h + h_b)/beta(N_h*RR_h + h_a, N_h - N_h*RR_h + h_b)
      w2 <- (1-w)*beta(suc_c + c_a, N_c - suc_c + c_b)/beta(c_a, c_b)
      w1n <- w1/(w1+w2)
      w2n <- w2/(w1+w2)

      m_a <- w1n*(suc_c + N_h*RR_h + h_a) + w2n*(suc_c + c_a)
      m_b <- w1n*(N_h + N_c - suc_c - N_h*RR_h + h_b) + w2n*(N_c - suc_c + c_b)

      crit <- 0
      e_a_temp <- e_a
      e_b_temp <- e_b + N_e

      while(integrate(function(y) { dbeta(y, e_a_temp, e_b_temp)*
          sapply(y, function(y) {
            pbeta(y-delta, m_a, m_b)
          })
      }, delta, 1)$value <= confidence){
        e_a_temp <- e_a + crit
        e_b_temp <- e_b - crit + N_e
        crit <- crit + 1
        if(crit>N_e){
          crit <- N_e + 1
          break
        }
      }

      resp[i+1] <- crit
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }


  oc_rct_dyn <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w, trues = seq(0,1,0.01), plot = T, legend = T, legend.pos){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct_dyn(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w)
    resp <- resp_cl
    resp[is.na(resp)] <- N_e+1
    for(j in 1:length(trues)){
      alphas <- pbinom(resp-1, N_e, trues[j], lower.tail = FALSE)
      weighted_alphas <- alphas * dbinom(0:N_c, N_c, trues[j])

      powers <- suppressWarnings(pbinom(resp-1, N_e, trues[j]+delta_power, lower.tail = FALSE))
      powers[is.nan(powers)] <- 0
      weighted_powers <- powers * dbinom(0:N_c, N_c, trues[j])

      alpha <- sum(weighted_alphas)
      power <- sum(weighted_powers)

      ocs[j,1] <- alpha
      ocs[j,2] <- power
    }

    if(plot){
      plot(0:N_c,resp_cl, xlab = "Successes in Control Group", ylab = "Successes in Experimental Group",
           main = "Minimum required number of successes in Exp Group to make \n  GO decision w/r to Successes in Con Group", type = "s")

      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      lines(trues, ocs[,2], type = "l")
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
    }

    return(ocs)
  }


  if(!is.null(RR_h) & !is.null(N_h) & !is.null(h_a) & !is.null(h_b) & !is.null(w)){


    if(length(N_e)>1){

      av_power <- numeric(length(N_e))
      av_alpha <- numeric(length(N_e))

      if(is.null(coresnum)){coresnum <- parallel::detectCores()-1}
      cores <- coresnum
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)

      "%dopar%" <- foreach::"%dopar%"
      i <- NULL
      res <- foreach::foreach(i = 1:length(N_e), .combine = 'rbind') %dopar%{

        ocs <- oc_rct_dyn(N_c[i], N_e[i], delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w, plot = F)

        av_power[i] <- mean(ocs[,2]*dbeta(trues,alpha_c,beta_c))
        av_alpha[i] <- mean(ocs[,1]*dbeta(trues,alpha_c,beta_c))

        return(c(av_alpha[i], av_power[i], N_c[i], N_e[i]))

      }

      doParallel::stopImplicitCluster()
      closeAllConnections()

      if(plot){
      par(mar = c(5,6,7,5))
      plot(N_e, res[,1], type = "l", xlab = "Sample size in experimental group", ylab = "Probability", lty = 2, ylim = c(0,1))
      title("Average decision power and decision type 1 error with respect to \n sample size and uncertainty about the true control response rate", line = 4.7)
      lines(N_e, res[,2])
      abline(h = 0.10, lty = 3)
      abline(h = 0.60, lty = 3)
      if(legend){
      legend(legend.pos, legend = c("Avg. Dec. Alpha", "Avg. Dec. Power"), bty = "n", lty = c(2,1))
      }
      par(new = T)
      plot(N_c, rep(5, length(N_c)), xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = NA_integer_)
      axis(side=3)
      mtext("Sample Size in the control group", side = 3, line = 2.5)
      }

      colnames(res) <- c("Avg. Dec. Alpha", "Avg. Dec. Power", "Ne", "Nc")
      return(round(res,2))

    }else{

      res <- oc_rct_dyn(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w, plot = plot, legend.pos=legend.pos)
      av_power <- mean(res[,2]*dbeta(trues,alpha_c,beta_c))
      av_alpha <- mean(res[,1]*dbeta(trues,alpha_c,beta_c))

      if(plot){
      plot(trues, res[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(res)), lty = 2, main = paste("Average Decision alpha, Decision power = \n", round(av_alpha,2),",", round(av_power,2)))
      lines(trues, res[,2], type = "l")
      x <- NULL
      curve(1/(2*max(alpha_c, beta_c))*dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power", "Beta density"), bty = "n", lty = c(2,1,3))
      }
      }

      ret <- round(c(av_alpha, av_power),2)
      names(ret) <- c("Avg. Dec. Alpha", "Avg. Dec. Power")
      return(ret)


    }


  }else{

  if(length(N_e)>1){

  av_power <- numeric(length(N_e))
  av_alpha <- numeric(length(N_e))

  if(is.null(coresnum)){coresnum <- parallel::detectCores()-1}
  cores <- coresnum
  cl <- parallel::makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)

  "%dopar%" <- foreach::"%dopar%"
  res <- foreach::foreach(i = 1:length(N_e), .combine = 'rbind') %dopar%{

    ocs <- oc_rct(N_c[i], N_e[i], delta, delta_power, confidence, e_a, e_b, c_a, c_b, plot = F)

    av_power[i] <- mean(ocs[,2]*dbeta(trues,alpha_c,beta_c))
    av_alpha[i] <- mean(ocs[,1]*dbeta(trues,alpha_c,beta_c))

    return(c(av_alpha[i], av_power[i], N_c[i], N_e[i]))

  }

  doParallel::stopImplicitCluster()
  closeAllConnections()

  if(plot){
    par(mar = c(5,6,7,5))
    plot(N_e, res[,1], type = "l", xlab = "Sample size in experimental group", ylab = "Probability", lty = 2, ylim = c(0,1))
    title("Average decision power and decision type 1 error with respect to \n sample size and uncertainty about the true control response rate", line = 4.7)
    lines(N_e, res[,2])
    abline(h = 0.10, lty = 3)
    abline(h = 0.60, lty = 3)
    if(legend){
    legend(legend.pos, legend = c("Avg. Dec. Alpha", "Avg. Dec. Power"), bty = "n", lty = c(2,1))
    }
    par(new = T)
    plot(N_c, rep(5, length(N_c)), xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = NA_integer_)
    axis(side=3)
    mtext("Sample Size in the control group", side = 3, line = 2.5)
  }

  colnames(res) <- c("Avg. Dec. Alpha", "Avg. Dec. Power", "Ne", "Nc")
  return(round(res,2))

  }else{

    res <- oc_rct(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, plot = plot, legend.pos=legend.pos)
    av_power <- mean(res[,2]*dbeta(trues,alpha_c,beta_c))
    av_alpha <- mean(res[,1]*dbeta(trues,alpha_c,beta_c))

    if(plot){
    plot(trues, res[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(res)), lty = 2, main = paste("Average Decision alpha, Decision power = \n", round(av_alpha,2),",", round(av_power,2)))
    lines(trues, res[,2], type = "l")
    curve(1/(2*max(alpha_c, beta_c))*dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
    if(legend){
    legend(legend.pos, legend = c("Decision type 1 error", "Decision power", "Beta density"), bty = "n", lty = c(2,1,3))
    }

    }

    ret <- round(c(av_alpha, av_power),2)
    names(ret) <- c("Avg. Dec. Alpha", "Avg. Dec. Power")
    return(ret)

  }
  }

}
