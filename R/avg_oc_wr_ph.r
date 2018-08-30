
avg_oc_wr_ph <- function(N_e, delta, delta_power, confidence, e_a=0.5, e_b=0.5, alpha_c,
                         beta_c, trues=seq(0,1,0.01), adapt=1, plot = T, legend = T, legend.pos="topleft"){
  average_oc <- function(true_RR_e, N_e, delta, delta_power, confidence,
                         e_a, e_b, hist_RR_c, adapt, alpha_c, beta_c, trues, plot = T, legend = T, legend.pos){

    true_RR_e <- trues # Probability to GO will be evaluated along the same sequence as power and type 1 error
    true_RR_c <- hist_RR_c

    req_resp_wr_true <- function(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt){
      crit <- 0
      e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
      e_b_temp <- e_b + N_e
      while(qbeta(1-confidence, e_a_temp, e_b_temp) <= hist_RR_c+delta){ # iteratively solving the inequation
        crit <- crit + 1
        e_a_temp <- e_a + adapt*crit
        e_b_temp <- e_b - adapt*crit + N_e
        if(crit > N_e){ # if impossible to achieve number of responders with given sample size, then "NA"
          crit <- NA
          break
        }
      } # calculate critical value of "successes"
      return(crit)
    }

    # Function for Setting 1 in order to calculate power and type 1 error for a given setting and a fixed \theta_c
    # ret == should there be a return value
    # other parameters as above
    oc <- function(true_RR_e, N_e, true_RR_c, hist_RR_c, delta, confidence, e_a, e_b, delta_power,
                   adapt, plot = F, ret = T){

      # firstly calculate how many "successes" are needed for "GO" Decision
      crit <- req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt)

      if(plot){
        P_go <- pbinom(crit-1, N_e, true_RR_e, lower.tail = FALSE)
      }# depending on how many successes are needed for GO, what is probability
      #   with respect to different experimental response rates
      if(true_RR_c+delta_power<=1){ # power at true_RR_c + delta power, unless true_RR_c + delta power > 1, in which case power = 1
        power_ach <- pbinom(crit-1, N_e, true_RR_c+delta_power, lower.tail = FALSE)
      }else{power_ach <- 0}
      alpha_ach <- pbinom(crit-1, N_e, true_RR_c, lower.tail = FALSE) # alpha at true_RR_c
      if(plot) {
        plot(true_RR_e, P_go, type = "l", xlab = "True experimental response rate",
             ylab = "Prob(GO)", main = "Prob(GO) w/r to True experimental response rate")
        abline(v = true_RR_e[which(round(true_RR_e,5) == round(true_RR_c+delta_power, 5))],
               h = P_go[which(round(true_RR_e,5) == round(true_RR_c+delta_power, 5))])
        abline(v = true_RR_e[which(round(true_RR_e,5) == round(true_RR_c, 5))],
               h = P_go[which(round(true_RR_e,5) == round(true_RR_c, 5))])

      }
      if(ret){return(c(alpha_ach, power_ach, crit))}
      #return achieved power, alpha and critical value
    }

      n_resp_wr_trues <- numeric(length(trues))
      if(plot){
        for(k in 1:length(trues)){
          n_resp_wr_trues[k] <- req_resp_wr_true(N_e, trues[k], delta, confidence, e_a, e_b, adapt)
        } # firstly compute minimum required number of successes with respect to control RR
        if(plot){
          plot(trues, n_resp_wr_trues, ylab = "N Successes Req for GO", xlab = "Historical Control response rate", type = "l", main = "Required Number of successes with \n respect to the historical control response rate")
          abline(h = req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt),
                 v = hist_RR_c)
        }
      }

      if(plot){
        # Now draw only the plots (== ret=F)
        oc(true_RR_e, N_e, true_RR_c, hist_RR_c, delta, confidence, e_a, e_b, delta_power, adapt, plot = plot, ret = F)
      }

      # Now compute OCs with respect to true_control_RR (i.e. \theta_c is now varying)
      oc_wr_trues <- matrix(nrow = length(trues), ncol = 2)
      crit <- req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt)
      for(j in 1:length(trues)){
        if(trues[j]+delta_power<=1){# power at true_RR_c + delta power, unless true_RR_c + delta power > 1, in which case power = 1
          oc_wr_trues[j,2] <- pbinom(crit-1, N_e, trues[j]+delta_power, lower.tail = FALSE)
        }else{oc_wr_trues[j,2] <- 0}
        oc_wr_trues[j,1] <- pbinom(crit-1, N_e, trues[j], lower.tail = FALSE) # alpha at true_RR_c
      }

      if(plot){
        plot(trues,oc_wr_trues[,1], type = "l", lty = 2, ylab = "Probability", xlab = "True control response rate", main = "Power and Type 1 error \n w/r to True Control response rate")
        lines(trues,oc_wr_trues[,2], type = "l", lty = 1)
        if(legend){
        legend(legend.pos, legend = c("Type 1 error", "Power"), bty = "n", lty = c(2,1))
        }
      }

      # calculate average OCs as mean of OC*pdf of Beta for true control RR
      av_power <- mean(oc_wr_trues[,2]*dbeta(trues,alpha_c,beta_c))
      av_alpha <- mean(oc_wr_trues[,1]*dbeta(trues,alpha_c,beta_c))

      if(plot){
        plot(trues,oc_wr_trues[,1], type = "l", lty = 2, ylab = "Probability", xlab = "True control response rate",
             main = paste("Average Alpha, Power = \n ", round(av_alpha,2), ",", round(av_power,2)))
        lines(trues,oc_wr_trues[,2], type = "l", lty = 1)
        abline(v = hist_RR_c)
        x <- NULL
        curve(1/max(alpha_c, beta_c)*dbeta(x,alpha_c,beta_c), from = 0, to = 1, n=101, add = T, lty = 3)
        if(legend){
        legend(legend.pos, legend = c("Type 1 error", "Power", "Beta density"), bty = "n", lty = c(2,1,3))
        }
      }

    return(c(av_alpha, av_power))

  }


  ph <- trues
  true_RR_e <- trues
  av_oc_wr_ph <- matrix(nrow = length(ph), ncol = 2)
  for(i in 1:length(ph)){

  av_oc_wr_ph[i,] <- average_oc(true_RR_e=true_RR_e, N_e=N_e, delta=delta,
                                  delta_power=delta_power, confidence=confidence, e_a=e_a, e_b=e_b,
                                  hist_RR_c=ph[i], adapt=adapt,
                                  alpha_c=alpha_c, beta_c=beta_c, trues=trues, plot = F)
  }
  if(plot){
    plot(trues, av_oc_wr_ph[,1], xlab = "Target ph", ylab = "Probability", type = "l", lty = 2,
         main = "Average OCs w/r to target")
    lines(trues, av_oc_wr_ph[,2])
    if(legend){
    legend(legend.pos, legend = c("Decision power", "Decision type 1 Error"), lty = c(1,2))
    }
  }
  res <- cbind(av_oc_wr_ph, trues)
  colnames(res) <- c("Dec.Alpha", "Dec.Power", "p_h")
  res
}
