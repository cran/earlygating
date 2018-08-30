oc_rct <- function(N_c, N_e, delta, delta_power, confidence, e_a=0.5, e_b=0.5, c_a=0.5,
                   c_b=0.5, h_a=0.5, h_b=0.5, RR_h=NULL, N_h=NULL, w=NULL,
                   trues=seq(0,1,0.01), plot=T, legend = T, legend.pos="topleft"){

  oc_rct_nodyn <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b,
                           trues = seq(0,1,0.01), plot = T, legend = T, legend.pos="topleft"){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, plot=plot)[,2]
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

      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      lines(trues, ocs[,2], type = "l")
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
      }

    return(ocs)
  }


  oc_rct_dyn <- function(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h,
                         N_h, h_a, h_b, w, trues = seq(0,1,0.01), plot = T, legend = T, legend.pos="topleft"){

    ocs <- matrix(nrow = length(trues), ncol = 2)
    resp_cl <- req_resp_rct(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, h_a, h_b, RR_h, N_h, w, plot)[,2]
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
      plot(trues, ocs[,1], xlab = "True control response rate", ylab = "Probability", type = "l", ylim = c(0,max(ocs)), lty = 2,
           main = "Decision power and Decision type 1 error \n with respect to true control response rate")
      lines(trues, ocs[,2], type = "l")
      if(legend){
      legend(legend.pos, legend = c("Decision type 1 error", "Decision power"), bty = "n", lty = c(2,1))
      }
    }

    return(ocs)
  }

  if(!is.null(RR_h)&!is.null(N_h)&!is.null(w)){
    ret <- oc_rct_dyn(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, RR_h,
                      N_h, h_a, h_b, w, trues, plot, legend, legend.pos)
  }else{
    ret <- oc_rct_nodyn(N_c, N_e, delta, delta_power, confidence, e_a, e_b, c_a, c_b, trues, plot, legend, legend.pos)
  }

  res <- cbind(ret, trues)
  colnames(res) <- c("Dec. Alpha", "Dec. Power", "True control RR")
  round(res,4)

}
