
req_resp <- function(N_e, delta, confidence, e_a=0.5, e_b=0.5, h_a=0.5, h_b=0.5,
                     RR_h=NULL, N_h=NULL, hist_RR_c=NULL, adapt = 1){

  req_resp_wr_true <- function(N_e, hist_RR_c, delta, confidence, e_a=1, e_b=1, adapt=1){
    crit <- 0
    e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
    e_b_temp <- e_b + N_e
    a <- try(while(qbeta(1-confidence, e_a_temp, e_b_temp) <= hist_RR_c+delta){ # iteratively solving the inequation
      crit <- crit + 1
      e_a_temp <- e_a + adapt*crit
      e_b_temp <- e_b - adapt*crit + N_e
      if(crit > N_e){ # if impossible to achieve number of responders with given sample size, then "NA"
        crit <- NA
        break
      }
    }, silent = T) # calculate critical value of "successes"
    if(class(a) == "try-error"){crit <- NA}
    return(crit)
  }

  req_resp_wr_hists <- function(N_e, delta, confidence, e_a=1, e_b=1, h_a=1,
                                h_b=1, RR_h, N_h, adapt=1){

    h_a_new <- h_a + RR_h*N_h
    h_b_new <- h_b - RR_h*N_h + N_h # it is assumed that "RR_h*N_h" responders were observed in historical control group
    crit <- 0
    e_a_temp <- e_a # starting values for Beta, i.e. 0 responders
    e_b_temp <- e_b + N_e
    a <- try(while(integrate(function(y) { dbeta(y, e_a_temp, e_b_temp)*
        sapply(y, function(y) {
          pbeta(y-delta, h_a_new, h_b_new)
        })
    }, delta, 1)$value <= confidence){ # iteratively solving the inequation
      crit <- crit + 1
      e_a_temp <- e_a + adapt*crit
      e_b_temp <- e_b - adapt*crit + N_e
      if(crit>N_e){# if impossible to achieve number of responders with given sample size, then "NA"
        crit <- NA
        break
      }
    }, silent = T)

    if(class(a) == "try-error"){crit <- NA}
    return(crit)
  }

  if(!is.null(RR_h) & !is.null(N_h)){
    ret <- req_resp_wr_hists(N_e, delta, confidence, e_a, e_b, h_a,
                      h_b, RR_h, N_h, adapt)
  }else{
    ret <- req_resp_wr_true(N_e, hist_RR_c, delta, confidence, e_a, e_b, adapt)
  }
  ret
}
