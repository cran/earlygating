
req_resp_rct <- function(N_c, N_e, delta, confidence, e_a=0.5, e_b=0.5, c_a=0.5, c_b=0.5,
                         h_a=0.5, h_b=0.5, RR_h=NULL, N_h=NULL, w=NULL, plot=T){

  req_resp_rct_nodyn <- function(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b){
    resp <- rep(NA, N_c+1)
    for(i in 0:N_c){

      suc_c <- i
      c_a_new <- c_a + suc_c
      c_b_new <- c_b + N_c - suc_c

      crit <- 0
      e_a_temp <- e_a
      e_b_temp <- e_b + N_e

      a <- try(
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
        }, silent = T)


      if(class(a) == "try-error"){resp[i+1] <- NA}else{
        resp[i+1] <- crit
      }
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
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

      a <- try(
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
        }, silent = T)

      if(class(a) == "try-error"){resp[i+1] <- NA}else{
        resp[i+1] <- crit
      }
    }
    resp_cl <- resp
    resp_cl[resp_cl == N_e+1] <- NA

    return(resp_cl)
  }

  if(!is.null(RR_h)&!is.null(N_h)&!is.null(w)){
    ret <- req_resp_rct_dyn(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b, RR_h, N_h, h_a, h_b, w)
  }else{
    ret <- req_resp_rct_nodyn(N_c, N_e, delta, confidence, e_a, e_b, c_a, c_b)
  }
if(plot){
  plot(0:N_c,ret, xlab = "Successes in Control Group", ylab = "Successes in Experimental Group",
       main = "Minimum required number of successes in Exp Group to make \n  GO decision w/r to Successes in Con Group", type = "s")
  }
  res <- cbind(0:N_c,ret)
  colnames(res) <- c("Suc. in Con.", "Req. Suc. in Exp.")
  res
}
