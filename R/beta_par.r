
beta_par <- function(mu_cov, phi_cov=NULL, orr, data, newdata, link = NULL, weights = NULL, plot = T){ # names in cov and orr

#require(betareg)

  fmla <- paste(orr, " ~ ", paste(mu_cov, collapse= "+"))

  if(!is.null(phi_cov)){
    fmla <- paste(fmla, "|", paste(phi_cov, collapse= "+"))
  }
  fmla <- as.formula(fmla)

  if(is.null(link)){
    if(is.null(weights)){
      modx <- betareg::betareg(fmla, data = data)
    }else{
      modx <- betareg::betareg(fmla, data = data, weights = weights)

    }
  a <- sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
         function(x) logLik(update(modx, link = x)))
  link <- names(which.max(a))
  }

  if(is.null(weights)){
  mod1 <- betareg::betareg(fmla, data = data, link = link)
  }else{
    mod1 <- betareg::betareg(fmla, data = data, link = link, weights = weights)

  }

  alpha <- mod1$coefficients$precision *
                predict(mod1, newdata = newdata)
  beta <- mod1$coefficients$precision - alpha

  if(plot){
    x <- NULL
     curve(dbeta(x,alpha, beta), xlab = "Response Rate", ylab = "Density",
           main = paste("Probability density function of the Beta (", round(alpha,2),";", round(beta,2),") distribution", sep = ""))
    }

  ret <- c(alpha, beta)
  names(ret) <- c("Alpha", "Beta")
  ret

}

