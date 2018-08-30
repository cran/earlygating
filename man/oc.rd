\name{oc}
\alias{oc}
\title{Single Arm Operating Characteristics}
\description{
Function for calculating the operating characteristics of the single arm bayesian designs in setting 1 and 2 for early gating.
}
\usage{
oc(N_e, delta, delta_power, confidence, e_a=0.5,
   e_b=0.5, h_a=0.5, h_b=0.5, RR_h=NULL, N_h=NULL,
   hist_RR_c=NULL, trues = seq(0,1,0.01),
   adapt = 1, plot = T, legend = T, legend.pos="topleft")
}

\arguments{
\item{N_e}{
Sample Size in the experimental group. Can be either a single value or a vector.
}
\item{delta}{
Required superiority to make a "GO" decision. Corresponds to \eqn{\delta}.
}
\item{delta_power}{
Superiority, at which decision power will be evaluated. Corresponds to \eqn{\bar{\delta}}{\delta bar}.
}
\item{confidence}{
Required confidence to make "GO" decision. Corresponds to \eqn{\gamma}{\gamma}.
}
\item{e_a}{
Alpha parameter of Beta Prior Distribution for the experimental response rate. Corresponds to \eqn{\alpha_e}{\alpha e}. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{e_b}{
Beta parameter of Beta Prior Distribution for the experimental response rate. Corresponds to \eqn{\beta_e}{\beta e}. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{h_a}{
Alpha parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\alpha_h}{\alpha h}. Only needs to be specified, if RR_h and N_h are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{h_b}{
Beta parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\beta_h}{\beta h}. Only needs to be specified, if RR_h and N_h are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{RR_h}{
Historical control response rate. Corresponds to \eqn{p_h}{p h}. If specified together with N_h, function will use setting 2 from pdf.
}
\item{N_h}{
Historical control sample size. Corresponds to \eqn{n_h}{n h}. If specified together with RR_h, function will use setting 2 from pdf.
}
\item{hist_RR_c}{
Point estimate of historical control repsonse rate. Corresponds to \eqn{\hat{p_h}}{p h hat}. If specified, while RR_h and N_h are not specified, function will use setting 1 from pdf.
}
\item{trues}{
Sequence of true control response rates and experimental response rates, at which the Probability to Go will be computed. Default is seq(0,1,0.01) to ensure continuous plots and accurate results.
}
\item{adapt}{
Level of adapting of experimental control rate to account for patient selection bias from phase II to phase III. Corresponds to \eqn{\xi}{\xi}. Default is 1, so no adapting.
}
\item{plot}{
Plots yes or no. Default is TRUE.
}
\item{legend}{
Logical; whether or not to include legend in plot. Default is TRUE.
}
\item{legend.pos}{
Position of legend. Default is "topleft".
}
}
\value{
A matrix containing the decision power and decision alpha with respect to the true control response rate.
}

\author{
Elias Laurin Meyer
}

\seealso{
{\code{\link{avg_oc_wr_ne}}},
{\code{\link{req_resp}}},
{\code{\link{oc_rct}}}
}
\examples{

# Setting 1
oc(N_e=50, delta=0.08, delta_power=0.13,
   confidence=0.6, hist_RR_c=0.5)

# Setting 2
oc(N_e=50, delta=0.08, delta_power=0.13,
  confidence=0.6, RR_h=0.5, N_h=50)
}

