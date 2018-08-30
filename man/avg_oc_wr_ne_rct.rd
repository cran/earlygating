\name{avg_oc_wr_ne_rct}
\alias{avg_oc_wr_ne_rct}
\title{RCT Average Operating Characteristics}
\description{
Function for calculating the average operating characteristics of two RCT bayesian designs for early gating with respect to the sample size in the experimental group, the sample size in the control group and possible historical data.
}
\usage{
avg_oc_wr_ne_rct(N_c, N_e, delta, delta_power, confidence,
                 e_a=0.5, e_b=0.5, c_a=0.5, c_b=0.5, h_a=0.5,
                 h_b=0.5, N_h = NULL, RR_h = NULL, w=NULL,
                 alpha_c, beta_c, trues = seq(0,1,0.01),
                 plot = T, coresnum = NULL, legend = T,
                 legend.pos = "topleft")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N_c}{
Sample Size in the control group. Can be either a single value or a vector, but needs to be the same length as N_e.
}
\item{N_e}{
Sample Size in the experimental group. Can be either a single value or a vector, but needs to be the same length as N_c.
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
\item{c_a}{
Alpha parameter of Beta Prior Distribution for the control response rate. Corresponds to \eqn{\alpha_c}{\alpha c}. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{c_b}{
Beta parameter of Beta Prior Distribution for the control response rate. Corresponds to \eqn{\beta_c}{\beta c}. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{h_a}{
Alpha parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\alpha_h}{\alpha h}. Only needs to be specified, if RR_h, N_h and w are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{h_b}{
Beta parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\beta_h}{\beta h}. Only needs to be specified, if RR_h, N_h and w are also specified. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{RR_h}{
Historical control response rate. Corresponds to \eqn{p_h}{p h}. If specified together with N_h and w, function will use setting 4 from pdf.
}
\item{N_h}{
Historical control sample size. Corresponds to \eqn{n_h}{n h}. If specified together with RR_h and w, function will use setting 4 from pdf.
}
\item{w}{
Level of dynmaic borrowing. Corresponds to \eqn{w}{w}.
}
\item{alpha_c}{
Alpha parameter of Beta Distribution for the control response rate used to calculate average operating characteristics. Corresponds to \eqn{\alpha_c}{\alpha c}.
}
\item{beta_c}{
Beta parameter of Beta Distribution for the control response rate used to calculate average operating characteristics. Corresponds to \eqn{\beta_c}{\beta c}.
}
\item{trues}{
Sequence of true control response rates and experimental response rates, at which the Probability to Go will be computed. Default is seq(0,1,0.01) to ensure continuous plots and accurate results.
}
\item{plot}{
Plots yes or no. Default is TRUE.
}
\item{coresnum}{
Number of cores used for parallel computing, in case N_e is a vector. Default is the number of total cores - 1.
}
\item{legend}{
Logical; whether or not to include legend in plot. Default is TRUE.
}
\item{legend.pos}{
Position of legend. Default is "topleft".
}
}
\value{
\item{length(N_e)=1}{
A vector containing the average decision power and average alpha.
}
\item{length(N_e)>1}{
A matrix containing the average decision power and average decision alpha. Every row corresponds to one value of N_e.
}
}

\author{
Elias Laurin Meyer
}
\seealso{
{\code{\link{avg_oc_wr_ne}}},
{\code{\link{oc_rct}}},
{\code{\link{req_resp_rct}}}
}
\examples{
# Setting 3
avg_oc_wr_ne_rct(N_c=25, N_e=25, delta=0.08,
delta_power=0.13, confidence=0.6, alpha_c=15, beta_c=13)

# Setting 4
avg_oc_wr_ne_rct(N_c=25, N_e=25, delta=0.08,
delta_power=0.13, confidence=0.6, alpha_c=15,
beta_c=13, RR_h = 0.5, N_h = 100, w = 0.3)

}
