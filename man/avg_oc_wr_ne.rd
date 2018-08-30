\name{avg_oc_wr_ne}
\alias{avg_oc_wr_ne}
\title{Single Arm Average Operating Characteristics}
\description{
Function for calculating the average operating characteristics of two single arm bayesian designs for early gating with respect to the sample size in the experimental group and possible historical data.
}
\usage{
avg_oc_wr_ne(N_e, true_RR_c=NULL, delta, delta_power,
             confidence, e_a=0.5, e_b=0.5, h_a=0.5,
             h_b=0.5, RR_h=NULL, N_h=NULL, hist_RR_c=NULL,
             alpha_c, beta_c, trues = seq(0,1,0.001),
             adapt = 1, plot = T, coresnum = NULL,
             legend = T, legend.pos = "topleft")
}

\arguments{
  \item{N_e}{
Sample Size in the experimental group. Can be either a single value or a vector.
}
  \item{true_RR_c}{
Default value is NULL. If specified, will be used in the generated plots, indicating the true achieved decision power and decision type 1 error. If not specified, will be set to either RR_h or hist_RR_c, depending on which was specified by the user.
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
\item{alpha_c}{
Alpha parameter of Beta Distribution for the control response rate used to calculate average operating characteristics. Corresponds to \eqn{\alpha_c}{\alpha c}.
}
\item{beta_c}{
Beta parameter of Beta Distribution for the control response rate used to calculate average operating characteristics. Corresponds to \eqn{\beta_c}{\beta c}.
}
\item{trues}{
Sequence of true control response rates and experimental response rates, at which the Probability to Go will be computed. Default is seq(0,1,0.001) to ensure continuous plots and accurate results.
}
\item{adapt}{
Level of adapting of experimental control rate to account for patient selection bias from phase II to phase III. Corresponds to \eqn{\xi}{\xi}. Default is 1, so no adapting.
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
{\code{\link{avg_oc_wr_ph}}},
{\code{\link{oc}}},
{\code{\link{req_resp}}},
{\code{\link{avg_oc_wr_ne_rct}}}
}
\examples{

# Setting 1
avg_oc_wr_ne(N_e=50, delta=0.08, delta_power=0.13,
confidence=0.6, hist_RR_c=0.5, alpha_c=15, beta_c=13)


# Setting 2
avg_oc_wr_ne(N_e=50, delta=0.08, delta_power=0.13,
confidence=0.6, RR_h=0.5, N_h = 50, alpha_c=15, beta_c=13)

}

