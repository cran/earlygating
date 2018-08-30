\name{avg_oc_wr_ph}
\alias{avg_oc_wr_ph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Average operating characteristics with respect to historic target}
\description{
Function for calculating the average operating characteristics of a single arm bayesian designs for early gating with respect to the historic target.
}
\usage{
avg_oc_wr_ph(N_e, delta, delta_power, confidence, e_a=0.5,
             e_b=0.5, alpha_c, beta_c, trues=seq(0,1,0.01),
             adapt=1, plot = T, legend = T, legend.pos="topleft")
}
%- maybe also 'usage' for other objects documented here.
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
\item{legend}{
Logical; whether or not to include legend in plot. Default is TRUE.
}
\item{legend.pos}{
Position of legend. Default is "topleft".
}
}
\value{
A matrix containing information about the decision power and the decision alpha with respect to p_h.
}

\author{
Elias Laurin Meyer
}

\seealso{
{\code{\link{avg_oc_wr_ne}}},
{\code{\link{avg_oc_wr_ne_rct}}}
}
\examples{
avg_oc_wr_ph(N_e=50, delta=0.08, delta_power=0.13,
confidence=0.6, alpha_c=15, beta_c=13)
}

