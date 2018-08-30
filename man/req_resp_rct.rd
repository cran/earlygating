\name{req_resp_rct}
\alias{req_resp_rct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Required Responders for GO decision RCT}
\description{
Function for calculating the minimum required number of responders in the experimental group to make a GO decision in Settings 3 and 4.
}
\usage{
req_resp_rct(N_c, N_e, delta, confidence, e_a=0.5, e_b=0.5,
             c_a=0.5, c_b=0.5, h_a=0.5, h_b=0.5, RR_h=NULL,
             N_h=NULL, w=NULL, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{N_c}{
Sample Size in the control group.
}
\item{N_e}{
Sample Size in the experimental group.
}
\item{delta}{
Required superiority to make a "GO" decision. Corresponds to \eqn{\delta}.
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
Alpha parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\alpha_h}{\alpha h}. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{h_b}{
Beta parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\beta_h}{\beta h}. Default is \eqn{\frac{1}{2}}{1/2}.
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
\item{plot}{
Plots yes or no. Default is TRUE.
}
}
\value{
Matrix containing pairs of successes in control group and respective required successes in experimental group.
}
\author{
Elias Laurin Meyer
}

\seealso{
{\code{\link{avg_oc_wr_ne_rct}}},
{\code{\link{oc_rct}}},
{\code{\link{req_resp}}},
}
\examples{
# Setting 3
req_resp_rct(N_c=25, N_e=25, delta=0.08, confidence=0.6)

# Setting 4
req_resp_rct(N_c=25, N_e=25, delta=0.08, confidence=0.6,
RR_h=0.5, N_h=50, w=0.3)
}

