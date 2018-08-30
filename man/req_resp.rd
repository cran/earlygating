\name{req_resp}
\alias{req_resp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Required Responders for GO decision Single Arm}
\description{
Function for calculating the minimum required number of responders in the experimental group to make a GO decision in Settings 1 and 2.
}
\usage{
req_resp(N_e, delta, confidence, e_a=0.5, e_b=0.5,
         h_a=0.5, h_b=0.5, RR_h=NULL, N_h=NULL,
         hist_RR_c=NULL, adapt = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
\item{h_a}{
Alpha parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\alpha_h}{\alpha h}. Default is \eqn{\frac{1}{2}}{1/2}.
}
\item{h_b}{
Beta parameter of Beta Prior Distribution for the historical control response rate. Corresponds to \eqn{\beta_h}{\beta h}. Default is \eqn{\frac{1}{2}}{1/2}.
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
\item{adapt}{
Level of adapting of experimental control rate to account for patient selection bias from phase II to phase III. Corresponds to \eqn{\xi}{\xi}. Default is 1, so no adapting.
}
}
\value{
Integer.
}
\author{
Elias Laurin Meyer
}

\seealso{
{\code{\link{avg_oc_wr_ne}}},
{\code{\link{oc}}},
{\code{\link{req_resp_rct}}}
}
\examples{
# Setting 1
req_resp(N_e=50, delta=0.08, confidence=0.6, hist_RR_c=0.5)

# Setting 2
req_resp(N_e=50, delta=0.08, confidence=0.6, RR_h=0.5, N_h = 50)
}

