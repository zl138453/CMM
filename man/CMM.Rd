\name{CMM}
\alias{CMM}
\title{Compositional Mediation Model}
\description{
A compositional mediation model for continuous outcome and binary outcomes to deal with mediators that are compositional data
}
\usage{
CMM(data,outcome,med,pred,cov_con=NULL,cov_cat=NULL,weight=NULL,family="identity",boot=5000)
}
\arguments{
  \item{data}{an input dataframe}
  \item{outcome}{column number that locate continuous or binary outcome variable in \code{data}}
  \item{med}{a vector of column numbers that locate the compositional mediators in \code{data}}
  \item{pred}{column number that locate continuous or binary exposure in \code{data}}
  \item{cov_con}{a vector of column numbers that locate the continuous covariates in \code{data} (default NULL).}
  \item{cov_cat}{a vector of column numbers that locate the categorical covariates in \code{data} (default NULL).}
  \item{weight}{column number that locate weights in \code{data} (default NULL).}
  \item{family}{If your outcome variable is continuous, then family="identity"; if your outcome variable is binary, then family="logistic" (default "identity") }
  \item{boot}{Number of bootstrap (default 5000)}

}
\value{
  An object of class \code{CMM}, which is a list with the following components:
  \item{Indirect.effect}{Indirect effects of exposure on an outcome variable (with 95\% bootstrap confidence intervals)}
  \item{Direct.effect}{Direct effects of exposure on an outcome variable (with 95\% bootstrap confidence intervals)}
  \item{Total.effect}{Total effects of exposure on an outcome variable (with 95\% bootstrap confidence intervals)}
  \item{Mediation.effect.plot}{A plot shows mediation effect of exposure on an outcome variables (mediation effect with with 95\% bootstrap confidence intervals)}
  \item{Relative.Effects.plot}{A plot shows relative effect of exposure on an outcome variables (relative effect with with 95\% bootstrap confidence intervals)}

}
\details{
  This code can be used to model with a situation when the mediators are compositional data.
}

\references{
  Lin Z, Zhu S, Cheng J, Lin Q, Lawrence WR, Zhang W, Huang Y, Chen Y, Gao Y. The mediating effect of engagement in physical activity over a 24-hour period on chronic disease and depression: using compositional mediation model. \emph{J Affect Disord}. 2021 Dec 10:S0165-0327(21)01337-9. doi: 10.1016/j.jad.2021.12.019.
}


\examples{
data(CMM_test_data)
result=CMM(uu,1,3:22,2,cov_con=23:24,cov_cat=NULL,weight=NULL,boot=100)
}
