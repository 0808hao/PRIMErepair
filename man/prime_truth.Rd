\name{prime_truth}
\alias{prime_truth}
\title{Evaluate the complete finite-pool loss gaps}
\description{Evaluate the complete finite-pool loss gaps.}
\usage{
prime_truth(pool, y, epsilon = .002)
}
\arguments{
\item{pool}{A locked object returned by \code{prime_pool()}.}
\item{y}{Complete binary labels. Named labels are matched by unit ID; unnamed labels follow pool row order.}
\item{epsilon}{Nonnegative Brier-loss gap margin, chosen before querying labels. Positive gaps mean the candidate is worse. The default 0.002 is an engineering setting, not a clinically validated margin.}
}
\value{A list with candidate feasibility, exact group gaps, minimum feasible cost and first feasible candidate.}
\details{Retrospective evaluation only. Do not use full-label truth to tune candidate definitions, thresholds or audit acquisition and then describe those choices as held out.}
\examples{
ex <- prime_example(40)
prime_truth(ex$pool, ex$y)
}
