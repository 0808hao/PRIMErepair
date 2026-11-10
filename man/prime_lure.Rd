\name{prime_lure}
\alias{prime_lure}
\title{Compute separate fixed-budget LURE estimates}
\description{Compute separate fixed-budget LURE estimates.}
\usage{
prime_lure(pool, reveal, budgets = .1, epsilon = .002,
  method = "hara_ensemble", seed = 1L)
}
\arguments{
\item{pool}{A locked object returned by \code{prime_pool()}.}
\item{reveal}{Function taking numeric row indices and returning binary numeric/logical labels in request order. If labels are named, names must equal requested unit IDs. Use a fresh \code{prime_oracle(y)} for each independent retrospective run.}
\item{budgets}{Distinct prespecified fractions in (0, 1]. Effective caps are sorted unique values of \code{max(1, floor(n * fraction))}; fractions below one label therefore request one label. Caps and batch boundaries are part of the query protocol.}
\item{epsilon}{Nonnegative Brier-loss gap margin, chosen before querying labels. Positive gaps mean the candidate is worse. The default 0.002 is an engineering setting, not a clinically validated margin.}
\item{method}{hara_ensemble, hara_reference, active_testing_reference or srs_lure.}
\item{seed}{Integer seed in [0, .Machine$integer.max]. The caller's RNG kind and state are restored, including after errors.}
}
\value{A \code{prime_lure} list containing group-gap and overall-risk estimates, plugin decisions, sampling diagnostics, query log and protocol. All decisions explicitly have no certificate.}
\details{Prespecified fixed budgets use conditional without-replacement query probabilities and LURE Eq.2, with known group-size scaling for group gaps. Weights are not normalized and estimates are not clipped. Overall risk estimates are repeated on each group's row and remain overall, not group-specific risks. Query weights have a fixed relative floor of 1e-10 times their mean after scaling; modified-unit counts are reported. No optional-stopping or post-selection confidence guarantee is claimed. This is a task adaptation of Hara et al. (2024) and Active Testing, not a replication of their original experiments.}
\examples{
ex <- prime_example(80)
z <- prime_lure(ex$pool, prime_oracle(ex$y), budgets = c(.1, .5))
print(z)
}
