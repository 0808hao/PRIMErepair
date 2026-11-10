\name{prime_audit}
\alias{prime_audit}
\title{Audit and certify minimum input cost}
\description{Audit and certify minimum input cost.}
\usage{
prime_audit(pool, reveal, epsilon = .002, method = "decision_focused",
  budgets = 1, batch = 64L, seed = 1L, alpha = NULL)
}
\arguments{
\item{pool}{A locked object returned by \code{prime_pool()}.}
\item{reveal}{Function taking numeric row indices and returning binary numeric/logical labels in request order. If labels are named, names must equal requested unit IDs. Use a fresh \code{prime_oracle(y)} for each independent retrospective run.}
\item{epsilon}{Nonnegative Brier-loss gap margin, chosen before querying labels. Positive gaps mean the candidate is worse. The default 0.002 is an engineering setting, not a clinically validated margin.}
\item{method}{Exactly one method name from \code{prime_methods()}. Default: frozen decision_focused deterministic-support acquisition.}
\item{budgets}{Distinct prespecified fractions in (0, 1]. Effective caps are sorted unique values of \code{max(1, floor(n * fraction))}; fractions below one label therefore request one label. Caps and batch boundaries are part of the query protocol.}
\item{batch}{Positive integer query batch size. It affects query schedules and stopping counts; keep it fixed for comparisons.}
\item{seed}{Integer seed in [0, .Machine$integer.max]. The caller's RNG kind and state are restored, including after errors.}
\item{alpha}{For statistical methods: prespecified distinct error probabilities strictly between zero and one; default 0.05. Must be NULL for deterministic methods.}
}
\value{A \code{prime_audit} list with decision summaries, bounds, the unique-query log, diagnostics and a reproducibility protocol. No unqueried outcomes are stored. Deterministic bounds include the proxy-completed mean; it is a heuristic, not an unbiased estimator. Statistical bounds are for the final requested cap, for each configuration at its stopping point.}
\details{A feasible candidate has upper bounds at most epsilon (numerical tolerance 1e-12) in every included group. A candidate is excluded if any lower bound exceeds epsilon. Minimum cost is certified when a feasible candidate exists and every strictly cheaper candidate is excluded. Equal-cost ranking is not required.

Statistical modes use the empirical-Bernstein without-replacement confidence sequence of Waudby-Smith and Ramdas (2020), with fixed-proxy residuals and public-cell aggregation where indicated. Bonferroni allocation covers all nonconstant candidate-region constraints. Coverage is simultaneous over time for each alpha separately, not jointly across alpha values, reused pools or future populations. Conflicting statistical intervals are recorded and that alpha falls back to deterministic support bounds. No outcome-based false-certificate correction is performed.

Deterministic runs accept one epsilon. Statistical runs can share a query path across several epsilon/alpha pairs. The decision-cell priority first targets epsilon 0.002 and the largest supplied alpha if present, then the first unresolved configuration. The query log counts all labels acquired; a configuration's labels_used may be lower because it stopped earlier. Prespecify the entire configuration, including ordering, before seeing labels.

On callback failure, a \code{prime_query_error} contains completed_queries, pending_indices and pool_fingerprint. Pending service charges cannot be inferred by the package. No durable resume API is provided in this version.}
\examples{
ex <- prime_example(80)
a <- prime_audit(ex$pool, prime_oracle(ex$y), batch = 16)
print(a)
b <- prime_audit(ex$pool, prime_oracle(ex$y),
  method = "srs_wor_EB", alpha = .05, batch = 16)
print(b)
}
