\name{prime_export}
\alias{prime_export}
\title{Export protocol and results to separate table folders}
\description{Export protocol and results to separate table folders.}
\usage{
prime_export(x, directory)
}
\arguments{
\item{x}{A prime_audit or prime_lure result.}
\item{directory}{A new destination directory. Existing paths are rejected, preserving prior results.}
}
\value{Invisibly, the normalized export path.}
\details{Writes audit_result.rds, 01_methodology_tables, 02_result_tables and an MD5 file manifest. Results include queried IDs and outcomes and should be handled as research data. Staging is finalized by renaming a directory in the same parent. No unqueried outcomes are exported.}
\examples{
ex <- prime_example(40)
a <- prime_audit(ex$pool, prime_oracle(ex$y))
out <- tempfile("prime_example_")
prime_export(a, out)
unlink(out, recursive = TRUE)
}
