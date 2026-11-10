# Public API: explicit validation, label accounting and reproducible RNG scope.
.scalar_number <- function(x, name, lower=-Inf, upper=Inf, integer=FALSE) {
  need(is.numeric(x) && length(x)==1L && is.finite(x) && x>=lower && x<=upper &&
    (!integer || x==floor(x)), paste0(name," has an invalid value."))
  invisible(x)
}
.numbers <- function(x, name, lower, upper, strict_lower=FALSE) {
  need(is.numeric(x) && length(x)>0L && all(is.finite(x)) &&
    all(if(strict_lower) x>lower else x>=lower) && all(x<=upper) &&
    !anyDuplicated(x), paste0(name," must contain distinct finite values in the documented range."))
  invisible(x)
}
.names_ok <- function(x, n) is.character(x) && length(x)==n &&
  !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
.align <- function(x, ids, name) {
  need(length(x)==length(ids),paste0(name," length does not match the corresponding rows/columns."))
  if(!is.null(names(x))) {
    need(.names_ok(names(x),length(ids)) && setequal(names(x),ids),
      paste0(name," names must match the corresponding IDs exactly."))
    x <- x[match(ids,names(x))]
  }
  unname(x)
}
.hash <- function(x) {
  f <- tempfile("primerepair_",fileext=".rds")
  on.exit(unlink(f),add=TRUE)
  saveRDS(x,f,compress=FALSE,version=2)
  unname(tools::md5sum(f))
}
.pool_hash <- function(x) {
  x$fingerprint <- NULL
  .hash(unclass(x))
}
.check_pool <- function(pool) {
  need(inherits(pool,"prime_pool") && identical(pool$schema,"PRIMErepair_pool_1"),
    "Use prime_pool() to construct a supported pool.")
  need(identical(pool$fingerprint,.pool_hash(pool)),
    "Pool contents have changed after construction. Create a new pool and a new audit protocol.")
  invisible(pool)
}
.rng_scope <- function(expr) {
  kinds <- RNGkind()
  had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old <- if(had) get(".Random.seed",envir=.GlobalEnv,inherits=FALSE) else NULL
  on.exit({
    do.call(RNGkind,as.list(kinds))
    if(had) assign(".Random.seed",old,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  force(expr)
}

prime_pool <- function(predictions, cost, reference, groups=NULL, strata=NULL,
                       proxy=NULL, ids=NULL, cost_unit="input channels") {
  need(is.matrix(predictions) && is.numeric(predictions) &&
    nrow(predictions)>0L && ncol(predictions)>=2L &&
    all(is.finite(predictions)) && all(predictions>=0 & predictions<=1),
    "predictions must be a finite numeric probability matrix with at least two columns.")
  P <- predictions; n <- nrow(P); K <- ncol(P); cn <- colnames(P)
  need(.names_ok(cn,K),"Prediction columns need unique nonempty candidate names.")
  need(is.character(reference) && length(reference)==1L && !is.na(reference) && reference %in% cn,
    "reference must name exactly one prediction column.")
  need(is.numeric(cost),"cost must be numeric.")
  cost <- .align(cost,cn,"cost")
  need(all(is.finite(cost)) && all(cost>=0),"Costs must be finite and nonnegative.")
  ref <- match(reference,cn)
  need(cost[ref]==max(cost),"The reference must have maximal cost in this version; remove more expensive candidates.")
  need(is.character(cost_unit) && length(cost_unit)==1L && !is.na(cost_unit) && nzchar(cost_unit),
    "cost_unit must be a nonempty character scalar.")
  if(is.null(ids)) ids <- if(is.null(rownames(P))) as.character(seq_len(n)) else rownames(P)
  need(is.atomic(ids) && is.null(dim(ids)),"ids must be a vector.")
  ids <- as.character(ids)
  need(.names_ok(ids,n),"Unit IDs must be unique and nonempty, with one per row.")
  if(!is.null(rownames(P))) need(identical(rownames(P),ids),"Prediction row names and ids differ in value or order.")
  if(is.null(groups)) groups <- matrix(1,n,1,dimnames=list(NULL,"All"))
  need(is.matrix(groups) && (is.numeric(groups)||is.logical(groups)) &&
    nrow(groups)==n && ncol(groups)>0L && !anyNA(groups) &&
    all(groups %in% c(0,1)) && all(colSums(groups)>0L),
    "groups must be a binary membership matrix with matching rows and no empty groups.")
  need(.names_ok(colnames(groups),ncol(groups)),"Group columns need unique nonempty names.")
  if(!is.null(rownames(groups))) need(identical(rownames(groups),ids),"Group row names and ids differ in value or order.")
  storage.mode(groups) <- "double"
  if(is.null(strata)) strata <- rep("All",n)
  need(is.atomic(strata) && is.null(dim(strata)),"strata must be a vector.")
  strata <- .align(strata,ids,"strata")
  need(!anyNA(strata) && all(nzchar(as.character(strata))),"strata cannot contain missing or empty values.")
  if(!is.null(proxy)) {
    need(is.numeric(proxy) && is.null(dim(proxy)),"proxy must be a numeric vector.")
    proxy <- .align(proxy,ids,"proxy")
    need(all(is.finite(proxy)) && all(proxy>=0 & proxy<=1),"proxy must be finite probabilities.")
  }
  # Stable cost ordering; ties retain input order except that the reference is last.
  o <- order(cost,seq_len(K)==ref,seq_len(K))
  pub <- make_public(P[,o,drop=FALSE],groups,strata,cost[o],proxy,ids)
  z <- list(schema="PRIMErepair_pool_1",public=pub,reference=reference,cost_unit=cost_unit,
    candidates=data.frame(candidate=cn[o],cost=cost[o],input_column=o,
      reference=cn[o]==reference,stringsAsFactors=FALSE))
  class(z) <- "prime_pool"
  z$fingerprint <- .pool_hash(z)
  z
}

prime_methods <- function() {
  data.frame(method=c(METHODS,CS_METHODS),
    mode=c(rep("deterministic",length(METHODS)),rep("statistical",length(CS_METHODS))),
    randomized=c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,rep(TRUE,4),rep(TRUE,4)),
    note=c("Frozen decision-focused acquisition","Frozen width-reduction anchor",
      "Two frozen virtual lanes with shared label charging","Static proxy pairwise ranking",
      "Locked-reference risk >= 0.20 first","Uniform without replacement",
      "Proportional interleaving of supplied strata","Hara Eq.11; ensemble surrogate adaptation",
      "Hara Eq.11; reference surrogate adaptation","Reference expected-Brier acquisition adaptation",
      "Uniform query order; same support stopping","Raw loss-gap WoR confidence sequences",
      "Fixed-proxy residual WoR confidence sequences","Public-cell residual confidence sequences",
      "Decision-focused cell choice; uniform within each public cell"),stringsAsFactors=FALSE)
}

prime_oracle <- function(y) {
  need((is.numeric(y)||is.logical(y)) && is.null(dim(y)) && length(y)>0L &&
    !anyNA(y) && all(y %in% c(0,1)),"y must be a nonempty binary numeric/logical vector.")
  if(!is.null(names(y))) need(.names_ok(names(y),length(y)),"Named labels need unique nonempty unit IDs.")
  y <- y + 0L; used <- rep(FALSE,length(y)); force(y)
  function(indices) {
    need(is.numeric(indices) && length(indices)>0L && !anyNA(indices) &&
      all(is.finite(indices)) && all(indices==floor(indices)) &&
      all(indices>=1 & indices<=length(y)) && !anyDuplicated(indices) && !any(used[indices]),
      "Duplicate or invalid oracle query. Use a fresh oracle for each independent audit.")
    used[indices] <<- TRUE
    y[indices]
  }
}

.query_gateway <- function(pool,reveal) {
  need(is.function(reveal),"reveal must be a function(indices); use prime_oracle(y) for retrospective data.")
  e <- new.env(parent=emptyenv())
  e$known <- rep(NA_real_,pool$public$n); e$indices <- integer(); e$labels <- numeric()
  e$batch <- integer(); e$calls <- 0L; e$pending <- integer()
  query <- function(indices) {
    need(is.numeric(indices) && length(indices)>0L && all(is.finite(indices)) &&
      all(indices==floor(indices)) && all(indices>=1 & indices<=pool$public$n) &&
      !anyDuplicated(indices) && all(is.na(e$known[indices])),"Internal duplicate or invalid query.")
    e$pending <- indices
    yy <- reveal(indices)
    need((is.numeric(yy)||is.logical(yy)) && is.null(dim(yy)) &&
      length(yy)==length(indices) && !anyNA(yy) && all(yy %in% c(0,1)),
      "reveal must return one binary numeric/logical label per requested index, in request order.")
    if(!is.null(names(yy))) need(identical(names(yy),pool$public$ids[indices]),
      "Named labels returned by reveal do not match the requested unit IDs in order.")
    yy <- as.numeric(yy); e$calls <- e$calls+1L
    e$known[indices] <- yy; e$indices <- c(e$indices,indices); e$labels <- c(e$labels,yy)
    e$batch <- c(e$batch,rep(e$calls,length(indices))); e$pending <- integer()
    yy
  }
  log <- function() data.frame(query=seq_along(e$indices),batch=e$batch,
    index=e$indices,id=pool$public$ids[e$indices],label=e$labels,stringsAsFactors=FALSE)
  guard <- function(expr) tryCatch(force(expr),error=function(err) {
    stop(structure(list(message=paste0("Audit stopped: ",conditionMessage(err)),call=NULL,
      completed_queries=log(),pending_indices=e$pending,pool_fingerprint=pool$fingerprint,
      cause=err),class=c("prime_query_error","error","condition")))
  })
  list(query=query,log=log,guard=guard)
}
.settings <- function(epsilon,budgets,batch,seed,method) {
  .numbers(epsilon,"epsilon",0,Inf)
  .numbers(budgets,"budgets",0,1,TRUE)
  .scalar_number(batch,"batch",1,.Machine$integer.max,TRUE)
  .scalar_number(seed,"seed",0,.Machine$integer.max,TRUE)
  need(is.character(method) && length(method)==1L && !is.na(method) &&
    method %in% c(METHODS,CS_METHODS),"Unknown method; see prime_methods().")
}
.rename_cost <- function(x) {
  names(x)[names(x)=="selected_input_channels"] <- "selected_cost"
  x
}
.protocol <- function(pool,method,epsilon,budgets,batch,seed,alpha,mode) {
  list(package_version="0.1.0",kernel_version="Step133_v0.4.1",mode=mode,method=method,
    epsilon=epsilon,alpha=alpha,budget_fractions=sort(budgets),batch=batch,seed=seed,
    pool_fingerprint=pool$fingerprint,n=pool$public$n,reference=pool$reference,
    cost_unit=pool$cost_unit,candidates=pool$candidates,
    groups=data.frame(group=colnames(pool$public$G),n=pool$public$ng),
    rng=c("Mersenne-Twister","Inversion","Rejection"),
    R_version=as.character(getRversion()),numeric_tolerance=TOL,
    scope="Locked finite audit pool, binary Brier loss, specified candidate library and groups. No population or clinical safety guarantee.")
}

prime_audit <- function(pool, reveal, epsilon=.002, method="decision_focused",
                        budgets=1, batch=64L, seed=1L, alpha=NULL) {
  .check_pool(pool); .settings(epsilon,budgets,batch,seed,method)
  statistical <- method %in% CS_METHODS
  if(statistical) {
    if(is.null(alpha)) alpha <- .05
    .numbers(alpha,"alpha",0,1,TRUE); need(all(alpha<1),"alpha must be strictly below one.")
  } else {
    need(length(epsilon)==1L,"Deterministic audits use one fixed epsilon per run.")
    need(is.null(alpha),"alpha is only used by statistical methods; omit it for deterministic auditing.")
  }
  gateway <- .query_gateway(pool,reveal)
  path <- gateway$guard(.rng_scope({
    if(statistical) cs_path133(pool$public,gateway$query,method,seed,epsilon,alpha,budgets,batch)
    else audit_path(pool$public,gateway$query,epsilon,method,seed,budgets,batch)
  }))
  if(statistical) {
    s <- .rename_cost(path$runs); bounds <- path$final_bounds
    diag <- list(regions=path$regions,allocated_constraints=path$allocated_constraints,
      conflicts=data.frame(alpha=alpha,conflicts=path$conflicts),
      guarantee="For each prespecified alpha separately: simultaneous finite-pool coverage across candidates, groups and query times. Not joint coverage across alpha values or repeated audits.")
  } else {
    s <- .rename_cost(path$snapshots); s$epsilon <- epsilon
    s$guarantee <- "Deterministic finite-pool support certificate"
    bounds <- path$bounds
    diag <- list(query_composition=path$query_composition,events=path$events,trace=path$trace,
      virtual_labels=path$virtual_labels,guarantee=s$guarantee[1L])
  }
  s$method <- method; s$cost_unit <- pool$cost_unit
  result <- list(summary=s,bounds=bounds,queries=gateway$log(),diagnostics=diag,
    protocol=.protocol(pool,method,epsilon,budgets,batch,seed,alpha,
      if(statistical) "statistical" else "deterministic"))
  result$diagnostics$total_unique_labels <- nrow(result$queries)
  result$diagnostics$multi_configuration_cost <- "labels_used is the stopping count for that configuration along the shared query path; queries records all labels actually acquired across configurations."
  class(result) <- "prime_audit"
  result
}

prime_lure <- function(pool, reveal, budgets=.1, epsilon=.002,
                       method="hara_ensemble", seed=1L) {
  .check_pool(pool); .settings(epsilon,budgets,1L,seed,method)
  need(method %in% LITERATURE_METHODS,"LURE supports hara_ensemble, hara_reference, active_testing_reference or srs_lure.")
  gateway <- .query_gateway(pool,reveal)
  x <- gateway$guard(.rng_scope({
    draw <- weighted_order(pool$public,method,seed)
    fixed_budget_lure(pool$public,gateway$query,draw,epsilon,budgets)
  }))
  x$summary <- .rename_cost(x$decisions); x$decisions <- NULL
  x$summary$method <- method; x$summary$cost_unit <- pool$cost_unit
  x$queries <- gateway$log()
  x$protocol <- .protocol(pool,method,epsilon,budgets,NA_integer_,seed,NULL,"fixed_budget_point_estimation")
  x$protocol$scope <- paste(x$protocol$scope,"LURE is a fixed-budget estimator; plugin choices have NO certificate and no optional-stopping guarantee.")
  class(x) <- "prime_lure"
  x
}

prime_truth <- function(pool, y, epsilon=.002) {
  .check_pool(pool); .scalar_number(epsilon,"epsilon",0,Inf)
  need((is.numeric(y)||is.logical(y)) && is.null(dim(y)),"y must be a binary vector.")
  y <- .align(y,pool$public$ids,"y")
  need(!anyNA(y) && all(y %in% c(0,1)),"y must contain only zero and one.")
  tt <- truth(pool$public,y,epsilon)
  out <- pool$candidates
  out$feasible <- tt$feasible; out$worst_group_gap <- apply(tt$delta,1,max)
  out$minimum_cost_feasible <- tt$feasible & out$cost==tt$min_cost
  b <- expand.grid(candidate=colnames(pool$public$P),group=colnames(pool$public$G),stringsAsFactors=FALSE)
  b$mean_brier_gap <- as.vector(tt$delta)
  list(candidates=out,group_gaps=b,minimum_cost=tt$min_cost,
    first_feasible=colnames(pool$public$P)[tt$best],epsilon=epsilon,
    purpose="Retrospective full-label evaluation only; never an input to acquisition.")
}

prime_example <- function(n=240L,seed=42L) {
  .scalar_number(n,"n",4,1000000,TRUE); .scalar_number(seed,"seed",0,.Machine$integer.max,TRUE)
  .rng_scope({
    set.seed(seed); p <- stats::plogis(stats::rnorm(n,-1,1)); y <- stats::rbinom(n,1,p)
    predictions <- cbind(cheap=pmin(.99,p+.16),middle=pmin(.99,p+.025),reference=p)
    ids <- paste0("sim_",seq_len(n)); names(y) <- ids
    group <- rep(c("A","B"),length.out=n)
    pool <- prime_pool(predictions,c(1,2,3),"reference",
      groups=cbind(All=1,A=group=="A",B=group=="B"),strata=group,ids=ids)
    list(pool=pool,y=y,description="Synthetic engineering example; not clinical evidence.")
  })
}

print.prime_pool <- function(x,...) {
  cat("PRIMErepair locked finite pool\n",x$public$n," units; ",x$public$K,
    " candidates; ",ncol(x$public$G)," groups\n",sep="")
  cat("Reference:",x$reference,"| Cost unit:",x$cost_unit,"\n")
  print(x$candidates,row.names=FALSE); invisible(x)
}
print.prime_audit <- function(x,...) {
  cat("PRIMErepair ",x$protocol$mode," audit | ",x$protocol$method,"\n",sep="")
  cat("Unique labels acquired:",nrow(x$queries),"/",x$protocol$n,"\n")
  s <- x$summary[x$summary$cap_labels==max(x$summary$cap_labels),,drop=FALSE]
  print(s[,intersect(c("epsilon","alpha","labels_used","selected_candidate","selected_cost",
    "certified_minimum_cost","status"),names(s)),drop=FALSE],row.names=FALSE)
  cat("Scope: fixed pool and locked candidate library.\n"); invisible(x)
}
print.prime_lure <- function(x,...) {
  cat("PRIMErepair fixed-budget LURE point estimates | NO CERTIFICATE\n")
  print(x$summary,row.names=FALSE); invisible(x)
}

prime_export <- function(x,directory) {
  need(inherits(x,"prime_audit")||inherits(x,"prime_lure"),"x must be a prime_audit or prime_lure result.")
  need(is.character(directory) && length(directory)==1L && !is.na(directory) && nzchar(directory),
    "directory must be a nonempty path.")
  need(!file.exists(directory),"Output path already exists. Use a new directory to preserve earlier results.")
  parent <- dirname(directory)
  if(!dir.exists(parent)) need(dir.create(parent,recursive=TRUE),"Cannot create parent directory.")
  stage <- tempfile(".primerepair_export_",tmpdir=parent)
  need(dir.create(stage),"Cannot create export staging directory.")
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  m <- file.path(stage,"01_methodology_tables"); r <- file.path(stage,"02_result_tables")
  dir.create(m); dir.create(r)
  write_csv <- function(obj,path) utils::write.csv(obj,path,row.names=FALSE,na="",fileEncoding="UTF-8")
  saveRDS(x,file.path(stage,"audit_result.rds"))
  saveRDS(x$protocol,file.path(m,"protocol.rds"))
  write_csv(x$protocol$candidates,file.path(m,"candidates.csv"))
  write_csv(x$protocol$groups,file.path(m,"groups.csv"))
  simple <- x$protocol[!vapply(x$protocol,is.data.frame,TRUE)]
  write_csv(data.frame(item=names(simple),value=vapply(simple,function(v) paste(v,collapse="; "),"")),
    file.path(m,"protocol.csv"))
  write_csv(x$summary,file.path(r,"decisions.csv"))
  write_csv(x$queries,file.path(r,"queried_labels.csv"))
  if(inherits(x,"prime_audit")) {
    write_csv(x$bounds,file.path(r,"bounds.csv"))
    for(nm in names(x$diagnostics)) if(is.data.frame(x$diagnostics[[nm]]) && ncol(x$diagnostics[[nm]]))
      write_csv(x$diagnostics[[nm]],file.path(r,paste0(nm,".csv")))
  } else {
    write_csv(x$estimates,file.path(r,"lure_estimates.csv"))
    write_csv(x$diagnostics,file.path(r,"lure_diagnostics.csv"))
  }
  writeLines(c(x$protocol$scope,
    "queried_labels.csv contains supplied unit IDs and queried outcomes; handle as research data.",
    "proxy_completed_mean, when present, is an acquisition heuristic, not an unbiased point estimate."),
    file.path(stage,"README.txt"))
  files <- list.files(stage,recursive=TRUE,full.names=TRUE)
  write_csv(data.frame(file=substring(files,nchar(stage)+2L),md5=unname(tools::md5sum(files))),
    file.path(stage,"file_manifest.csv"))
  need(file.rename(stage,directory),"Could not finalize the export directory.")
  invisible(normalizePath(directory))
}
