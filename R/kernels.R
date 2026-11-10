# Frozen algorithm kernels extracted from Step133 v0.4.1.
# Numerical function bodies are retained; public validation lives in api.R.
TOL <- 1e-12
cfg <- list(probability_floor_relative=1e-10, primary_epsilon=.002)
LEGACY_METHODS <- c("decision_focused","shared_width","portfolio_replay","pairwise_proxy","highrisk_first","srs","age_stratified")
LITERATURE_METHODS <- c("hara_ensemble","hara_reference","active_testing_reference","srs_lure")
METHODS <- c(LEGACY_METHODS,LITERATURE_METHODS)
CS_METHODS <- c("srs_wor_EB","srs_residual_wor_EB","stratified_residual_wor_EB","decision_cell_residual_wor_EB")

need <- function (ok, msg) 
if (!isTRUE(ok)) stop(msg, call. = FALSE)

bind <- function (x) 
if (!length(x)) data.frame() else do.call(rbind, x)

make_public <- function (P, G, strata, costs, proxy = NULL, ids = NULL) 
{
    P <- as.matrix(P)
    G <- as.matrix(G)
    n <- nrow(P)
    K <- ncol(P)
    need(n > 0 && K >= 2 && is.numeric(P) && all(is.finite(P)) && all(P >= 0 & P <= 1), "\u5019\u9009\u9884\u6d4b\u5fc5\u987b\u662f\u6709\u96500-1\u6982\u7387\u77e9\u9635")
    need(!is.null(colnames(P)) && !anyDuplicated(colnames(P)), "\u5019\u9009\u540d\u79f0\u7f3a\u5931\u6216\u91cd\u590d")
    need(nrow(G) == n && ncol(G) > 0 && all(G %in% c(0, 1)) && all(colSums(G) > 0), "\u4e9a\u7ec4\u77e9\u9635\u5f02\u5e38\u6216\u6709\u7a7a\u4e9a\u7ec4")
    need(length(strata) == n && !anyNA(strata), "\u62bd\u6837\u5206\u5c42\u7f3a\u5931")
    need(length(costs) == K && all(is.finite(costs)) && all(diff(costs) >= 0), "\u5019\u9009\u987b\u4e8b\u5148\u6309\u6210\u672c\u548c\u56fa\u5b9a\u4f18\u5148\u987a\u5e8f\u6392\u5e8f")
    if (is.null(ids)) 
        ids <- as.character(seq_len(n))
    need(length(ids) == n && !anyNA(ids) && !anyDuplicated(ids), "\u5ba1\u8ba1\u5355\u4f4dID\u7f3a\u5931\u6216\u91cd\u590d")
    if (is.null(proxy)) 
        proxy <- P[, K]
    need(length(proxy) == n && all(is.finite(proxy)) && all(proxy >= 0 & proxy <= 1), "\u4ee3\u7406\u6982\u7387\u4e0d\u5408\u6cd5")
    D0 <- sweep(P^2, 1, P[, K]^2, "-")
    D1 <- sweep((P - 1)^2, 1, (P[, K] - 1)^2, "-")
    lo <- pmin(D0, D1)
    hi <- pmax(D0, D1)
    pr <- D0 + (D1 - D0) * proxy
    pair_score <- numeric(n)
    for (a in seq_len(K - 1L)) for (b in seq.int(a + 1L, K)) {
        d0 <- P[, a]^2 - P[, b]^2
        d1 <- (P[, a] - 1)^2 - (P[, b] - 1)^2
        pair_score <- pair_score + (1 - proxy) * abs(d0) + proxy * abs(d1)
    }
    list(P = P, G = G, strata = as.character(strata), costs = costs, ids = as.character(ids), n = n, K = K, ng = colSums(G), 
        D0 = D0, D1 = D1, lo = lo, hi = hi, pr = pr, pair_score = pair_score, initL = sweep(crossprod(lo, G), 2, 
            colSums(G), "/"), initU = sweep(crossprod(hi, G), 2, colSums(G), "/"), initE = sweep(crossprod(pr, 
            G), 2, colSums(G), "/"))
}

query_scores <- function (pub, L, U, E, eps, method, dec) 
{
    s <- numeric(pub$n)
    for (k in dec$contenders) {
        unresolved <- which(U[k, ] > eps + TOL)
        if (method == "shared_width") {
            for (g in unresolved) s <- s + (pub$hi[, k] - pub$lo[, k]) * pub$G[, g]/pub$ng[g]
        }
        else {
            bad <- which(E[k, ] > eps)
            if (length(bad)) {
                g <- bad[which.max((E[k, bad] - eps)/pmax(U[k, bad] - L[k, bad], 1e-12))]
                gap <- max(eps - L[k, g], 0.05 * (U[k, g] - L[k, g]), 1e-12)
                s <- s + (pub$pr[, k] - pub$lo[, k]) * pub$G[, g]/(pub$ng[g] * gap)
            }
            else {
                for (g in unresolved) {
                  gap <- max(U[k, g] - eps, 0.05 * (U[k, g] - L[k, g]), 1e-12)
                  s <- s + (pub$hi[, k] - pub$pr[, k]) * pub$G[, g]/(pub$ng[g] * gap)
                }
            }
        }
    }
    need(all(is.finite(s)), "\u67e5\u8be2\u4f18\u5148\u5206\u6570\u975e\u6709\u9650")
    s
}

stratified_order <- function (strata) 
{
    n <- length(strata)
    score <- numeric(n)
    for (ii in split(seq_len(n), strata)) score[ii[sample.int(length(ii))]] <- (seq_along(ii) - 0.5)/length(ii)
    order(score, stats::runif(n))
}

new_oracle <- function (y) 
{
    need(all(y %in% c(0, 1)), "oracle\u7ed3\u5c40\u4e0d\u662f0/1")
    used <- rep(FALSE, length(y))
    function(indices) {
        need(is.numeric(indices) && length(indices) > 0 && !anyNA(indices) && all(indices == as.integer(indices)) && 
            all(indices >= 1 & indices <= length(y)) && !anyDuplicated(indices) && !any(used[indices]), "\u91cd\u590d\u6536\u8d39\u6216\u975e\u6cd5\u67e5\u8be2\u7d22\u5f15")
        used[indices] <<- TRUE
        y[indices]
    }
}

decision <- function (S, eps, costs) 
{
    feasible <- apply(S$U <= eps + TOL, 1, all)
    infeasible <- apply(S$L > eps + TOL, 1, any)
    need(!any(feasible & infeasible), "\u8bc1\u4e66\u72b6\u6001\u51b2\u7a81")
    unknown <- !(feasible | infeasible)
    best <- if (any(feasible)) 
        which(feasible)[1L]
    else NA_integer_
    contenders <- if (is.na(best)) 
        which(unknown)
    else which(unknown & costs < costs[best])
    certified <- !is.na(best) && !length(contenders)
    preferred <- !is.na(best) && !any(unknown[seq_len(best - 1L)])
    list(feasible = feasible, infeasible = infeasible, unknown = unknown, best = best, contenders = contenders, 
        certified = certified, preferred = preferred, status = if (certified) "CERTIFIED_MINIMUM_COST" else if (!is.na(best)) "FEASIBLE_COST_UNRESOLVED" else "UNRESOLVED")
}

initial_state <- function (pub) 
list(L = pub$initL, U = pub$initU, E = pub$initE)

update_state <- function (S, pub, ids, y) 
{
    if (!length(ids)) 
        return(S)
    zz <- pub$D0[ids, , drop = FALSE] + (pub$D1[ids, , drop = FALSE] - pub$D0[ids, , drop = FALSE]) * y
    gg <- pub$G[ids, , drop = FALSE]
    S$L <- S$L + sweep(crossprod(zz - pub$lo[ids, , drop = FALSE], gg), 2, pub$ng, "/")
    S$U <- S$U + sweep(crossprod(zz - pub$hi[ids, , drop = FALSE], gg), 2, pub$ng, "/")
    S$E <- S$E + sweep(crossprod(zz - pub$pr[ids, , drop = FALSE], gg), 2, pub$ng, "/")
    need(all(S$L <= S$U + 1e-09), "\u754c\u9650\u987a\u5e8f\u9519\u8bef")
    S
}

fixed_order_131 <- function (pub, method, seed) 
{
    set.seed(seed)
    if (method == "srs") 
        return(sample.int(pub$n))
    if (method == "age_stratified") 
        return(stratified_order(pub$strata))
    if (method == "pairwise_proxy") 
        return(order(-pub$pair_score, seq_len(pub$n)))
    if (method == "highrisk_first") {
        return(order(-(pub$P[, pub$K] >= 0.2), stats::runif(pub$n)))
    }
    NULL
}

audit_path <- function (pub, reveal, eps, method, seed, budget_fractions, batch, observer = NULL) 
{
    need(method %in% METHODS && length(eps) == 1L && is.finite(eps) && eps >= 0, "\u7b56\u7565\u6216\u5bb9\u5fcd\u5ea6\u65e0\u6548")
    need(length(batch) == 1L && batch >= 1 && batch == as.integer(batch), "batch\u5fc5\u987b\u4e3a\u6b63\u6574\u6570")
    need(all(budget_fractions > 0 & budget_fractions <= 1), "\u9884\u7b97\u5fc5\u987b\u5728(0,1]")
    caps <- sort(unique(pmax(1L, as.integer(floor(pub$n * budget_fractions)))))
    anchor_methods <- if (method == "portfolio_replay") 
        c("decision_focused", "shared_width")
    else method
    lanes <- lapply(anchor_methods, function(m) list(method = m, S = initial_state(pub), used = rep(FALSE, pub$n), 
        count = 0L, order = fixed_order(pub, m, seed), pending = integer()))
    S <- initial_state(pub)
    dec <- decision(S, eps, pub$costs)
    known <- rep(NA_integer_, pub$n)
    paid <- 0L
    lane_id <- 1L
    iterations <- 0L
    queried <- integer()
    events <- list()
    trace <- list()
    snapshots <- list()
    bounds <- list()
    record_events <- function(before, after, at) {
        out <- list()
        for (k in which(after$infeasible & !before$infeasible)) {
            g <- which.max(S$L[k, ] - eps)
            out[[length(out) + 1L]] <- data.frame(candidate = colnames(pub$P)[k], event = "EXCLUDED", group = colnames(pub$G)[g], 
                labels_used = at)
        }
        for (k in which(after$feasible & !before$feasible)) out[[length(out) + 1L]] <- data.frame(candidate = colnames(pub$P)[k], 
            event = "FEASIBLE", group = "ALL_INCLUDED_GROUPS", labels_used = at)
        out
    }
    if (!is.null(observer)) 
        observer(S, queried)
    for (cap in caps) {
        while (paid < cap && !dec$certified) {
            iterations <- iterations + 1L
            need(iterations <= 4L * pub$n + length(caps), "\u8c03\u5ea6\u672a\u53d6\u5f97\u8fdb\u5c55")
            lane <- lanes[[lane_id]]
            if (!length(lane$pending)) {
                ld <- decision(lane$S, eps, pub$costs)
                need(!ld$certified, "\u865a\u62df\u7b56\u7565\u5df2\u8ba4\u8bc1\u4f46\u5168\u5c40\u4ecd\u672a\u8ba4\u8bc1\uff1a\u72b6\u6001\u4e00\u81f4\u6027\u9519\u8bef")
                count <- min(as.integer(batch), pub$n - lane$count)
                need(count > 0, "\u865a\u62df\u7b56\u7565\u5df2\u8017\u5c3d\u6807\u7b7e\u4f46\u672a\u8ba4\u8bc1")
                if (!is.null(lane$order)) 
                  ids <- lane$order[seq.int(lane$count + 1L, lane$count + count)]
                else {
                  score <- query_scores(pub, lane$S$L, lane$S$U, lane$S$E, eps, lane$method, ld)
                  score[lane$used] <- -Inf
                  ids <- head(order(-score, seq_len(pub$n)), count)
                }
                lane$pending <- ids
            }
            fresh <- lane$pending[is.na(known[lane$pending])]
            buy <- head(fresh, cap - paid)
            if (length(buy)) {
                yy <- reveal(buy)
                need(length(yy) == length(buy) && all(yy %in% c(0, 1)), "oracle\u8fd4\u56de\u975e\u6cd5\u6807\u7b7e")
                need(all(is.na(known[buy])) && !anyDuplicated(buy), "\u91cd\u590d\u6536\u8d39")
                known[buy] <- yy
                S <- update_state(S, pub, buy, yy)
                paid <- paid + length(buy)
                queried <- c(queried, buy)
                if (paid == pub$n) {
                  dd <- pub$D0 + (pub$D1 - pub$D0) * known
                  S$L <- S$U <- S$E <- sweep(crossprod(dd, pub$G), 2, pub$ng, "/")
                }
                old <- dec
                dec <- decision(S, eps, pub$costs)
                events <- c(events, record_events(old, dec, paid))
                trace[[length(trace) + 1L]] <- data.frame(labels_used = paid, active_lane = anchor_methods[lane_id], 
                  status = dec$status, best = dec$best, remaining_cheaper_candidates = length(dec$contenders))
                if (!is.null(observer)) 
                  observer(S, queried)
            }
            if (all(!is.na(known[lane$pending]))) {
                ids <- lane$pending
                lane$S <- update_state(lane$S, pub, ids, known[ids])
                lane$used[ids] <- TRUE
                lane$count <- lane$count + length(ids)
                lane$pending <- integer()
                lanes[[lane_id]] <- lane
                lane_id <- lane_id%%length(lanes) + 1L
            }
            else lanes[[lane_id]] <- lane
        }
        chosen <- dec$best
        snapshots[[length(snapshots) + 1L]] <- data.frame(cap_labels = cap, cap_fraction = cap/pub$n, labels_used = paid, 
            fraction_used = paid/pub$n, selected_index = chosen, selected_candidate = if (is.na(chosen)) 
                NA_character_
            else colnames(pub$P)[chosen], selected_input_channels = if (is.na(chosen)) 
                NA_real_
            else pub$costs[chosen], status = dec$status, certified_minimum_cost = dec$certified, fixed_preference_certified = dec$preferred, 
            remaining_cheaper_candidates = length(dec$contenders))
        b <- expand.grid(candidate = colnames(pub$P), group = colnames(pub$G), stringsAsFactors = FALSE)
        b$cap_labels <- cap
        b$labels_used <- paid
        b$epsilon <- eps
        b$lower <- as.vector(S$L)
        b$upper <- as.vector(S$U)
        b$proxy_completed_mean <- as.vector(S$E)
        bounds[[length(bounds) + 1L]] <- b
    }
    composition <- bind(lapply(seq_len(ncol(pub$G)), function(g) {
        nq <- if (length(queried)) 
            sum(pub$G[queried, g])
        else 0
        data.frame(group = colnames(pub$G)[g], group_n = pub$ng[g], queried_in_group = nq, fraction_of_group_queried = nq/pub$ng[g], 
            share_of_queries = if (paid) 
                nq/paid
            else NA_real_)
    }))
    list(snapshots = bind(snapshots), bounds = bind(bounds), events = bind(events), trace = bind(trace), queried_indices = queried, 
        query_composition = composition, finalS = S, virtual_labels = sum(vapply(lanes, function(v) v$count, 0L)), 
        unique_labels = paid)
}

truth <- function (pub, y, eps) 
{
    d <- pub$D0 + (pub$D1 - pub$D0) * y
    delta <- sweep(crossprod(d, pub$G), 2, pub$ng, "/")
    feasible <- apply(delta <= eps + TOL, 1, all)
    need(any(feasible), "\u9884\u5b9a\u53c2\u8003\u5019\u9009\u5e94\u59cb\u7ec8\u53ef\u884c")
    best <- which(feasible)[1L]
    list(delta = delta, feasible = feasible, best = best, min_cost = min(pub$costs[feasible]))
}

query_weights <- function (pub, method) 
{
    p <- pub$P
    n <- pub$n
    K <- pub$K
    if (method == "srs_lure") 
        return(list(w = rep(1, n), floored = 0L, min_raw = 1, proxy = "none"))
    pi <- if (method == "hara_ensemble") 
        rowMeans(p)
    else p[, K]
    if (method == "active_testing_reference") {
        w <- (1 - pi) * p[, K]^2 + pi * (p[, K] - 1)^2
    }
    else {
        need(method %in% c("hara_ensemble", "hara_reference"), "\u672a\u77e5\u6587\u732e\u67e5\u8be2\u89c4\u5219")
        w <- numeric(n)
        for (k in seq_len(K - 1L)) for (j in seq.int(k + 1L, K)) w <- w + (1 - pi) * abs(p[, k]^2 - p[, j]^2) + 
            pi * abs((p[, k] - 1)^2 - (p[, j] - 1)^2)
    }
    need(all(is.finite(w)) && all(w >= 0), "\u6587\u732e\u67e5\u8be2\u6743\u91cd\u975e\u6cd5")
    raw_min <- min(w)
    if (max(w) == 0) 
        return(list(w = rep(1, n), floored = n, min_raw = 0, proxy = "all_zero_uniform_fallback"))
    w <- w/max(w)
    floor <- cfg$probability_floor_relative * mean(w)
    count <- sum(w < floor)
    w <- pmax(w, floor)
    list(w = w, floored = count, min_raw = raw_min, proxy = if (method == "hara_ensemble") "candidate_ensemble" else "locked_reference")
}

weighted_order <- function (pub, method, seed) 
{
    ww <- query_weights(pub, method)
    w <- ww$w
    set.seed(seed)
    o <- order(stats::rexp(pub$n, rate = w), seq_len(pub$n))
    remaining <- rev(cumsum(rev(w[o])))
    q <- w[o]/remaining
    need(!anyDuplicated(o) && all(is.finite(q)) && all(q > 0 & q <= 1 + TOL), "\u62bd\u6837\u6982\u7387\u5f02\u5e38")
    list(order = o, q = pmin(1, q), weight = w[o], remaining_weight = remaining, floored = ww$floored, min_raw = ww$min_raw, 
        proxy = ww$proxy)
}

fixed_order <- function (pub, method, seed) 
{
    if (method %in% LITERATURE_METHODS) 
        weighted_order(pub, method, seed)$order
    else fixed_order_131(pub, method, seed)
}

lure_weights <- function (n, M, q) 
{
    need(M >= 1 && M <= n && length(q) >= M && all(q[seq_len(M)] > 0), "LURE\u9884\u7b97/\u6982\u7387\u65e0\u6548")
    if (M == n) 
        return(rep(1, n))
    m <- seq_len(M)
    v <- 1 + (n - M)/(n - m) * (1/((n - m + 1) * q[m]) - 1)
    need(all(is.finite(v)) && all(v > 0), "LURE\u6743\u91cd\u975e\u6709\u9650\u6216\u975e\u6b63")
    v
}

fixed_budget_lure <- function (pub, reveal, draw, eps_values, budgets) 
{
    caps <- sort(unique(pmax(1L, as.integer(floor(pub$n * budgets)))))
    known <- numeric()
    rows <- decisions <- diagnostics <- list()
    paid <- 0L
    for (M in caps) {
        idx_new <- draw$order[seq.int(paid + 1L, M)]
        known <- c(known, reveal(idx_new))
        paid <- M
        ii <- draw$order[seq_len(M)]
        v <- lure_weights(pub$n, M, draw$q)
        loss <- (pub$P[ii, , drop = FALSE] - known)^2
        est_risk <- colSums(loss * v)/M
        gap <- sweep(loss, 1, loss[, pub$K], "-")
        est <- sweep(crossprod(gap * v, pub$G[ii, , drop = FALSE]), 2, pub$ng, "/") * pub$n/M
        d <- expand.grid(candidate = colnames(pub$P), group = colnames(pub$G), stringsAsFactors = FALSE)
        d$cap_labels <- M
        d$cap_fraction <- M/pub$n
        d$estimate <- as.vector(est)
        d$estimated_overall_risk <- rep(est_risk, ncol(pub$G))
        rows[[length(rows) + 1L]] <- d
        for (e in eps_values) {
            eligible <- apply(est <= e + TOL, 1, all)
            k <- if (any(eligible)) 
                which(eligible)[1L]
            else NA_integer_
            decisions[[length(decisions) + 1L]] <- data.frame(cap_labels = M, cap_fraction = M/pub$n, epsilon = e, 
                selected_index = k, selected_candidate = if (is.na(k)) 
                  NA_character_
                else colnames(pub$P)[k], selected_input_channels = if (is.na(k)) 
                  NA_real_
                else pub$costs[k], minimum_mean_loss_candidate = colnames(pub$P)[which.min(est_risk)], certificate = "NONE_POINT_ESTIMATE_ONLY")
        }
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(cap_labels = M, cap_fraction = M/pub$n, min_conditional_q = min(draw$q[seq_len(M)]), 
            max_lure_weight = max(v), sum_lure_weights = sum(v), weight_ess_heuristic = sum(v)^2/sum(v^2), floored_units = draw$floored, 
            proxy = draw$proxy)
    }
    list(estimates = bind(rows), decisions = bind(decisions), diagnostics = bind(diagnostics))
}

col_cs <- function (x) 
{
    x <- as.matrix(x)
    ans <- x
    for (j in seq_len(ncol(x))) ans[, j] <- cumsum(x[, j])
    ans
}

col_cmax <- function (x) 
{
    x <- as.matrix(x)
    for (j in seq_len(ncol(x))) x[, j] <- cummax(x[, j])
    x
}

col_cmin <- function (x) 
{
    x <- as.matrix(x)
    for (j in seq_len(ncol(x))) x[, j] <- cummin(x[, j])
    x
}

row_shift <- function (x, first) 
rbind(first, if (nrow(x) > 1L) x[-nrow(x), , drop = FALSE] else matrix(numeric(), 0, ncol(x)))

eb_state <- function (N, lower, upper, offset, alpha, J) 
{
    K <- length(lower)
    A <- length(alpha)
    rr <- upper - lower
    need(N >= 1 && all(is.finite(lower)) && all(rr >= 0), "CS\u533a\u57df\u652f\u6301\u8303\u56f4\u5f02\u5e38")
    list(N = N, t = 0L, lower = lower, upper = upper, range = rr, offset = offset, alpha = alpha, J = max(1L, J), 
        sumx = numeric(K), sse = numeric(K), wor_num = numeric(K), wor_den = 0, weighted_num = matrix(0, K, A), 
        weighted_den = matrix(0, K, A), penalty = matrix(0, K, A), L = matrix(lower, K, A), U = matrix(upper, K, 
            A))
}

eb_update <- function (s, X) 
{
    X <- as.matrix(X)
    b <- nrow(X)
    K <- ncol(X)
    need(b > 0 && K == length(s$lower) && s$t + b <= s$N, "CS\u66f4\u65b0\u7ef4\u5ea6\u6216\u533a\u57df\u8ba1\u6570\u5f02\u5e38")
    need(all(is.finite(X)) && all(sweep(X, 2, s$lower, ">=") | abs(sweep(X, 2, s$lower, "-")) < 1e-10) && all(sweep(X, 
        2, s$upper, "<=") | abs(sweep(X, 2, s$upper, "-")) < 1e-10), "\u67e5\u8be2\u635f\u5931\u8d85\u51fa\u516c\u5f00\u652f\u6301\u96c6")
    j <- s$t + seq_len(b)
    remain <- s$N - j + 1
    total <- sweep(col_cs(X), 2, s$sumx, "+")
    prior <- row_shift(total, s$sumx)
    means <- total/j
    vadd <- (X - means)^2
    ssepath <- sweep(col_cs(vadd), 2, s$sse, "+")
    prior_sse <- row_shift(ssepath, s$sse)
    corr <- X + prior/remain
    denstep <- s$N/remain
    wn <- sweep(col_cs(corr), 2, s$wor_num, "+")
    wd <- s$wor_den + cumsum(denstep)
    pwn <- row_shift(wn, s$wor_num)
    pwd <- c(s$wor_den, head(wd, -1L))
    center <- pwn/pmax(pwd, 1e-300)
    if (s$t == 0L) 
        center[1, ] <- (s$lower + s$upper)/2
    c0 <- pmax(s$range, 1e-14)
    predvar <- sweep(prior_sse, 2, c0^2/4, "+")/j
    for (a in seq_along(s$alpha)) {
        ell <- log(2 * s$J/s$alpha[a])
        lam <- sqrt(2 * ell/(pmax(predvar, 1e-300) * (j * log(j + 1))))
        lam <- pmin(lam, matrix(0.5/c0, b, K, byrow = TRUE))
        eta <- sweep(lam, 2, c0, "*")
        psi <- (-log1p(-eta) - eta)
        penalties <- sweep((X - center)^2 * psi, 2, c0^2, "/")
        nn <- sweep(col_cs(lam * corr), 2, s$weighted_num[, a], "+")
        dd <- sweep(col_cs(lam * denstep), 2, s$weighted_den[, a], "+")
        pp <- sweep(col_cs(penalties), 2, s$penalty[, a], "+")
        rawL <- (nn - pp - ell)/dd
        rawU <- (nn + pp + ell)/dd
        s$L[, a] <- pmax(s$L[, a], col_cmax(rawL)[b, ])
        s$U[, a] <- pmin(s$U[, a], col_cmin(rawU)[b, ])
        s$weighted_num[, a] <- nn[b, ]
        s$weighted_den[, a] <- dd[b, ]
        s$penalty[, a] <- pp[b, ]
    }
    s$t <- s$t + b
    s$sumx <- total[b, ]
    s$sse <- ssepath[b, ]
    s$wor_num <- wn[b, ]
    s$wor_den <- wd[b]
    constant <- s$range == 0
    if (any(constant)) {
        s$L[constant, ] <- s$lower[constant]
        s$U[constant, ] <- s$upper[constant]
    }
    if (s$t == s$N) {
        exact <- matrix(s$sumx/s$N, K, length(s$alpha))
        s$L <- pmax(s$L, exact)
        s$U <- pmin(s$U, exact)
    }
    s
}

cs_setup <- function (pub, method, alpha, seed) 
{
    need(method %in% CS_METHODS, "\u672a\u77e5CS\u67e5\u8be2\u7b56\u7565")
    residual <- method != "srs_wor_EB"
    if (method %in% c("stratified_residual_wor_EB", "decision_cell_residual_wor_EB")) {
        signature <- apply(pub$G, 1, paste0, collapse = "")
        cells <- split(seq_len(pub$n), factor(signature, levels = sort(unique(signature))))
        R <- matrix(0, pub$n, length(cells))
        for (h in seq_along(cells)) R[cells[[h]], h] <- 1
        W <- sweep(crossprod(R, pub$G), 2, pub$ng, "/")
        set.seed(seed)
        order <- stratified_order(signature)
    }
    else {
        R <- pub$G
        W <- diag(ncol(pub$G))
        set.seed(seed)
        order <- sample.int(pub$n)
    }
    offset <- if (residual) 
        pub$pr
    else matrix(0, pub$n, pub$K)
    regions <- lapply(seq_len(ncol(R)), function(h) {
        ii <- which(R[, h] == 1)
        lo <- pub$lo[ii, , drop = FALSE] - offset[ii, , drop = FALSE]
        hi <- pub$hi[ii, , drop = FALSE] - offset[ii, , drop = FALSE]
        list(ii = ii, lower = apply(lo, 2, min), upper = apply(hi, 2, max), offset = colMeans(offset[ii, , drop = FALSE]))
    })
    J <- sum(vapply(regions, function(x) sum(x$upper > x$lower), 0L))
    states <- lapply(regions, function(x) eb_state(length(x$ii), x$lower, x$upper, x$offset, alpha, J))
    list(R = R, W = W, states = states, offset = offset, order = order, J = J, metadata = data.frame(region = seq_along(regions), 
        n = vapply(regions, function(x) length(x$ii), 0L), nonconstant_constraints = vapply(regions, function(x) sum(x$upper > 
            x$lower), 0L)))
}

aggregate_cs <- function (setup, a, K) 
{
    L <- vapply(setup$states, function(s) s$L[, a] + s$offset, numeric(K))
    U <- vapply(setup$states, function(s) s$U[, a] + s$offset, numeric(K))
    list(L = matrix(L, K, length(setup$states)) %*% setup$W, U = matrix(U, K, length(setup$states)) %*% setup$W)
}

cs_path133 <- function (pub, reveal, method, seed, epsilon, alpha, budgets, batch, observer = NULL, keep_bounds = TRUE) 
{
    need(all(alpha > 0 & alpha < 1) && !anyDuplicated(alpha), "alpha\u5e94\u5728(0,1)\u4e14\u4e0d\u91cd\u590d")
    setup <- cs_setup(pub, method, alpha, seed)
    caps <- sort(unique(pmax(1L, as.integer(floor(pub$n * budgets)))))
    points <- sort(unique(c(if (max(caps) >= batch) seq.int(batch, max(caps), by = batch) else integer(), caps)))
    combinations <- expand.grid(epsilon = epsilon, alpha = alpha, stringsAsFactors = FALSE)
    combos <- lapply(seq_len(nrow(combinations)), function(i) list(epsilon = combinations$epsilon[i], ai = match(combinations$alpha[i], 
        alpha), done = FALSE, stopped = NA_integer_, dec = NULL, S = NULL))
    det <- initial_state(pub)
    disabled <- rep(FALSE, length(alpha))
    conflicts <- integer(length(alpha))
    seen <- integer()
    paid <- 0L
    snaps <- bounds <- list()
    used <- rep(FALSE, pub$n)
    adaptive <- method == "decision_cell_residual_wor_EB"
    cell_orders <- if (adaptive) 
        lapply(seq_len(ncol(setup$R)), function(h) setup$order[setup$R[setup$order, h] == 1])
    else NULL
    prefer <- which(combinations$epsilon == cfg$primary_epsilon & combinations$alpha == max(alpha))
    assess <- function(at) {
        SS <- vector("list", length(alpha))
        for (a in seq_along(alpha)) {
            raw <- aggregate_cs(setup, a, pub$K)
            L <- pmax(det$L, raw$L)
            U <- pmin(det$U, raw$U)
            if (!is.null(observer)) 
                observer(list(L = L, U = U), seen, alpha[a])
            if (any(L > U + 1e-10)) {
                disabled[a] <<- TRUE
                conflicts[a] <<- conflicts[a] + 1L
            }
            if (disabled[a]) {
                L <- det$L
                U <- det$U
            }
            else {
                tiny <- L > U
                if (any(tiny)) {
                  mid <- (L[tiny] + U[tiny])/2
                  L[tiny] <- U[tiny] <- mid
                }
            }
            SS[[a]] <- list(L = L, U = U, E = det$E)
        }
        for (i in seq_along(combos)) if (!combos[[i]]$done) {
            cc <- combos[[i]]
            cc$S <- SS[[cc$ai]]
            cc$dec <- decision(cc$S, cc$epsilon, pub$costs)
            if (cc$dec$certified) {
                cc$done <- TRUE
                cc$stopped <- at
            }
            combos[[i]] <<- cc
        }
    }
    assess(0L)
    for (point in points) {
        while (paid < point && !all(vapply(combos, `[[`, TRUE, "done"))) {
            if (adaptive) {
                active <- which(!vapply(combos, `[[`, TRUE, "done"))
                target <- if (length(intersect(prefer, active))) 
                  intersect(prefer, active)[1L]
                else active[1L]
                cc <- combos[[target]]
                score <- query_scores(pub, cc$S$L, cc$S$U, det$E, cc$epsilon, "decision_focused", cc$dec)
                remaining <- lapply(cell_orders, function(o) o[!used[o]])
                value <- vapply(remaining, function(ii) if (length(ii)) 
                  mean(score[ii])
                else -Inf, 0)
                h <- which.max(value)
                need(is.finite(value[h]), "\u81ea\u9002\u5e94\u5206\u5c42\u6ca1\u6709\u53ef\u67e5\u8be2\u8bb0\u5f55")
                ids <- head(remaining[[h]], point - paid)
            }
            else ids <- setup$order[seq.int(paid + 1L, point)]
            yy <- reveal(ids)
            need(length(yy) == length(ids) && all(yy %in% c(0, 1)), "CS oracle\u8fd4\u56de\u6807\u7b7e\u5f02\u5e38")
            seen <- c(seen, ids)
            used[ids] <- TRUE
            paid <- paid + length(ids)
            det <- update_state(det, pub, ids, yy)
            d <- pub$D0[ids, , drop = FALSE] + (pub$D1[ids, , drop = FALSE] - pub$D0[ids, , drop = FALSE]) * yy
            for (h in seq_along(setup$states)) {
                ii <- which(setup$R[ids, h] == 1)
                if (length(ii)) 
                  setup$states[[h]] <- eb_update(setup$states[[h]], d[ii, , drop = FALSE] - setup$offset[ids[ii], 
                    , drop = FALSE])
            }
            if (paid == pub$n) {
                det$L <- det$U <- (det$L + det$U)/2
                det$E <- det$L
            }
            assess(paid)
        }
        if (point %in% caps) 
            for (i in seq_along(combos)) {
                cc <- combos[[i]]
                dec <- cc$dec
                k <- dec$best
                nused <- if (cc$done) 
                  cc$stopped
                else paid
                snaps[[length(snaps) + 1L]] <- data.frame(cap_labels = point, cap_fraction = point/pub$n, labels_used = nused, 
                  fraction_used = nused/pub$n, epsilon = cc$epsilon, alpha = alpha[cc$ai], selected_index = k, 
                  selected_candidate = if (is.na(k)) 
                    NA_character_
                  else colnames(pub$P)[k], selected_input_channels = if (is.na(k)) 
                    NA_real_
                  else pub$costs[k], certified_minimum_cost = dec$certified, remaining_cheaper_candidates = length(dec$contenders), 
                  status = if (dec$certified) 
                    "STATISTICAL_MINIMUM_COST_CERTIFICATE"
                  else "STATISTICAL_COST_UNRESOLVED", guarantee = "finite-pool simultaneous anytime coverage >= 1-alpha", 
                  cs_conflict_fallback = disabled[cc$ai])
                if (keep_bounds && point == max(caps)) {
                  b <- expand.grid(candidate = colnames(pub$P), group = colnames(pub$G), stringsAsFactors = FALSE)
                  b$lower <- as.vector(cc$S$L)
                  b$upper <- as.vector(cc$S$U)
                  b$epsilon <- cc$epsilon
                  b$alpha <- alpha[cc$ai]
                  b$labels_used <- nused
                  bounds[[length(bounds) + 1L]] <- b
                }
            }
    }
    list(runs = bind(snaps), final_bounds = bind(bounds), queried_indices = seen, regions = setup$metadata, allocated_constraints = setup$J, 
        conflicts = conflicts)
}
