MVMMNWmix.EM = function(Y,
                        M,
                        L,
                        S = NA,
                        Ps = NA,
                        nu,
                        pp,
                        fix.nu = FALSE,
                        class = NA,
                        const = c("1,1", "det.1", "diag.n"),
                        M.const = c("I", "Col", "Row", "UN"),
                        L.const = c("I", "Col", "Row", "UN"),
                        row.variance = c("AR", "CS", "UN"),
                        col.variance = c("AR", "CS", "UN"),
                        rho.row = NA,
                        sigma.row = NA,
                        rho.col = NA,
                        sigma.col = NA,
                        fix.rho.row = FALSE,
                        fix.rho.col = FALSE,
                        fix.sigma.row = FALSE,
                        fix.sigma.col = FALSE,
                        tol = 1e-5,
                        Stop.rule = c("Log.like", "Aitken"),
                        max.iter = 1000,
                        per = 100,
                        print = TRUE)
{
  #' --------------------------------------------- ###
  #' --------- Required functions ---------------- ###
  #' --------------------------------------------- ###

  massage1 = combine_ansi_styles("red", "bold")
  massage2 = combine_ansi_styles("yellow", "bold")

  All.ones = function(p, n) matrix(1, nrow = p, ncol = n)

  Stop.Rule <- function(log.lik) {
    if (length(log.lik) >= 3) {
      n = length(log.lik)
      l.new = log.lik[n]
      l.old  = log.lik[(n - 1)]
      l.old2  = log.lik[(n - 2)]
      ait = (l.new - l.old) / (l.old - l.old2)
      ln.Inf = l.old + (l.new - l.old) / (1 - ait)
      out = ln.Inf - l.old
      if (!is.na(out))
        out = abs(out)
      else
        out = 0

    } else
      out = (log.lik[2] - log.lik[1]) / abs(log.lik[1])
    return(out)
  }

  tr = function(A) sum(diag(A))

  ARCSgen = function(n, rho, variance) {
    #' Generate AR(1) correlation matrices
    #' @param n number of columns/rows
    #' @param rho correlation parameter
    nbld = as.character(rho)
    if (n <= 1)
      stop(massage1("Error in AR-CS-generator: n must be greater than 1."), call. = FALSE)
    if (is.na(rho))
      stop(massage1("Error in AR-CS-generator: rho is not specified correctly."), call. = FALSE)
    if (rho >= 1)
      stop(massage1("Error in AR-CS-generator: rho must be a correlation less than 1."), call. = FALSE)
    if (rho <= -1)
      stop(massage1("Error in AR-CS-generator: rho must be a correlation greater than -1."), call. = FALSE)
    if (abs(rho) == 0)
      cli_alert_warning(massage2("Rho = {.emph {nbld}} and should be greater than 0."))
    if (rho > 0.999)
      cli_alert_warning(massage2("Rho = {.emph {nbld}}, high correlation may cause numerical problems."))

    #' Generate a AR(1) correlation matrix
    ARgen = function(n, rho) {
      X = toeplitz(c(1, rho ^ (1:(n - 1))))
      return(X)
    }
    #' Generate a compound symmetric correlation matrix
    CSgen = function(n, rho) {
      A = matrix(rho, nrow = n, ncol = n)
      diag(A) = 1
      A
    }

    if (variance == "AR(1)") return(ARgen(n, rho))
    if (variance == "CS") return(CSgen(n, rho))
  }

  varinv <- function(n, rho, deriv = F, variance) {
    if (!(n > 1))
      stop(massage1("n needs to be greater than 1."), call. = FALSE)
    if (!(rho < 1 && rho > -1))
      stop(massage1("rho needs to be < 1."), call. = FALSE)

    invAR <- function(n, rho, deriv = FALSE) {
      if (deriv) {
        X = toeplitz(c(4 * rho,-(rho ^ 2 + 1), rep(0, n - 2)))
        X[1, 1] <-  2 * rho
        X[n, n] <-  2 * rho
        return((1 / (1 - rho ^ 2) ^ 2) * X)
      }
      X = toeplitz(c(1 + rho ^ 2,-rho, rep(0, n - 2)))
      X[1, 1] <- X[n, n] <- 1
      return((1 / (1 - rho ^ 2)) * X)
    }

    invCS <- function(n, rho, deriv = FALSE) {
      alpha = sqrt(1 - rho)

      if (deriv) {
        X = diag(1 / alpha ^ 4, n) -
          matrix(((n - 1) * rho ^ 2 + 1) / ((1 - rho) ^ 2 * ((n - 1) * rho + 1) ^ 2),
                 nrow = n, ncol = n)
        return(X)
      }

      X = diag(1 / alpha ^ 2, n) -
        matrix(rho / (alpha ^ 2 * (alpha ^ 2 + n * rho)), nrow = n, ncol = n)
      return(X)
    }

    if (variance == "AR(1)")
      return(invAR(n, rho, deriv))
    if (variance == "CS")
      return(invCS(n, rho, deriv))
  }

  Cluster.error.rate = function (clust1, clust2) {
    clust1 <- unclass(as.ordered(clust1))
    clust2 <- unclass(as.ordered(clust2))
    if ((n = length(clust1)) != length(clust2)) {
      stop(massage1("Error in clustering evaluation: the length Groups are not equal."), call. = FALSE)
    }
    if ((g = length(table(clust1))) != length(table(clust2))) {
      stop(massage1("Error in clustering evaluation: the number of clusters are not equal."), call. = FALSE)
    }
    permute <- function(a) {
      n <- length(a)
      if (n == 1)
        f <- a
      else {
        nm <- gamma(n)
        f <- array(0, c(n, n * nm))
        j <- 1
        for (i in a) {
          f[1, (j - 1) * nm + 1:nm] <- i
          f[-1, (j - 1) * nm + 1:nm] <- permute(setdiff(a, i))
          j <- j + 1
        }
      }
      f
    }
    id <- 1:n
    cmb <- permute(1:g)
    nperm <- ncol(cmb)
    rate <- rep(0, nperm)
    for (i in 1:nperm) {
      tmp <- rep(0, g)
      tc <- rep(0, n)
      for (j in 1:g)
        tc[clust2 == j] = cmb[j, i]
      for (j in 1:g) {
        tmp1 <- 0
        for (k in (1:g)[-j])
          tmp1 <- tmp1 + length(intersect(id[clust1 ==
                                               j], id[tc == k]))
        tmp[j] <- tmp1
      }
      rate[i] <- sum(tmp) / n
    }
    min(rate)
  }

  ARI.fun = function (LabelA, LabelB) {
    u <- unclass(as.ordered(LabelA))
    v <- unclass(as.ordered(LabelB))
    if ((N <- length(u)) != length(v))
      stop(massage1("Error in ARI computing: Labels of the groups do not match!"), call. = FALSE)
    row <- max(u)
    col <- max(v)
    nvect <- array(0, c(row, col))
    for (i in 1:row) {
      for (j in 1:col) {
        nvect[i, j] <- sum(u == i & v == j)
      }
    }
    SumsA <- rowSums(nvect)
    SumsB <- colSums(nvect)
    a = 0
    for (i in 1:row)
      a = a + choose(SumsA[i], 2)
    b = 0
    for (j in 1:col)
      b = b + choose(SumsB[j], 2)
    c <- a * b / choose(N, 2)
    d = 0
    for (i in 1:row) {
      for (j in 1:col) {
        d = d + choose(nvect[i, j], 2)
      }
    }
    ind <- (d - c) / ((a + b) / 2 - c)
    ind
  }
  HF = function(x) exp(dnorm(x, log = T) - pnorm(x, log.p = T))

  f.MMNW = function(Y, M, L, S, Ps, nu, Log = F) {
    p = dim(Y)[1]
    n = dim(Y)[2]
    Ps.inv = solve(Ps)
    S.inv = solve(S)
    eta = tr(Ps.inv %*% t(L) %*% S.inv %*% L) + 2 * nu ^ 2
    f.int = function(Y) {
      A = (tr(Ps.inv %*% t(L) %*% S.inv %*% (Y - M))) / sqrt(eta)
      out = log(2 * nu ^ 2 / eta) + 0.5 * log(2 * pi) + 0.5 * A ^ 2 + pnorm(A, log.p = T) +
        log(A + HF(A)) - 0.5 * tr(Ps.inv %*% t(Y - M) %*% S.inv %*% (Y - M)) -
        0.5 * n * p * log(2 * pi) - 0.5 * p * log(det(Ps)) - 0.5 * n * log(det(S))
      return(out)
    }
    YY = lapply(seq(dim(Y)[3]), function(x)
      Y[, , x])
    Out = unlist(lapply(YY, function(x)
      f.int(x)))

    ifelse(Log == F, return(exp(Out)), return(Out))
  }

  f.MMNW.mix = function(Y, M, L, S, Ps, nu, PI, Log = F) {
    N = dim(Y)[3]
    Mat =  matrix(0, length(PI), N)
    for (i in 1:length(PI))
      Mat[i, ] = PI[i] * f.MMNW(
        Y,
        M = M[, , i],
        L = L[, , i],
        S = S[, , i],
        Ps = Ps[, , i],
        nu = nu[i]
      )
    Out =  colSums(Mat)
    Out[which(Out == 0)] <- .Machine$double.xmin
    ifelse(Log == F, return(Out), return(log(Out)))
  }

  #' --------------------------------------------- ###
  #' --------------------------------------------- ###

  begin = proc.time()[1]
  G = length(pp)
  p = dim(Y)[1]
  n = dim(Y)[2]
  N = dim(Y)[3]

  Rpar = p * (p + 1) / 2
  Cpar = n * (n + 1) / 2

  if (all(is.na(row.variance)) || (length(row.variance) > 1)) row.variance = "UN"
  if (all(is.na(col.variance)) || (length(col.variance) > 1)) col.variance = "UN"

  if (grepl("^ar", x = row.variance, ignore.case = TRUE)) {
    row.variance = "AR(1)"
    Rpar = 2
  }
  if (grepl("^ar", x = col.variance, ignore.case = TRUE)) {
    col.variance = "AR(1)"
    Cpar = 2
  }
  if (grepl("^cs", x = row.variance, ignore.case = TRUE)) {
    row.variance = "CS"
    Rpar = 2
  }
  if (grepl("^cs", x = col.variance, ignore.case = TRUE)) {
    col.variance = "CS"
    Cpar = 2
  }
  if (grepl("^un", x = row.variance, ignore.case = TRUE))
    row.variance = "UN"
  if (grepl("^un", x = col.variance, ignore.case = TRUE))
    col.variance = "UN"

  if(all(is.na(const)) || (length(const) > 1)) const = "1,1"
  if (!any(const == c("1,1", "det.1", "diag.n")))
    stop(massage1("const is not specified correctly."))

  if(all(is.na(M.const)) || (length(M.const) > 1)) M.const = "UN"
  if (!any(M.const == c("I", "Col", "Row", "UN")))
    stop(massage1("M.const is not specified correctly."))

  if(all(is.na(L.const)) || (length(L.const) > 1)) L.const = "UN"
  if (!any(L.const == c("I", "Col", "Row", "UN")))
    stop(massage1("L.const is not specified correctly."))

  if (!any(row.variance == c("AR(1)", "CS", "UN")))
    stop(massage1("row.variance is not specified correctly."))

  if (!any(col.variance == c("AR(1)", "CS", "UN")))
    stop(massage1("col.variance is not specified correctly."))

  if (row.variance != "UN" &&
      !all(is.numeric(rho.row)))
    stop(massage1("rho.row must be a numeric value."))
  if (col.variance != "UN" &&
      !all(is.numeric(rho.col)))
    stop(massage1("rho.col must be a numeric value."))
  if (row.variance != "UN" &&
      !all(is.numeric(sigma.row)))
    stop(massage1("sigma.row must be a numeric value."))
  if (col.variance != "UN" &&
      !all(is.numeric(sigma.col)))
    stop(massage1("sigma.col must be a numeric value."))
  if (!all(is.numeric(nu)))
    stop(massage1("nu must be a numeric value."))

  if (fix.nu)
    nu = rep(nu[1], G)
  if (fix.rho.row)
    rho.row = rep(rho.row[1], G)
  if (fix.rho.col)
    rho.col = rep(rho.col[1], G)
  if (fix.sigma.row)
    sigma.row = rep(sigma.row[1], G)
  if (fix.sigma.col)
    sigma.col = rep(sigma.col[1], G)

  if (row.variance != "UN")
    for (j in 1:G)
      S[, , j] = sigma.row[j] * ARCSgen(p, rho.row[j], variance = row.variance)
  if (col.variance != "UN") {
    for (j in 1:G) {
      Ps[, , j] = sigma.col[j] * ARCSgen(n, rho.col[j], variance = col.variance)
      if (const == "1,1")
        stop(massage1(
          "It is impossible to set const = 1,1 when col.variance is not UN."
        ))
      aa = ifelse(det(Ps[, , j]) < 0 ,-1, 1)
      if (const == "det.1")
        Ps[, , j] = aa * Ps[, , j] / abs(det(Ps[, , j])) ^ (1 / n)
      if (const == "diag.n")
        Ps[, , j] = n * Ps[, , j] / tr(Ps[, , j])
    }
  }

  lk = logli.old = sum(f.MMNW.mix(
    Y,
    M = M,
    L = L,
    S = S,
    Ps = Ps,
    nu = nu,
    PI = pp,
    Log = T
  ))

  iter = 0
  if (isTRUE(print)) {
    cat(rule(line = 2, line_col = "orange"), "\n")
    cat(rule(
      center = col_green(
        "* Fitting Mixture of matrix-variate MMNW distributions with ",
        row.variance,
        " structure for Sigma and ",
        col.variance,
        " structure for Psi",
        " *"
      ),
      line_col = "deeppink3",
      line = 2
    ),
    "\n")
    cat(rule(
      center = col_blue("Initial Log-likelihood = ", lk),
      line_col = "deeppink3",
      line = 2
    ),
    "\n")
    cat(rule(line = 1, line_col = "lightcoral"), "\n")
  }

  repeat {
    iter = iter + 1

    ### E step
    Zhat = PLOG = w = t = matrix(NA, nrow = N, ncol = G)
    for (i in 1:G) {
      PLOG[, i] = log(pp[i]) + f.MMNW(
        Y,
        M = M[, , i],
        L = L[, , i],
        S = S[, , i],
        Ps = Ps[, , i],
        nu = nu[i],
        Log = T
      )
    }

    for (i in 1:G) {
      Zhat[, i] = 1 / rowSums(exp(PLOG - PLOG[, i]))
    }

    for (i in 1:G) {
      Ps.inv = solve(Ps[, , i])
      S.inv = solve(S[, , i])
      eta = tr(Ps.inv %*% t(L[, , i]) %*% S.inv %*% L[, , i]) + 2 * nu[i] ^
        2
      A = c()
      for (j in 1:N) {
        A[j] = tr(Ps.inv %*% t(L[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) /
          sqrt(eta)
      }
      SS = A / sqrt(eta) + HF(A) / sqrt(eta)
      w[, i] = sqrt(eta) * (A * SS / sqrt(eta) + 1 / eta) / (HF(A) + A)
      t[, i] = sqrt(eta) * (A * (A / sqrt(eta) * SS + 1 / eta) / sqrt(eta) +
                              2 * SS / eta) / (HF(A) + A)
    }

    TZ = t * Zhat
    WZ = w * Zhat

    STZ = colSums(TZ)
    SWZ = colSums(WZ)
    #### CM step ####
    ### step 1
    Ni = colSums(Zhat)
    pp = Ni / N

    nu = sqrt(Ni / STZ)
    if (fix.nu)
      nu = rep(nu[1], G)

    #### CM step ####
    ### step 1

    for (i in 1:G) {

      wbar = SWZ[i] / Ni[i]
      tbar = STZ[i] / Ni[i]

      if (M.const == "UN") {
        Ybar = matrix(0, p, n)
        for (j in 1:N) {
          Ybar = Ybar + Zhat[j, i] * Y[, , j] / Ni[i]
        }
        M[, , i] = Ybar - L[, , i] * wbar
      }

      if (M.const == "I") {
        Ps.inv = solve(Ps[, , i])
        S.inv = solve(S[, , i])
        SumY = 0
        for (j in 1:N) {
          SumY = SumY + Zhat[j, i] *
            tr(S.inv %*% (Y[, , j] - w[j, i] * L[, , i]) %*% Ps.inv %*% t(All.ones(p, n)))
        }
        #sum(diag(A %*% All.ones(p, n) %*% B %*% t(All.ones(p, n)))) == sum(kronecker(A, B))
        M[, , i] = (SumY / Ni[i] / sum(kronecker(Ps.inv, S.inv))) * All.ones(p, n)
      }

      if (M.const == "Col") {
        Ps.inv = solve(Ps[, , i])
        SumY = matrix(0, p, n)
        for (j in 1:N) {
          SumY = SumY + Zhat[j, i] * ((Y[, , j] - w[j, i] * L[, , i]) %*% Ps.inv)
        }
        # All.ones(1, p) %*% A %*% All.ones(p, 1) == sum(A)
        M[, , i] = (SumY / Ni[i] / sum(Ps.inv)) %*% All.ones(n, n)
      }

      if (M.const == "Row") {
        S.inv = solve(S[, , i])
        SumY = matrix(0, p, n)
        for (j in 1:N) {
          SumY = SumY + Zhat[j, i] * (S.inv %*% (Y[, , j] - w[j, i] * L[, , i]))
        }
        M[, , i] = All.ones(p, p) %*% (SumY / Ni[i] / sum(S.inv))
      }

      if (L.const == "UN") {
        Ybar = matrix(0, p, n)
        Ywbar = matrix(0, p, n)
        for (j in 1:N) {
          Ywbar = Ywbar + WZ[j, i] * Y[, , j] / Ni[i]
          Ybar = Ybar + Zhat[j, i] * Y[, , j] / Ni[i]
        }
        L[, , i] = (Ywbar - Ybar * wbar) / (tbar - wbar ^ 2)
      }

      if (L.const == "I") {
        Ps.inv = solve(Ps[, , i])
        S.inv = solve(S[, , i])
        SumY = 0
        for (j in 1:N) {
          SumY = SumY + WZ[j, i] *
            tr(S.inv %*% (Y[, , j] - M[, , i]) %*% Ps.inv %*% t(All.ones(p, n)))
        }
        L[, , i] = (SumY / STZ[i] / sum(kronecker(Ps.inv, S.inv))) * All.ones(p, n)
      }

      if (L.const == "Col") {
        Ps.inv = solve(Ps[, , i])
        SumY = matrix(0, p, n)
        for (j in 1:N) {
          SumY = SumY + WZ[j, i] * ((Y[, , j] - M[, , i]) %*% Ps.inv)
        }
        # All.ones(1, p) %*% A %*% All.ones(p, 1) == sum(A)
        L[, , i] = (SumY / STZ[i] / sum(Ps.inv)) %*% All.ones(n, n)
      }

      if (L.const == "Row") {
        S.inv = solve(S[, , i])
        SumY = matrix(0, p, n)
        for (j in 1:N) {
          SumY = SumY + WZ[j, i] * (S.inv %*% (Y[, , j] - M[, , i]))
        }
        L[, , i] = All.ones(p, p) %*% (SumY / STZ[i] / sum(S.inv))
      }

      ### step 2
      if (row.variance == "UN") {
        int = matrix(0, p, p)
        Ps.inv = solve(Ps[, , i])

        for (j in 1:N) {
          int = int + (
            Zhat[j, i] * (Y[, , j] - M[, , i]) %*% Ps.inv %*% t(Y[, , j] - M[, , i]) +
              TZ[j, i] * L[, , i] %*% Ps.inv %*% t(L[, , i]) -
              WZ[j, i] * ((L[, , i]) %*% Ps.inv %*% t(Y[, , j] - M[, , i]) +
                            (Y[, , j] - M[, , i]) %*% Ps.inv %*% t(L[, , i])
              )
          )
        }
        S[, , i] = int / Ni[i] / n
        S.inv = solve(S[, , i])
      }

      ### step 3
      if (col.variance == "UN") {
        int = matrix(0, n, n)
        for (j in 1:N) {
          int = int + (
            Zhat[j, i] * t(Y[, , j] - M[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i]) +
              TZ[j, i] * t(L[, , i]) %*% S.inv %*% L[, , i] -
              WZ[j, i] * (t(L[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i]) +
                            t(Y[, , j] - M[, , i]) %*% S.inv %*% L[, , i])
          )
        }
        Ps[, , i] = int / Ni[i] / p


        if (const == "1,1")
          Ps[, , i] = Ps[, , i] / Ps[1, 1, i]
        aa = ifelse(det(Ps[, , i]) < 0 ,-1, 1)
        if (const == "det.1")
          Ps[, , i] = aa * Ps[, , i] / abs(det(Ps[, , i])) ^ (1 / n)
        if (const == "diag.n")
          Ps[, , i] = n * Ps[, , i] / tr(Ps[, , i])
      }

      if (fix.rho.row == F) {
        if (fix.sigma.row == F) {
          if (row.variance == "AR(1)" || row.variance == "CS") {
            ROW.RHO = ARCSgen(p, rho.row[i], variance = row.variance)
            SS = 0
            Ps.inv = solve(Ps[, , i])
            ROW.RHO.inv = varinv(p,
                                 rho.row[i],
                                 deriv = F,
                                 variance = row.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                    tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (L[, , i]))
                )
            }
            sigma.row[i] = SS / (n * p) / Ni[i]

            F.Rho1 = function(rho) {
              S = ARCSgen(p, rho, variance = row.variance)
              - sum(Zhat[, i] * f.MMNW(
                Y,
                M = M[, , i],
                L = L[, , i],
                S = sigma.row[i] * S,
                Ps = Ps[, , i],
                nu = nu[i],
                Log = T
              ))
            }
            Estim = try(nlminb(rho.row[i],
                               F.Rho1,
                               lower = -0.999,
                               upper = 0.999))
            if(!('try-error' %in% class(Estim))){
              rho.row[i] = Estim$par
            } else {rho.row[i] = rho.row[i]}

            S[, , i] = sigma.row[i] * ARCSgen(p, rho.row[i], variance = row.variance)
          }
        }
      }

      if (fix.rho.col == F) {
        if (fix.sigma.col == F) {
          if (col.variance == "AR(1)" || col.variance == "CS") {
            COL.RHO = ARCSgen(n, rho.col[i], variance = col.variance)
            SS = 0
            S.inv = solve(S[, , i])
            COL.RHO.inv = varinv(n,
                                 rho.col[i],
                                 deriv = F,
                                 variance = col.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                    tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (L[, , i]))
                )
            }
            sigma.col[i] = SS / (n * p) / Ni[i]

            F.Rho2 = function(rho) {
              Ps = ARCSgen(n, rho, variance = col.variance)
              - sum(Zhat[, i] * f.MMNW(
                Y,
                M = M[, , i],
                L = L[, , i],
                S = S[, , i],
                Ps = sigma.col[i] * Ps,
                nu = nu[i],
                Log = T
              ))
            }
            Estim = try(nlminb(rho.col[i],
                               F.Rho2,
                               lower = -0.999,
                               upper = 0.999))
            if(!('try-error' %in% class(Estim))){
              rho.col[i] = Estim$par
            } else {rho.col[i] = rho.col[i]}

            Ps[, , i] = sigma.col[i] * ARCSgen(n, rho.col[i], variance = col.variance)
            aa = ifelse(det(Ps[, , i]) < 0 ,-1, 1)
            if (const == "det.1")
              Ps[, , i] = aa * Ps[, , i] / abs(det(Ps[, , i])) ^ (1 / n)
            if (const == "diag.n")
              Ps[, , i] = n * Ps[, , i] / tr(Ps[, , i])
          }
        }
      }
    }

    if (fix.rho.row) {
      if (fix.sigma.row == F) {
        if (row.variance == "AR(1)" || row.variance == "CS") {
          for (i in 1:G) {
            SS = 0
            Ps.inv = solve(Ps[, , i])
            ROW.RHO.inv = varinv(p,
                                 rho.row[1],
                                 deriv = F,
                                 variance = row.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                    tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (L[, , i]))
                )
            }
            sigma.row[i] = SS / (n * p) / Ni[i]
          }
          F.Rho3 = function(Rrho) {
            for (i in 1:G)
              S[, , i] = sigma.row[i] * ARCSgen(p, Rrho, variance = row.variance)
            - sum(f.MMNW.mix(
              Y,
              M = M,
              L = L,
              S = S,
              Ps = Ps,
              nu = nu,
              PI = pp,
              Log = T
            ))
          }
          Estim = try(nlminb(rho.row[1], F.Rho3, lower = -0.999, upper = 0.999))
          if(!('try-error' %in% class(Estim))){
            rho.row = Estim$par
          } else { rho.row = rho.row }

          for (i in 1:G)
            S[, , i] = sigma.row[i] * ARCSgen(p, rho.row, variance = row.variance)

        }
      }

      if (fix.sigma.row) {
        if (row.variance == "AR(1)" || row.variance == "CS") {
          SS = 0
          for (i in 1:G) {
            Ps.inv = solve(Ps[, , i])
            ROW.RHO.inv = varinv(p,
                                 rho.row[1],
                                 deriv = F,
                                 variance = row.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                    tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (L[, , i]))
                )
            }
          }
          sigma.row = SS / (n * p) / N / G
          F.Rho4 = function(Rrho) {
            for (i in 1:G)
              S[, , i] = sigma.row[1] * ARCSgen(p, Rrho, variance = row.variance)
            - sum(f.MMNW.mix(
              Y,
              M = M,
              L = L,
              S = S,
              Ps = Ps,
              nu = nu,
              PI = pp,
              Log = T
            ))
          }
          Estim = try(nlminb(rho.row[1], F.Rho4, lower = -0.999, upper = 0.999))
          if(!('try-error' %in% class(Estim))){
            rho.row = Estim$par
          } else { rho.row = rho.row }

          for (i in 1:G)
            S[, , i] = sigma.row[1] * ARCSgen(p, rho.row, variance = row.variance)
        }
      }
    }

    if (fix.rho.row == F) {
      if (fix.sigma.row) {
        if (row.variance == "AR(1)" || row.variance == "CS") {
          SS = 0
          for (i in 1:G) {
            Ps.inv = solve(Ps[, , i])
            ROW.RHO.inv = varinv(p,
                                 rho.row[1],
                                 deriv = F,
                                 variance = row.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(Ps.inv %*% t(L[, , i]) %*% ROW.RHO.inv %*% (Y[, , j] - M[, , i])) +
                    tr(Ps.inv %*% t(Y[, , j] - M[, , i]) %*% ROW.RHO.inv %*% (L[, , i]))
                )
            }
          }
          sigma.row = SS / (n * p) / N / G
          for (i in 1:G) {
            F.Rho5 = function(rho) {
              S = ARCSgen(p, rho, variance = row.variance)
              - sum(Zhat[, i] * f.MMNW(
                Y,
                M = M[, , i],
                L = L[, , i],
                S = sigma.row * S,
                Ps = Ps[, , i],
                nu = nu[i],
                Log = T
              ))
            }
            Estim = try(nlminb(rho.row[i],
                               F.Rho5,
                               lower = -0.999,
                               upper = 0.999))
            if(!('try-error' %in% class(Estim))){
              rho.row[i] = Estim$par
            } else { rho.row[i] = rho.row[i] }
          }
          for (i in 1:G)
            S[, , i] = sigma.row[1] * ARCSgen(p, rho.row[i], variance = row.variance)

        }
      }
    }

    if (fix.rho.col) {
      if (fix.sigma.col == F) {
        if (col.variance == "AR(1)" || col.variance == "CS") {
          for (i in 1:G) {
            COL.RHO = ARCSgen(n, rho.col[1], variance = col.variance)
            SS = 0
            S.inv = solve(S[, , i])
            COL.RHO.inv = varinv(n,
                                 rho.col[1],
                                 deriv = F,
                                 variance = col.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                    tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (L[, , i]))
                )
            }
            sigma.col[i] = SS / (n * p) / Ni[i]
          }
          F.Crho1 = function(Crho) {
            for (i in 1:G)
              Ps[, , i] = sigma.col[i] * ARCSgen(n, Crho, variance = col.variance)
            - sum(f.MMNW.mix(
              Y,
              M = M,
              L = L,
              S = S,
              Ps = Ps,
              nu = nu,
              PI = pp,
              Log = T
            ))
          }
          Estim = try(nlminb(rho.col[1],
                             F.Crho1,
                             lower = -0.999,
                             upper = 0.999))
          if(!('try-error' %in% class(Estim))){
            rho.col = Estim$par
          } else { rho.col = rho.col }

          for (i in 1:G) {
            Ps[, , i] = sigma.col[i] * ARCSgen(n, rho.col, variance = col.variance)
            aa = ifelse(det(Ps[, , i]) < 0 ,-1, 1)
            if (const == "det.1")
              Ps[, , i] = aa * Ps[, , i] / abs(det(Ps[, , i])) ^ (1 / n)
            if (const == "diag.n")
              Ps[, , i] = n * Ps[, , i] / tr(Ps[, , i])
          }
        }
      }

      if (fix.sigma.col) {
        if (col.variance == "AR(1)" || col.variance == "CS") {
          SS = 0
          for (i in 1:G) {
            COL.RHO = ARCSgen(n, rho.col[1], variance = col.variance)
            S.inv = solve(S[, , i])
            COL.RHO.inv = varinv(n,
                                 rho.col[1],
                                 deriv = F,
                                 variance = col.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                    tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (L[, , i]))
                )
            }
          }
          sigma.col = SS / (n * p) / Ni / G
          F.Crho2 = function(Crho) {
            for (i in 1:G)
              Ps[, , i] = sigma.col[1] * ARCSgen(n, Crho, variance = col.variance)
            - sum(f.MMNW.mix(
              Y,
              M = M,
              L = L,
              S = S,
              Ps = Ps,
              nu = nu,
              PI = pp,
              Log = T
            ))
          }
          Estim = try(nlminb(rho.col[1],
                             F.Crho2,
                             lower = -0.999,
                             upper = 0.999))
          if(!('try-error' %in% class(Estim))){
            rho.col = Estim$par
          } else { rho.col = rho.col }

          for (i in 1:G){
            Ps[, , i] = sigma.col[1] * ARCSgen(n, rho.col, variance = col.variance)

            aa = ifelse(det(Ps[, , i]) < 0 , -1, 1)
            if (const == "det.1")
              Ps[, , i] = aa * Ps[, , i] / abs(det(Ps[, , i])) ^ (1 / n)
            if (const == "diag.n")
              Ps[, , i] = n * Ps[, , i] / tr(Ps[, , i])
          }
        }
      }
    }

    if (fix.rho.col == F) {
      if (fix.sigma.col) {
        if (col.variance == "AR(1)" || col.variance == "CS") {
          SS = 0
          for (i in 1:G) {
            COL.RHO = ARCSgen(n, rho.col[1], variance = col.variance)
            S.inv = solve(S[, , i])
            COL.RHO.inv = varinv(n,
                                 rho.col[1],
                                 deriv = F,
                                 variance = col.variance)
            for (j in 1:N) {
              SS = SS + Zhat[j, i] * tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                TZ[j, i] * tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (L[, , i])) -
                WZ[j, i] * (
                  tr(COL.RHO.inv %*% t(L[, , i]) %*% S.inv %*% (Y[, , j] - M[, , i])) +
                    tr(COL.RHO.inv %*% t(Y[, , j] - M[, , i]) %*% S.inv %*% (L[, , i]))
                )
            }
          }
          sigma.col = SS / (n * p) / Ni / G

          for (i in 1:G) {
            F.Rho6 = function(rho) {
              Ps = ARCSgen(n, rho, variance = col.variance)
              - sum(Zhat[, i] * f.MMNW(
                Y,
                M = M[, , i],
                L = L[, , i],
                S = S[, , i],
                Ps = sigma.col * Ps,
                nu = nu[i],
                Log = T
              ))
            }
            Estim = try(nlminb(rho.col[i],
                               F.Rho6,
                               lower = -0.999,
                               upper = 0.999))
            if(!('try-error' %in% class(Estim))){
              rho.col[i] = Estim$par
            } else {rho.col[i] = rho.col[i]}

            Ps[, , i] = sigma.col * ARCSgen(n, rho.col[i], variance = col.variance)
            aa = ifelse(det(Ps[, , i]) < 0 ,-1, 1)
            if (const == "det.1")
              Ps[, , i] = aa * Ps[, , i] / abs(det(Ps[, , i])) ^ (1 / n)
            if (const == "diag.n")
              Ps[, , i] = n * Ps[, , i] / tr(Ps[, , i])
          }
        }
      }
    }

    logli.new = sum(f.MMNW.mix(
      Y,
      M = M,
      L = L,
      S = S,
      Ps = Ps,
      nu = nu,
      PI = pp,
      Log = T
    ))

    if(is.infinite(logli.new)) logli.new = logli.old

    # Aikten's method
    lk = c(lk, logli.new)
    if (Stop.rule == "Log.like")
      epsilon = (logli.new - logli.old) / abs(logli.old)
    else {
      epsilon = Stop.Rule(lk)
    }

    diff = logli.new - logli.old
    if (iter %% per == 0 | is.na(diff))
    {
      if (isTRUE(print)) {
        cat(rule(
          left = col_cyan(
            'iter = ',
            iter,
            ', logli = ',
            logli.new,
            ", log.Like's diff = ",
            diff,
            ", Criteria.'s diff = ",
            epsilon
          ),
          line_col = "deeppink3",
          line = 1
        ),
        "\n")
        cat(rule(line = 1, line_col = "lightcoral"), "\n")
      }
    }

    if (is.na(diff)) {
      logli.new = logli.old
      break
    }
    if (is.na(epsilon))
    {
      logli.new = logli.old
      break
    }
    if (iter > 1)
      if (epsilon < tol | iter == max.iter)
        break
    logli.old = logli.new
  }

  if (!is.na(epsilon)) {
    for (i in 1:G) {
      PLOG[, i] = log(pp[i]) + f.MMNW(
        Y,
        M = M[, , i],
        L = L[, , i],
        S = S[, , i],
        Ps = Ps[, , i],
        nu = nu[i],
        Log = T
      )
    }

    for (i in 1:G) {
      Zhat[, i] = 1 / rowSums(exp(PLOG - PLOG[, i]))
    }
  }

  Z = Zhat

  Mpar = ifelse(M.const == "I", 1, ifelse(M.const == "Col", p, ifelse(M.const == "Row", n, p * n )))
  Lpar = ifelse(L.const == "I", 1, ifelse(L.const == "Col", p, ifelse(L.const == "Row", n, p * n )))

  m = G * (2 + Mpar + Lpar + Rpar + Cpar) - 1
  if (fix.nu)
    m = m - G
  if (fix.sigma.row)
    m = m - G + 1
  if (fix.sigma.col)
    m = m - G + 1
  if (fix.rho.row)
    m = m - G + 1
  if (fix.rho.col)
    m = m - G + 1

  BIC = -2 * logli.new + m * log(N)
  AIC = 2 * (-logli.new + m)

  ent = -sum(Z  * log(Z))
  ICL = BIC + 2 * ent
  CLC = 2 * (-logli.new + ent)
  EDC = -2 * logli.new + 0.2 * m * sqrt(N)

  Group = apply(Z, 1, which.max)
  di = table(Group)

  end = proc.time()[1]
  time = end - begin
  if (isTRUE(print)) {
    cat(rule(
      center = col_green('** Convergence ', epsilon < tol, ' **'),
      line_col = "deeppink3",
      line = 1
    ),
    "\n")
    cat(rule(
      center = col_green("iter = ", iter, " logli = ",  logli.new),
      line_col = "deeppink3",
      line = 1
    ),
    "\n")
    cat(rule(
      center = col_green(
        '** The process of fitting MIX-MV-MMNW distributions with the specified ',
        'structures for Sigma and Psi takes ',
        time,
        ' seconds **'
      ),
      line_col = "deeppink3",
      line = 1
    ),
    "\n")
  }

  if (is.numeric(class) == T) {
    true.clus = class
    km.clus = Group

    tab = table(true.clus, km.clus)
    MCR = Cluster.error.rate(km.clus, true.clus)
    ARI = ARI.fun(km.clus, true.clus)
    #variation.information = mcclust::vi.dist(true.clus, km.clus)

    Out = list(
      para.len = m,
      ni = di,
      Group = Group,
      p = pp,
      M = M,
      Lambda = L,
      Sigma = S,
      Psi = Ps,
      nu = nu,
      sigma.row = sigma.row,
      sigma.col = sigma.col,
      rho.row = rho.row,
      rho.col = rho.col,
      logli = logli.new,
      AIC = AIC,
      BIC = BIC,
      ICL = ICL,
      EDC = EDC,
      CLC = CLC,
      MCR = MCR,
      #variation.information = variation.information,
      Cross.Tab = tab,
      ARI = ARI,
      time = time
    )
  }
  else{
    Out = list(
      para.len = m,
      ni = di,
      Group = Group,
      p = pp,
      M = M,
      Lambda = L,
      Sigma = S,
      Psi = Ps,
      nu = nu,
      sigma.row = sigma.row,
      sigma.col = sigma.col,
      rho.row = rho.row,
      rho.col = rho.col,
      logli = logli.new,
      AIC = AIC,
      BIC = BIC,
      ICL = ICL,
      EDC = EDC,
      CLC = CLC,
      time = time
    )
  }
  return(Out)
}
