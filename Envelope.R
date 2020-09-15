gwer.envelope <- function (object, B = 100, arg, ident = NULL,...) 
{
  initial <- NULL ; data = object$lm$data
  X <- model.matrix(object$lm$terms)
  Xd <- as.matrix(object$lm$Xmodel)
  n <- nrow(Xd)
  p <- ncol(Xd)
  W <- object$gweights
  l.fitted = object$flm
  family <- object$family
  control <- object$lm$control
  #  ro <- object$resid
  #  tdf <- ro/sqrt(object$scalevariance)
  resid <- (object$y - object$fitted)/sqrt(object$dispersion)
  H <- matrix(0,n,n)
  for(i in 1:n)
    H[i,] <- Xd[i,] %*% solve(t(Xd) %*% diag(W[i,]) %*% Xd) %*% t(Xd) %*% diag(W[i,])  
  
  scale <- 4 * family$g2(resid, df = family$df, 
                         r = family$r, s = family$s, alpha = family$alpha, 
                         mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                         k = family$k) 
  scalevariance <- family$g4(resid, df = family$df, 
                             r = family$r, s = family$s, alpha = family$alpha, 
                             mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                             k = family$k)
  H1 <- (1/(scalevariance * scale)) * H
  varr <- scalevariance * object$dispersion * (diag(1, n) - H1)
  varr <- diag(varr)
  ri <- object$lm$y - object$fitted
  tdf <- ri/sqrt(varr)
  e <- e.i <- matrix(0, n, B) ; med <- matrix(0, n, n)
  resp <- NULL
  method <- "gwer.fit.envel"
  elliptical.fitter <- get(method)
  #  offset = object$offset
  #  if (length(offset) == 1 && offset == 0) 
  #    offset <- rep(0, nobs)
  
  for(i in 1:n){
    mu <- object$flm[i,]
    phi <- object$dispersion[i]
    w.i <- W[i,]
    for (j in 1:B) {
      dist <- object$family[[1]]
      if (charmatch(dist, "Normal", F)) {
        resp <- rnorm(n, 0, 1)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Normal(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Cauchy", F)) {
        resp <- rcauchy(n, 0, 1)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Cauchy(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Student", F)) {
        resp <- rt(n, arg)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Student(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Gstudent", F)) {
        resp <- rgstudent(n, arg[1], arg[2])
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Gstudent(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "LogisI", F)) {
        stop(paste("not implemented yet"))
        resp <- rlogisI(n, 0, 1)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = LogisI(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "LogisII", F)) {
        resp <- rlogisII(n)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = LogisII(), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Glogis", F)) {
        stop(paste("not implement yet"))
        resp <- rglogis(n, arg[1], arg[2])
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Glogis(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Cnormal", F)) {
        stop(paste("not implemented yet"))
        resp <- rcnormal(n, arg[1], arg[2])
        #fit <- elliptical(resp ~ X + (-1), family = Cnormal(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      else if (charmatch(dist, "Powerexp", F)) {
        resp <- rpowerexp(n, arg)
        resp <- mu + sqrt(phi) * resp
        #fit <- elliptical(resp ~ X + (-1), family = Powerexp(arg), 
        #                  control = glm.control(maxit = 1000))
      }
      
      lm.i <- elliptical.fitter(X = X, Y = resp, gweights = w.i, family = family, offset = NULL,
                                dispersion = NULL, maxit = control$maxit, epsilon = control$epsilon,
                                trace = control$trace, ...)
      ro.i <- lm.i$resid
      td.i <- ro.i/sqrt(lm.i$scalevariance)
      Xd.i <- as.matrix(lm.i$Xmodel)
      H.i <- Xd.i %*% solve(t(Xd.i) %*% diag(w.i) %*% Xd.i) %*% t(Xd.i) %*% diag(w.i) 
      H1.i <- (1/(lm.i$scalevariance * lm.i$scale)) * H.i
      varr.i <- lm.i$scalevariance * lm.i$dispersion * (diag(1, n) - H1.i)
      varr.i <- diag(varr.i)
      ri.i <- resp - lm.i$fitted
      td.i <- ri.i/sqrt(varr.i)
      e.i[, j] <- td.i #sort(td.i)
    }
    e[i, ] <- e.i[i, ]
  }
  e <- apply(e, 2, sort)
  
  e1 <- numeric(n)
  e2 <- numeric(n)
  e3 <- numeric(n)
  e4 <- numeric(n)
  e5 <- numeric(n)
  e6 <- numeric(n)
  e7 <- numeric(n)
  for (i in 1:n) {
    eo <- sort(e[i, ])
    e1[i] <- eo[ceiling(B * 0.05)]
    e2[i] <- eo[ceiling(B * 0.95)]
  }
  e3 <- t(t(apply(e, 2, mean)))
  e4 <- t(t(apply(e, 2, vari)))
  e5 <- t(t(apply(e, 2, skewn)))
  e6 <- t(t(apply(e, 2, kurt)))
  e7 <- cbind(e3, e4, e5, e6)
  desc <- apply(e, 2, mean)
  med <- apply(e, 1, mean)
  faixa <- range(tdf, e1, e2)
  screen(4)
  par(pty = "s")
  points.p <- qqnorm(tdf, xlab = "Quantiles of N(0,1)", ylab = "Standardized residual", 
                     ylim = faixa, pch = 16, main = '')
  par(new = TRUE)
  qqnorm(e1, axes = F, xlab = "", ylab = "", type = "l", ylim = faixa, 
         lty = 1, main = '')
  par(new = TRUE)
  qqnorm(e2, axes = F, xlab = "", ylab = "", type = "l", ylim = faixa, 
         lty = 1, main = '')
  par(new = TRUE)
  qqnorm(med, axes = F, xlab = "", ylab = "", type = "l", 
         ylim = faixa, lty = 2, main = '')
  if(!is.null(ident))
    identify(points.p$x, points.p$y, n = ident)
  x <- list(mean = desc[1], var = desc[2], skewness = desc[3], 
            kurtosis = desc[4], tdf = tdf)
  invisible(x)
}



gwer.fit.envel <- function (X, Y, gweights=NULL, offset, family, dispersion, 
                      maxit, epsilon, trace) 
{
  n <- nrow(X)
  if(is.null(gweights))
    gweights <- rep(1,n)
  if (is.null(offset)) 
    offset <- rep(0, n)
  
  p <- ncol(X)
  aux.model <- glm.fit(x = X, y = Y, offset = offset, weights = gweights,
                       family = gaussian())
  attr(aux.model, "class") <- c("glm", "lm")
  start <- aux.model$coef
  
  is.null.disp <- is.null(dispersion)
  elliptical.disp <- !is.null.disp && !is.numeric(dispersion)
  if (is.null.disp){
    options(warn = -1)
    dispersion <- (summary(aux.model)$dispersion)
    options(warn = 0) 
  }
  if (elliptical.disp){
    options(warn = -1)
    dispersion <- (summary(aux.model)$dispersion)
    options(warn = 0)
  }
  args <- resid(aux.model)/sqrt(dispersion)
  
  if (any(nas <- is.na(start))) {
    names(nas) <- dimnames(X)[[2]]
    X <- X[, !nas]
    aux.model <- glm.fit(x = X, y = Y, offset = offset, weights = gweights,
                         family = gaussian())
    attr(aux.model, "class") <- c("glm", "lm")
    start <- aux.model$coef
    options(warn = -1)
    dispersion <- (summary(aux.model)$dispersion)
    options(warn = 0)
  }
  
  
  linearconv <- TRUE
  iter <- 1
  error2 <- error3 <- 0
  repeat {
    if (trace) 
      cat("\n iteration", iter, ":")
    {
      w.1 <- family$g1(args, df = family$df, r = family$r, 
                       s = family$s, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, 
                       k = family$k)
      dg <- family$g2(args, df = family$df, r = family$r, 
                      s = family$s, alpha = family$alpha, mp = family$mp, 
                      epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      fg <- family$g3(args, df = family$df, r = family$r, 
                      s = family$s, alpha = family$alpha, mp = family$mp, 
                      epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      
      y.aux <- Y - offset
      w.h <- as.vector(-2 * w.1 * gweights)
      
      aux.model <- glm.fit(x = X, y = y.aux, weights = w.h, 
                           family = gaussian())
      attr(aux.model, "class") <- c("glm", "lm")
      new.start <- coef(aux.model)
      }
    error1 <- max(abs((new.start - start)/start))
    start <- new.start
    
    abs.res <- Y - X %*% start - offset
    
    if (is.null.disp) {
      aux.dispersion <- dispersion
      new.dispersion <- sum((-2 * w.1 * gweights) * abs.res^2)/(sum(gweights))
      error2 <- abs((new.dispersion - dispersion)/dispersion)
      dispersion <- new.dispersion
    }
    
    old.args <- args
    args <- abs.res/sqrt(dispersion)
    if (trace) {
      loglik <- -0.5 * length(abs.res) * log((dispersion)) + 
        sum(family$g0(abs.res/sqrt(dispersion), df = family$df, 
                      s = family$s, r = family$r, alpha = family$alpha, 
                      mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k))
      cat(" log-likelihood =", signif(loglik, 6))
    }
    error3 <- sqrt(sum((args - old.args)^2)/max(1e-20, sum(old.args^2)))
    if ((iter == maxit) || (max(error1, error2, error3, 
                                na.rm = TRUE) < epsilon)) 
      break
    iter <- iter + 1
  }
  
  if (trace) 
    cat("\n")
  if (maxit > 1 && iter == maxit){
    linearconv <- F
    warning(paste("\n linear convergence not obtained in", 
                  maxit, "iterations"))
  }
  coefs <- rep(NA, length(nas))
  coefs[!nas] <- start
  names(coefs) <- names(nas)
  names(dispersion) <- "dispersion"
  
  fitted <- as.vector(X %*% start + offset)
  
  residuals <- (Y - fitted)/sqrt(dispersion)
  w.1 <- family$g1(residuals, df = family$df, s = family$s, 
                   r = family$r, alpha = family$alpha, mp = family$mp, 
                   epsi = family$epsi, sigmap = family$sigmap, k = family$k)
  w.2 <- -2 * w.1 * gweights
  if (any(w.2 < 0)) 
    cat("\n --- negative iterative weights returned! --- \n")
  
  if (is.null.disp) {
    rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    Xd <- cbind(X, residuals)
  }
  dimnames(Xd)[[2]] <- c(Rnames, "scale")
  nn <- is.null(Rnames)
  Rnames <- list(dimnames(Xd)[[2]], dimnames(Xd)[[2]])
  R <- t(Xd) %*% diag(gweights) %*% Xd
  if (is.null.disp) 
    R[rank + 1, rank + 1] <- R[rank + 1, rank + 1] + length(residuals)
  attributes(R) <- list(dim = dim(R))
  if (!nn) 
    attr(R, "dimnames") <- Rnames
  loglik <- -0.5 * length(residuals) * log((dispersion)) + 
    sum(family$g0(residuals, df = family$df, s = family$s, 
                  r = family$r, alpha = family$alpha, mp = family$mp, 
                  epsi = family$epsi, sigmap = family$sigmap, k = family$k))
  names(loglik) <- NULL
  fit <- list(coefficients = coefs, dispersion = dispersion, gweights = gweights,
              fixed = !is.null.disp, residuals = residuals, fitted.values = fitted, 
              loglik = loglik, convergence = linearconv, Wg = family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), Wgder = family$g5(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), v = -2 * family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), rank = rank, R = as.matrix(R),  iter = iter - 
              1, scale = 4 * family$g2(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), scaledispersion = -1 + 4 * family$g3(args, 
              df = family$df, r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), scalevariance = family$g4(args, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), df = if (charmatch(family$family, 
              "Student", F)) family$df, s = if (charmatch(family$family, 
              "Gstudent", F)) family$s, r = if (charmatch(family$family, 
              "Gstudent", F)) family$r, alpha = if (charmatch(family$family, 
              "Glogis", F)) family$alpha, mp = if (charmatch(family$family, 
              "Glogis", F)) family$m, epsi = if (charmatch(family$family, 
              "Cnormal", F)) family$epsi, sigmap = if (charmatch(family$family, 
              "Cnormal", F)) family$sigmap, k = if (charmatch(family$family, 
              "Powerexp", F)) family$k, Xmodel = matrix(Xd[, (1:rank)],nrow(Xd), rank))
  fit
}



elliptical.envelope <- function (object, B = 100, arg,...) 
{
  initial <- NULL
  X <- model.matrix(object$terms)
  Xd <- as.matrix(object$Xmodel)
  n <- nrow(Xd)
  p <- ncol(Xd)
  family <- object$family
  control <- object$control
  #  ro <- object$resid
  #  tdf <- ro/sqrt(object$scalevariance)
  H <- Xd %*% solve(t(Xd) %*% Xd) %*% t(Xd)
  H1 <- (1/(object$scalevariance * object$scale)) * H
  varr <- object$scalevariance * object$dispersion * (diag(1, n) - H1)
  varr <- diag(varr)
  ri <- object$y - object$fitted
  tdf <- ri/sqrt(varr)
  e <- matrix(0, n, B)
  mu <- object$fitted
  phi <- object$dispersion
  resp <- NULL
  method <- "elliptical.fit.envel"
  elliptical.fitter <- get(method)
  #  offset = object$offset
  #  if (length(offset) == 1 && offset == 0) 
  #    offset <- rep(0, nobs)
  
  for (i in 1:B) {
    dist <- object$family[[1]]
    if (charmatch(dist, "Normal", F)) {
      resp <- rnorm(n, 0, 1)
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = Normal(), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "Cauchy", F)) {
      resp <- rcauchy(n, 0, 1)
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = Cauchy(), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "Student", F)) {
      resp <- rt(n, arg)
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = Student(arg), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "Gstudent", F)) {
      resp <- rgstudent(n, arg[1], arg[2])
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = Gstudent(arg), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "LogisI", F)) {
      stop(paste("not implemented yet"))
      resp <- rlogisI(n, 0, 1)
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = LogisI(), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "LogisII", F)) {
      resp <- rlogisII(n)
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = LogisII(), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "Glogis", F)) {
      stop(paste("not implement yet"))
      resp <- rglogis(n, arg[1], arg[2])
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = Glogis(arg), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "Cnormal", F)) {
      stop(paste("not implemented yet"))
      resp <- rcnormal(n, arg[1], arg[2])
      #fit <- elliptical(resp ~ X + (-1), family = Cnormal(arg), 
      #                  control = glm.control(maxit = 1000))
    }
    else if (charmatch(dist, "Powerexp", F)) {
      resp <- rpowerexp(n, arg)
      resp <- mu + sqrt(phi) * resp
      #fit <- elliptical(resp ~ X + (-1), family = Powerexp(arg), 
      #                  control = glm.control(maxit = 1000))
    }
    
    fit <- elliptical.fitter(X = X, Y = resp, family = family, dispersion = NULL, offset = 
                               NULL, maxit = control$maxit, epsilon = control$epsilon,
                             trace = control$trace, ...)
    ro <- fit$resid
    td <- ro/sqrt(fit$scalevariance)
    Xd <- as.matrix(fit$Xmodel)
    H <- Xd %*% solve(t(Xd) %*% Xd) %*% t(Xd)
    H1 <- (1/(fit$scalevariance * fit$scale)) * H
    varr <- fit$scalevariance * fit$dispersion * (diag(1, n) - H1)
    varr <- diag(varr)
    ri <- resp - fit$fitted
    td <- ri/sqrt(varr)
    e[, i] <- sort(td)
  }
  
  e1 <- numeric(n)
  e2 <- numeric(n)
  e3 <- numeric(n)
  e4 <- numeric(n)
  e5 <- numeric(n)
  e6 <- numeric(n)
  e7 <- numeric(n)
  for (i in 1:n) {
    eo <- sort(e[i, ])
    e1[i] <- eo[ceiling(B * 0.05)]
    e2[i] <- eo[ceiling(B * 0.95)]
  }
  e3 <- t(t(apply(e, 2, mean)))
  e4 <- t(t(apply(e, 2, vari)))
  e5 <- t(t(apply(e, 2, skewn)))
  e6 <- t(t(apply(e, 2, kurt)))
  e7 <- cbind(e3, e4, e5, e6)
  desc <- apply(e7, 2, mean)
  med <- apply(e, 1, mean)
  faixa <- range(tdf, e1, e2)
  screen(4)
  par(pty = "s")
  qqnorm(tdf, xlab = "Quantiles of N(0,1)", ylab = "Standardized residual", 
         ylim = faixa, pch = 16, main = '')
  par(new = TRUE)
  qqnorm(e1, axes = F, xlab = "", ylab = "", type = "l", ylim = faixa, 
         lty = 1, main = '')
  par(new = TRUE)
  qqnorm(e2, axes = F, xlab = "", ylab = "", type = "l", ylim = faixa, 
         lty = 1, main = '')
  par(new = TRUE)
  qqnorm(med, axes = F, xlab = "", ylab = "", type = "l", 
         ylim = faixa, lty = 2, main = '')
  x <- list(mean = desc[1], var = desc[2], skewness = desc[3], 
            kurtosis = desc[4])
  invisible(x)
}





elliptical.fit.envel <- function (X, Y, offset, family, dispersion, 
                            maxit, epsilon, trace) 
{
  n <- nrow(X)
  if (is.null(offset)) {
    offset <- rep(0, n)
  }
  
  p <- ncol(X)
  aux.model <- glm.fit(x = X, y = Y, offset = offset, 
                       family = gaussian())
  attr(aux.model, "class") <- c("glm", "lm")
  start <- aux.model$coef
  
  
  is.null.disp <- is.null(dispersion)
  elliptical.disp <- !is.null.disp && !is.number(dispersion)
  if (is.null.disp) 
    dispersion <- (summary(aux.model)$dispersion)
  if (elliptical.disp) 
    dispersion <- (summary(aux.model)$dispersion)
  args <- resid(aux.model)/sqrt(dispersion)
  
  if (any(nas <- is.na(start))) {
    names(nas) <- dimnames(X)[[2]]
    X <- X[, !nas]
    aux.model <- glm.fit(x = X, y = Y, offset = offset, 
                         family = gaussian())
    attr(aux.model, "class") <- c("glm", "lm")
    start <- aux.model$coef
    dispersion <- (summary(aux.model)$dispersion)
  }
  
  
  iter <- 1
  error2 <- error3 <- 0
  repeat {
    if (trace) 
      cat("\n iteration", iter, ":")
    {
      w.1 <- family$g1(args, df = family$df, r = family$r, 
                       s = family$s, alpha = family$alpha, mp = family$mp, 
                       epsi = family$epsi, sigmap = family$sigmap, 
                       k = family$k)
      dg <- family$g2(args, df = family$df, r = family$r, 
                      s = family$s, alpha = family$alpha, mp = family$mp, 
                      epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      fg <- family$g3(args, df = family$df, r = family$r, 
                      s = family$s, alpha = family$alpha, mp = family$mp, 
                      epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k)
      
      y.aux <- Y - offset
      w.h <- as.vector(-2 * w.1)
      aux.model <- glm.fit(x = X, y = y.aux, weights = w.h, 
                           family = gaussian())
      attr(aux.model, "class") <- c("glm", "lm")
      new.start <- coef(aux.model)
      }
    error1 <- max(abs((new.start - start)/start))
    start <- new.start
    abs.res <- Y - X %*% start - offset
    
    if (is.null.disp) {
      aux.dispersion <- dispersion
      new.dispersion <- mean((-2 * w.1) * abs.res^2)
      error2 <- abs((new.dispersion - dispersion)/dispersion)
      dispersion <- new.dispersion
    }
    old.args <- args
    args <- abs.res/sqrt(dispersion)
    if (trace) {
      loglik <- -0.5 * length(abs.res) * log((dispersion)) + 
        sum(family$g0(abs.res/sqrt(dispersion), df = family$df, 
                      s = family$s, r = family$r, alpha = family$alpha, 
                      mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
                      k = family$k))
      cat(" log-likelihood =", signif(loglik, 6))
    }
    error3 <- sqrt(sum((args - old.args)^2)/max(1e-20, sum(old.args^2)))
    if ((iter == maxit) || (max(error1, error2, error3, 
                                na.rm = TRUE) < epsilon)) 
      break
    iter <- iter + 1
  }
  if (trace) 
    cat("\n")
  if (maxit > 1 && iter == maxit) 
    warning(paste("\n linear convergence not obtained in", 
                  maxit, "iterations"))
  coefs <- rep(NA, length(nas))
  coefs[!nas] <- start
  names(coefs) <- names(nas)
  names(dispersion) <- "dispersion"
  fitted <- as.vector(X %*% start + offset)
  
  
  residuals <- (Y - fitted)/sqrt(dispersion)
  w.1 <- family$g1(residuals, df = family$df, s = family$s, 
                   r = family$r, alpha = family$alpha, mp = family$mp, 
                   epsi = family$epsi, sigmap = family$sigmap, k = family$k)
  w.2 <- -2 * w.1
  if (any(w.2 < 0)) 
    cat("\n --- negative iterative weights returned! --- \n")
  
  if (is.null.disp) {
    Xd <- cbind(X, residuals) ; rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    dimnames(Xd)[[2]] <- c(Rnames, "scale")
  } else {
    Xd <- X ; rank <- dim(X)[2]
    Rnames <- dimnames(X)[[2]]
    dimnames(Xd)[[2]] <- Rnames
  }
  
  nn <- is.null(Rnames)
  Rnames <- list(dimnames(Xd)[[2]], dimnames(Xd)[[2]])
  R <- t(Xd) %*% Xd
  if (is.null.disp) 
    R[rank + 1, rank + 1] <- R[rank + 1, rank + 1] + length(residuals)
  attributes(R) <- list(dim = dim(R))
  if (!nn) 
    attr(R, "dimnames") <- Rnames
  loglik <- -0.5 * length(residuals) * log((dispersion)) + 
    sum(family$g0(residuals, df = family$df, s = family$s, 
                  r = family$r, alpha = family$alpha, mp = family$mp, 
                  epsi = family$epsi, sigmap = family$sigmap, k = family$k))
  names(loglik) <- NULL
  fit <- list(coefficients = coefs, dispersion = dispersion, 
              fixed = !is.null.disp, residuals = residuals, fitted.values = fitted, 
              loglik = loglik, Wg = family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), Wgder = family$g5(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), v = -2 * family$g1(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), rank = rank, R = as.matrix(R), iter = iter - 
              1, scale = 4 * family$g2(residuals, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), scaledispersion = -1 + 4 * family$g3(args, 
              df = family$df, r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), scalevariance = family$g4(args, df = family$df, 
              r = family$r, s = family$s, alpha = family$alpha, 
              mp = family$mp, epsi = family$epsi, sigmap = family$sigmap, 
              k = family$k), df = if (charmatch(family$family, 
              "Student", F)) family$df, s = if (charmatch(family$family, 
              "Gstudent", F)) family$s, r = if (charmatch(family$family, 
              "Gstudent", F)) family$r, alpha = if (charmatch(family$family, 
              "Glogis", F)) family$alpha, mp = if (charmatch(family$family, 
              "Glogis", F)) family$m, epsi = if (charmatch(family$family, 
              "Cnormal", F)) family$epsi, sigmap = if (charmatch(family$family, 
              "Cnormal", F)) family$sigmap, k = if (charmatch(family$family, 
              "Powerexp", F)) family$k, Xmodel = matrix(Xd[, (1:rank)], 
              nrow(Xd), rank))
  fit
}



vari <- function (x) 
{
  wnas <- x[!is.na(x)]
  var(x, na.rm = TRUE) * (length(wnas) - 1)/length(wnas)
}



skewn <- function (x, na.rm = F, method = "fisher") 
{
  method <- char.expand(method, c("fisher", "moment"), stop("argument 'method' must match either \"fisher\" or \"moment\""))
  if (na.rm) {
    wnas <- x[!is.na(x)]
    if (length(wnas)) 
      x <- wnas
  }
  else if (any(is.na(x[!is.na(x)]))) 
    return(NA)
  n <- length(x)
  if (method == "fisher" && n < 3) 
    return(NA)
  x <- x - mean(x)
  if (method == "moment") 
    (sum(x^3)/n)/(sum(x^2)/n)^1.5
  else ((sqrt(n * (n - 1))/(n - 2)) * (sum(x^3)/n))/((sum(x^2)/n)^1.5)
}

kurt <- function (x, na.rm = F, method = "fisher") 
{
  method <- char.expand(method, c("fisher", "moment"), stop("argument 'method' must match either \"fisher\" or \"moment\""))
  if (na.rm) {
    wnas <- x[!is.na(x)]
    if (length(wnas)) 
      x <- wnas
  }
  else if (any(is.na(x[!is.na(x)]))) 
    return(NA)
  n <- length(x)
  if (method == "fisher" && n < 4) 
    return(NA)
  x <- x - mean(x)
  if (method == "moment") 
    (sum(x^4)/n)/(sum(x^2)/n)^2 - 3
  else ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - 
                               (3 * (n - 1))/(n + 1)))/((n - 2) * (n - 3))
}

