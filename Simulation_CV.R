#----------#
# Packages #
#----------#
library(spgwr)  # Package required to obtain dataset map
library(gwer)   # Package for aplication of GWER models



#--------------------------#
# Simulation Specification #
#--------------------------#
# Seed used for this simulation #
set.seed(1111)

# Preliminary definitions and dataset construction #
data(georgia) ; map <- gSRDF ; n <- length(map)
map.coord <- coordinates(map) ; order.map <- map[order(map.coord[,2]),]
order.coord <- map.coord[order(map.coord[,2], decreasing = FALSE),]
max.coord <- apply(abs(map.coord), 2, max) ; min.coord <- apply(abs(map.coord), 2, min)
range.coord <- ((max.coord[1] - min.coord[1]) + (max.coord[2] - min.coord[2]))
repl  <- 1000 ; df.s  <- 3 ; p <- 1 ; length.coef <- range.coord/4 ; band.method = 'cv' ; outlier.pos = floor(n/2)
i <- cont.n <- cont.t <- 0 ; x <- rnorm(n, 0, 1)

# Definition of the real parameters #
beta.0 <- 1 + (1/length.coef)*(abs(order.coord[, 1]) + abs(order.coord[, 2]))
beta.1 <- 1 + (1/length.coef)*(abs(order.coord[, 1]) + rev(abs(order.coord[, 2])))
beta <- cbind(beta.0, beta.1) 
sp.beta <- SpatialPolygonsDataFrame(order.map, data = as.data.frame(beta))
spplot(sp.beta, scales = list(draw = TRUE))

# The output matrix #
Estat <- matrix(0, repl, 15)
colnames(Estat) <- c('RMSE.B0.n', 'RMSE.B0.t', 'RMSE.B1.n', 'RMSE.B1.t', 'RSS', 'RSS_gn', 'RSS_n', 'RSS_t', 'RSS_gt',
                     'CE.n', 'CE.t', 'band.n', 'band.t', 'Outlier.Lat', 'Outlier.Long')



#------------------------#
# Monte Carlo Simulation #
#------------------------#
while(i < repl){
  # Specification for the i-th simulation #
  bw.n <- bw.t <- NULL
  timings <- list()
  
  # Generation of response variable #
  error <- rnorm(n, 0, 1)
  y <- beta.0 + beta.1*x + error

  # Generation of response variable #
  y[outlier.pos] <- max(y) + 3*sd(y)
  outlier.coord <- map.coord[y == max(y)]

  # Spatial dataset structure #
  data <- cbind(y, x) ; data <- data.frame(data) ; colnames(data) <- c('y','x')
  spdata <- SpatialPolygonsDataFrame(order.map, data = data)
  
  # Fit elliptical classical models #
  ajuste.gn <- elliptical(y ~ x, family = Normal(), data = spdata)
  ajuste.gt <- elliptical(y ~ x, family = Student(df=df.s), data = spdata)

  # Fit GWER models (with timer) #
  timings[["startn"]] <- Sys.time()
  bw.n <- try(gwer.sel(y ~ x, data = spdata, family = Normal(), method = band.method, verbose = F, adapt = T, longlat = T), silent = T)
  ajuste.n <- try(gwer(y ~ x, family = Normal(), adapt = bw.n, data=spdata, longlat = T, hatmatrix = TRUE), silent = T)
  timings[["stopn"]] <- Sys.time()
  
  timings[["startt"]] <- Sys.time()
  bw.t <- try(gwer.sel(y ~ x, data = spdata, family = Student(df=df.s), method = band.method, verbose = F, adapt = T, longlat = T), silent = T)
  ajuste.t <- try(gwer(y ~ x, family = Student(df=df.s), adapt = bw.t, data=spdata, longlat = T, hatmatrix = TRUE), silent = T)
  timings[["stopt"]] <- Sys.time()
  
  # Convergence errors count  #
  if(attr(ajuste.n,"class") == "try-error" || attr(ajuste.t,"class") == "try-error"){
    if(attr(ajuste.n,"class") == "try-error")
      cont.n <- cont.n + 1
    
    if(attr(ajuste.cv.t,"class") == "try-error")
      cont.t <- cont.t + 1
    
    next
  }
  i=i+1
  
  # Estimated coefficients #
  beta.est <- cbind(beta[, 1], ajuste.n$coef$est[, 1], ajuste.t$coef$est[, 1], beta[, 2], ajuste.n$coef$est[, 2], ajuste.t$coef$est[, 2])
  sp.beta.est <- SpatialPolygonsDataFrame(order.map, data = as.data.frame(beta.est))
  names.plot <- c('Real intercept', 'GWR intercept', 'GWER intercept', 'Real coefficient of x', 'GWR coefficient of x', 'GWER coefficient of x')
  
  # Prints the spatial surface of estimated coefficients #
  postscript(paste("robust_g_georgia_n_", n, "_band_", band.method, "_model_", i, "_beta.eps", sep = ''))
  print(spplot(sp.beta.est, names.attr = names.plot, as.table = TRUE, scales = list(draw = TRUE)))
  dev.off()
  
  # Evaluation methods #
  RMSE.n <- sqrt(apply((beta-ajuste.n$coef$est)^2, 2, mean))
  RMSE.t <- sqrt(apply((beta-ajuste.t$coef$est)^2, 2, mean))
  RSS <- sum(error^2) ; V.n <- rep(1, n) ; V.t <- (df.s+1)/(df.s+residuals(ajuste.t, type = 'response')^2)
  RSS.gn <- sum((y - ajuste.gn$fitted)^2*ajuste.gn$v) ; RSS.gt <- sum((y - ajuste.gt$fitted)^2*ajuste.gt$v)
  RSS.n <- sum((y - ajuste.n$fitted)^2*V.n) ; RSS.t <- sum((y - ajuste.t$fitted)^2*V.t)
  CE.n <- difftime(timings[["stopn"]], timings[["startn"]], units = 'hour') ; CE.t <- difftime(timings[["stopt"]], timings[["startt"]], units = 'hour')
  
  # Writing the results matrix #
  Estat[i, ] <- c(RMSE.n[1], RMSE.t[1], RMSE.n[2], RMSE.t[2], RSS, RSS.gn, RSS.n, RSS.gt, RSS.t, CE.n, CE.t, bw.n, bw.t, outlier.coord)
  file <- paste("robust_g_georgia_n_", n, "_band_", band.method, "matrix.txt", sep = '')
  write.table(round(Estat, 4), file, sep=" ", col.names = TRUE)
}

# Writing number of no convergênce for the GWER models #
nconv <- c(cont.n, cont.t) ; names(nconv) <- c('nconv.n', 'nconv.t')
file.nconv <- paste("robust_g_georgia_n_", n, "_band_", band.method, "_nconv.txt", sep = '')
write.table(nconv, file.nconv, sep=" ", col.names = TRUE)

# Prints boxplots of the RMSE for the GWER models #
est.plot <- Estat[, 1:4] ; names.est.plot <- c(expression(GWR~beta[0]), expression(GWER~beta[0]), expression(GWR~beta[1]), expression(GWER~beta[1]))
postscript(paste("robust_g_georgia_n_", n, "_rmse_beta.eps", sep = ''))
boxplot(est.plot, names = names.est.plot, ylab = 'RMSE')
dev.off()

# Prints boxplots of the RSS for the classical models and GWER models #
pred.plot <- Estat[, 5:9] ; names.pred.plot <- c(expression(sum(epsilon^2)), 'LM', 'GWR', 't-LM', 'GWER')
postscript(paste("robust_g_georgia_n_", n, "_rss_pred1.eps", sep = ''))
boxplot(pred.plot, names = names.pred.plot, ylab = 'RSS')
dev.off()

# Prints boxplots of the RSS only for the GWER models #
pred.plot2 <- Estat[, c(5, 7, 9)] ; names.pred2.plot <- c(expression(sum(epsilon^2)), 'GWR', 'GWER')
postscript(paste("robust_g_georgia_n_", n, "_rss_pred2.eps", sep = ''))
boxplot(pred.plot2, names = names.pred2.plot, ylab = 'RSS')
dev.off()

# Prints boxplots of the optimal bandwidths obtained in GWER models #
band.plot <- Estat[, 12:13] ; names.band.plot <- c('GWR', 'GWER')
postscript(paste("robust_g_georgia_n_", n, "_band_model.eps", sep = ''))
boxplot(band.plot, names = names.band.plot, ylab = 'Bandwidth')
dev.off()

# Prints boxplots of the computational time for the GWER models #
time.plot <- Estat[, 10:11] ; names.time.plot <- c('GWR', 'GWER')
postscript(paste("robust_g_georgia_n_", n, "_time_model.eps", sep = ''))
boxplot(time.plot, names = names.time.plot, ylab = 'Time')
dev.off()