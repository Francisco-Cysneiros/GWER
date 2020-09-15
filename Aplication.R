#----------#
# Packages #
#----------#
library(spgwr)          # Package required to obtain dataset map
library(gwer)           # Package for aplication of GWER models
library(maptools)       # Package of spatial analysis functions
library(e1071)          # Package for measures of skewness and kurtosis
library(spdep)          # Package for measure the Moran index and spatial autocorrelation tests
library(xtable)         # Package for generation of tables in LaTeX format
library(robustbase)     # Package for aplied adjusted boxplot



#--------------------------#
# Georgia Spatial Database #
#--------------------------#
# Dataset in shapefile format #
data(georgia)
spdata.aux <- gSRDF

# Response variable used #
PctBach <- spdata.aux$PctBach

# Explanatory variables used #
TotPop90 <- spdata.aux$TotPop90/100000
PctRural <- spdata.aux$PctRural
PctEld <- spdata.aux$PctEld
PctFB <- spdata.aux$PctFB
PctPov <- spdata.aux$PctPov
PctBlack <- spdata.aux$PctBlack



#----------------------------#
# Preliminary Specifications #
#----------------------------#
# Building the spatial dataset used #
varc <- data.frame(PctBach, TotPop90, PctRural, PctEld, PctFB, PctPov, PctBlack)
varnames <- names(varc) ; varf <- varc[, -c(4,7)] ; varnames.f <- names(varf)
spdata <- SpatialPolygonsDataFrame(spdata.aux, varc, match.ID = F)

# Coordinates specification #
coord <- coordinates(spdata) ; colnames(coord) <- c("Latitude","Longitude") 
distcoord <- knn2nb(knearneigh(coord, longlat = T)) ; col.test <- nb2listw(distcoord, style="W")

# Model fit specifications #
formula.ajuste <- PctBach ~ TotPop90 + PctRural + PctFB + PctPov
n <- nrow(spdata) ; p = ncol(spdata) ; band.method <- "aic" ; gweight.fit = gwr.Gauss
residuals.names <- 'Standardized residuals'

# Specification for prints of the results #
colres <- colorRampPalette(c("black", "darkblue", "lightblue", "darkblue", "black")) ; cortres <- seq(-8, 8, 2)
colpv <- colorRampPalette(c("white", "black")) ; cortpv <- c(0.00, 0.05, 1.00)

# Additional functions used #
panel.cor <- function(x, y, ...){
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 1)
}



#----------------------#
# Descriptive Analysis #
#----------------------#
# Obtaining a summary of the variables in study #
sumvar <- t(apply(varc, 2, summary))
cv <- round(rapply(varc,sd)/rapply(varc,mean),4)
skewn <- round(rapply(varc,skewness),4)
kurt <- round(rapply(varc,kurtosis),4)+3
corv <- cor(varc) ; sort(corv[1,-1]) ; corv > 0.75
summaryv <- cbind(sumvar[,-c(2,5)],cv,skewn,kurt)
colnames(summaryv) <- c('Min.', 'Median', 'Mean', 'Max.', 'CV', 'Skewness', 'Kurtosis')

# Prints summary of variables in study #
xtable(corv, caption = "Correlation of variables", auto=T, digits = 2)
xtable(summaryv, caption = "Descriptive analysis of variables", auto=T, digits = 2)

# Explanatory graphical analysis #
adjbox(varc[,1], ylab=varnames[1])
pairs(varc, lower.panel = panel.cor)

# Spatial distribution of the variables in study #
for(i in 1:length(varf)){
  print(spplot(spdata, varnames.f[i], main = varnames.f[i], xlab = colnames(coord)[1], ylab = colnames(coord)[2], colorkey = list(space = "right", height = 1), as.table = T, scales = list(draw = TRUE)))
}

# Moran indexes and spatial autocorrelation tests of the variables #
for(i in 1:length(varc)){
  moranvi <- moran.test(varc[,i], col.test) ; moranvi
  moran.plot(varc[,i], col.test, xlab = varnames[i], ylab = paste("Spatially lagged", varnames[i]))
  text(max(varc[,i]) - 1.5*sd(varc[,i]), min(varc[,i]) + 1.5*sd(varc[,i]), substitute(Moran.I == idex, list(idex = round(moranvi$estimate[1], 4))),cex = .8)
  text(max(varc[,i]) - 1.5*sd(varc[,i]), min(varc[,i]) + 1.2*sd(varc[,i]), substitute(p-Value == pv, list(pv = round(moranvi$p.value, 4))),cex = .8)
}



#-------------------------------------------------#
# Ellptical regression model (with normal errors) #
#-------------------------------------------------#
# Fit of the model #
ajuste.gn <- elliptical(formula.ajuste, family = Normal(), data = spdata, method = "elliptical.fit")

# Summary of the fitted model #
summary.gn <- summary(ajuste.gn) ; summary.gn 

# Calulating diagnostic measures #
diag.gn <- elliptical.diag(ajuste.gn)
res.gn <- diag.gn$rs ; fv.gn <- ajuste.gn$fitted.values

# Moran indexes and spatial autocorrelation tests of the standardized residuals #
moranr.gn <- moran.test(res.gn, col.test) ; moranr.gn
moran.plot(res.gn, col.test, main = 'Classical linear regression model (Normal)', xlab = residuals.names, ylab = paste(residuals.names, "Spatially lagged"))
text(max(res.gn) - sd(res.gn), min(res.gn) + 2.4*sd(res.gn), substitute(Moran.I == idex, list(idex = round(moranr.gn$estimate[1], 4))), cex = .8)
text(max(res.gn) - sd(res.gn), min(res.gn) + 2.0*sd(res.gn), substitute(p-Value == pv, list(pv = format(round(moranr.gn$p.value, 4), scientific = F))),cex = .8)

# Total local influence index plot for case-weight perturbation #
Cic.gn = diag.gn$Cic/sum(diag.gn$Cic)
plot(Cic.gn, ylim = c(0,0.6), main = 'Local influence for case-weight perturbation', ylab = expression(C[i]), xlab = 'Index')
abline(2/n,0) ; identify(Cic.gn,n=2)

# Total local influence index plot for scale perturbation #
Cih.gn = diag.gn$Cih/sum(diag.gn$Cih)
plot(Cih.gn, ylim = c(0,1), main = 'Local influence for perturbation in the scale', ylab = expression(C[i]), xlab='Index')
abline(2/n,0) ; identify(Cih.gn,n=2)



#----------------------------------------------------#
# Ellptical regression model (with t-Student errors) #
#----------------------------------------------------#
# Selection for degree of freedom #
dft1 <- rep(0,30)
for(i in 3:30){
  ajuste.gt <- elliptical(formula.ajuste, family = Student(df=i), data = spdata, method = "elliptical.fit")
  loglikgt <- ajuste.gt$loglik ; nobs <- dim(spdata)[1] ; rank <- dim(spdata)[2] 
  dft1[i] <- (2*rank - 2*loglikgt + 2*(rank)*(rank + 1)/(nobs - rank - 1))
}
dfgt <- which(dft1==min(dft1[-c(1,2)])) ; dfgt

# Fit of the model #
ajuste.gt <- elliptical(formula.ajuste, family = Student(df=dfgt), data = spdata, method = "elliptical.fit")

# Summary of the fitted model #
summary.gt <- summary(ajuste.gt) ; summary.gt

# Calulating diagnostic measures #
diag.gt <- elliptical.diag(ajuste.gt)
res.gt <- diag.gt$rs ; fv.gt <- ajuste.gt$fitted.values

# Moran indexes and spatial autocorrelation tests of the standardized residuals #
moranr.gt <- moran.test(res.gt, col.test) ; moranr.gt
moran.plot(res.gt, col.test, main = 'Elliptical linear regression model (t-Student)', xlab = residuals.names, ylab = paste(residuals.names, "Spatially lagged"))
text(max(res.gt) - 1.5*sd(res.gt), min(res.gt) + 2.5*sd(res.gt), substitute(Moran.I == idex, list(idex = round(moranr.gt$estimate[1], 4))),cex = .8)
text(max(res.gt) - 1.5*sd(res.gt), min(res.gt) + 2*sd(res.gt), substitute(p-Value == pv, list(pv = format(round(moranr.gt$p.value, 4), scientific = F))),cex = .8)

# Total local influence index plot for case-weight perturbation #
Cic.gt = diag.gt$Cic/sum(diag.gt$Cic)
plot(Cic.gt, ylim = c(0,0.6), main = 'Local influence for case-weight perturbation', ylab = expression(C[i]), xlab='Index')
abline(2/n,0)

# Total local influence index plot for scale perturbation #
Cih.gt = diag.gt$Cih/sum(diag.gt$Cih)
plot(Cih.gt, ylim = c(0,1), main = 'Local influence for perturbation in the scale', ylab = expression(C[i]), xlab='Index')
abline(2/n,0)



#---------------------------------#
# GWER model (with normal errors) #
#---------------------------------#
# Fit of the model #
gwr.bw <- gwer.sel(formula.ajuste, data = spdata, family = Normal(), method = band.method, longlat = TRUE, adapt = TRUE, gweight = gweight.fit)
ajuste.gwr <- gwer(formula.ajuste, family = Normal(), adapt = gwr.bw, parplot = FALSE, longlat = TRUE, gweight = gweight.fit,
                   hatmatrix = TRUE, spdisp = TRUE, data = spdata, method = "gwer.fit")

# Summary of the fitted model #
xtable(t(apply(ajuste.gwr$coef$est,2,summary)), caption = "Summary of local parameters estimates in normal GWER model", auto=T, digits = 2)
xtable(t(apply(ajuste.gwr$coef$se,2,summary)), caption = "Summary of local parameters standard errors in t-Student GWER model", auto=T, digits = 2)
print(ajuste.gwr)

# Calulating diagnostic measures #
diag.ln <- gwer.diag(ajuste.gwr)
res.ln = diag.ln$rs ; fv.ln <- ajuste.gwr$fitted

# Moran indexes and spatial autocorrelation tests of the standardized residuals #
moranr.ln <- moran.test(res.ln, col.test) ; moranr.ln
moran.plot(res.ln, col.test, main = 'Normal GWER model', xlab = residuals.names, ylab = paste("Spatially lagged", residuals.names))
text(max(res.ln) - 1.5*sd(res.ln), min(res.ln) + 2.5*sd(res.ln), substitute(Moran.I == idex, list(idex = format(round(moranr.ln$estimate[1], 4), scientific=F))), cex = .8)
text(max(res.ln) - 1.5*sd(res.ln), min(res.ln) + 2*sd(res.ln), substitute(Valor-p == pv, list(pv = format(round(moranr.ln$p.value, 4), scientific=F))), cex = .8)

# Total local influence index plot for case-weight perturbation #
Cic.ln <- diag.ln$Cic/sum(diag.ln$Cic)
plot(Cic.ln, ylim = c(0, 0.2), main = 'Local influence for case-weight perturbation', ylab = expression(C[i]), xlab = 'Idex')
abline(2/n, 0) ; identify(Cic.ln, n=2)

# Total local influence index plot for scale perturbation #
Cih.ln <- diag.ln$Cih/sum(diag.ln$Cih)
plot(Cih.ln, ylim = c(0, 0.2), main = 'Local influence for perturbation in the scale', ylab = expression(C[i]), xlab = 'Index')
abline(2/n, 0) ; identify(Cih.ln, n=2)



#-------------------------------#
# GWER model (t-Student errors) #
#-------------------------------#
# Selection for degree of freedom #
dft2 <- rep(0,30)
for(i in 3:30){
  gwer.bw <- gwer.sel(formula.ajuste, data=spdata, family = Student(df=i), method = band.method, longlat = TRUE, adapt = TRUE, gweight = gweight.fit)
  ajuste.gwer <- gwer(formula.ajuste, family = Student(df=i), adapt = gwer.bw, parplot = FALSE, longlat = TRUE, gweight = gweight.fit,
                      hatmatrix = TRUE, spdisp = TRUE, data=spdata, method = "gwer.fit", control = glm.control(epsilon = 1e-04, trace = F))
  dft2[i] <- ajuste.gwer$results$AICc
}
dflt <- which(dft2==min(dft2[-c(1,2)])) ; dflt

# Fit of the model #
gwer.bw <- gwer.sel(formula.ajuste, data=spdata, family = Student(df=dflt), method = band.method, longlat = TRUE, adapt = TRUE, gweight = gweight.fit)
ajuste.gwer <- gwer(formula.ajuste, family = Student(df=dflt), adapt = gwer.bw, parplot = FALSE, longlat = TRUE,gweight = gweight.fit,
                    hatmatrix = TRUE, spdisp = TRUE, data = spdata, method = "gwer.fit")

# Summary of the fitted model #
xtable(t(apply(ajuste.gwer$coef$est,2,summary)), caption = "Summary of local parameters estimates in t-Student GWER model", auto=T, digits = 2)
xtable(t(apply(ajuste.gwer$coef$se,2,summary)), caption = "Summary of local parameters standard errors in t-Student GWER model", auto=T, digits = 2)
print(ajuste.gwer)

# Calulating diagnostic measures #
diag.lt <- gwer.diag(ajuste.gwer)
res.lt <- diag.lt$rs ; fv.lt <- ajuste.gwer$fitted

# Moran indexes and spatial autocorrelation tests of the standardized residuals #
moranr.lt <- moran.test(res.lt, col.test) ; moranr.lt
moran.plot(res.lt, col.test, main = 't-Student GWER model', xlab = residuals.names, ylab = paste("Spatially lagged", residuals.names))
text(max(res.lt) - 1.5*sd(res.lt), min(res.lt) + 2.5*sd(res.lt), substitute(Moran.I == idex, list(idex = format(format(round(moranr.lt$estimate[1], 4)), scientific=F))), cex = .8)
text(max(res.lt) - 1.5*sd(res.lt), min(res.lt) + 2*sd(res.lt), substitute(Valor-p == pv, list(pv = format(round(moranr.lt$p.value, 4), scientific=F))), cex = .8)

# Total local influence index plot for case-weight perturbation #
Cic.lt <- diag.lt$Cic/sum(diag.lt$Cic)
plot(Cic.lt, ylim = c(0, 0.2), main = 'Local influence for case-weight perturbation', ylab = expression(C[i]), xlab = 'Index')
abline(2/n,0)

# Total local influence index plot for scale perturbation #
Cih.lt <- diag.lt$Cih/sum(diag.lt$Cih)
plot(Cih.lt, ylim = c(0, 0.2), main = 'Local influence for perturbation in the scale', ylab = expression(C[i]), xlab = 'Index')
abline(2/n,0)



#------------------------------------------------------#
# Preliminary Specifications for Confirmatory Analysis #
#------------------------------------------------------#
# Fit of the GWER models without observation #108 #
exc <- 108 ; obsc1 <- list('sp.polygons',spdata[exc,], density = 3)

gwr.bw.c1 <- gwer.sel(formula.ajuste, data=spdata[-exc,], family = Normal(), method = band.method, longlat = TRUE, adapt = TRUE, gweight = gweight.fit)
ajuste.gwr.c1 <- gwer(formula.ajuste, family = Normal(), adapt = gwr.bw.c1, parplot = FALSE, longlat = TRUE, gweight = gweight.fit,
                      hatmatrix = TRUE, spdisp = TRUE, data=spdata[-exc,], method = "gwer.fit")
print(ajuste.gwr.c1)  # Normal GWER model

gwer.bw.c1 <- gwer.sel(formula.ajuste, data=spdata[-exc,], family = Student(df=dflt), method = band.method, longlat = TRUE, adapt = TRUE, gweight = gweight.fit)
ajuste.gwer.c1 <- gwer(formula.ajuste, family = Student(df=dflt), adapt = gwer.bw.c1, parplot = FALSE, longlat = TRUE, gweight = gweight.fit,
                       hatmatrix = TRUE, spdisp = TRUE, data=spdata[-exc,], method = "gwer.fit")
print(ajuste.gwer.c1) # t-Student GWER model


# Standardizing the scales of the spatial surfaces of the estimated parameters #
cort.b <- matrix(0, (length(varf) + 1), 13)
for(i in 1:(length(varf) + 1)){
  beta.lower <- min(ajuste.gwr$coef$est[,i],ajuste.gwr.c1$coef$est[,i],ajuste.gwer$coef$est[,i],ajuste.gwer.c1$coef$est[,i]) - .001
  beta.upper <- max(ajuste.gwr$coef$est[,i],ajuste.gwr.c1$coef$est[,i],ajuste.gwer$coef$est[,i],ajuste.gwer.c1$coef$est[,i]) + .001
  cort.b[i, ] <- seq(beta.lower, beta.upper, (beta.upper-beta.lower)/12) 
}


#-----------------------------------------------------#
# Confirmatory Analysis for GwER (with normal errors) #
#-----------------------------------------------------#
# Spatial surface of the parameters estimations for normal errors with all observation #
est.mapce0.ln <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwr$coef$est[,1]), match.ID = F)
spplot(est.mapce0.ln, zcol = 1, main = 'The intercept coeficients of normal GWER model', at = cort.b[1,], colorkey = list(space = "right", height = 1, at = cort.b[1,]), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  est.mapce1.ln <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwr$coef$est[,i]), match.ID = F) 
  print(spplot(est.mapce1.ln, zcol = 1, main = paste('The ', varnames.f[i], ' coeficients of normal GWR model', sep = ''), at = cort.b[i,], colorkey = list(space = "right", height = 1, at = cort.b[i,]), as.table = T, scales = list(draw = TRUE)))
}

est.mapcp.ln <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwr$coef$est[,length(varf)+1]), match.ID = F)
spplot(est.mapcp.ln, zcol = 1, main = 'The dispersion estimates of normal GWER model', at = cort.b[length(varf)+1,], colorkey = list(space = "right", height = 1, at = cort.b[length(varf)+1,]), as.table = T, scales = list(draw = TRUE)) 

# Spatial surface for the p-value of the parameters for normal errors with all observation #
pv.mapce0.ln <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwr$coef$pvalue[,1]), match.ID = F) 
spplot(pv.mapce0.ln, zcol = 1, main = 'p-Value of the intercept coeficients for normal GWER model', col.regions = colpv(4), at = cortpv, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  pv.mapce1.ln <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwr$coef$pvalue[,i]), match.ID = F) # Complete data
  print(spplot(pv.mapce1.ln, zcol = 1, main = paste('p-Value of the ', varnames.f[i], ' coeficients for normal GWR model', sep = ''), col.regions = colpv(4), at = cortpv, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)))
}


# Spatial surface of the parameters estimations for normal errors without observation #108 #
est.mapc0.ln <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwr.c1$coef$est[,1]), match.ID = F) 
spplot(est.mapc0.ln, zcol = 1, main = paste('The intercept coeficients of normal GWER model (-', exc, ')', sep = ''), at = cort.b[1,], sp.layout = obsc1, colorkey = list(space = "right", height = 1, at = cort.b[1,]), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  est.mapc1.ln <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwr.c1$coef$est[,i]), match.ID = F) # Data without observation #27
  print(spplot(est.mapc1.ln, zcol = 1, main = paste('The ', varnames.f[i], ' coeficients of normal GWR model (-', exc, ')', sep = ''), at = cort.b[i,], sp.layout = obsc1, colorkey = list(space = "right", height = 1, at = cort.b[i,]), as.table = T, scales = list(draw = TRUE)))
}

est.mapcp.ln <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwr.c1$coef$est[,length(varf)+1]), match.ID = F)
spplot(est.mapcp.ln, zcol = 1, main = paste('The dispersion estimates of normal GWER model (-', exc, ')', sep = ''), at = cort.b[length(varf)+1,], sp.layout = obsc1, colorkey = list(space = "right", height = 1, at = cort.b[length(varf)+1,]), as.table = T, scales = list(draw = TRUE)) 

# Spatial surface for the p-value of the parameters for normal errors without observation #108 #
pv.mapc0.ln <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwr.c1$coef$pvalue[,1]), match.ID = F) 
spplot(pv.mapc0.ln, zcol = 1, main = paste('p-Value of the intercept coeficients for normal GWER model (-', exc, ')', sep = ''), col.regions = colpv(4), at = cortpv, sp.layout = obsc1, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  pv.mapc1.ln <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwr.c1$coef$pvalue[,i]), match.ID = F) # Data without observation #27
  print(spplot(pv.mapc1.ln, zcol = 1, main = paste('p-Value of the ', varnames.f[i], ' coeficients for normal GWR model (-', exc, ')', sep = ''), col.regions = colpv(4), at = cortpv, sp.layout = obsc1, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)))
}




#--------------------------------------------------------#
# Confirmatory Analysis for GwER (with t-Student errors) #
#--------------------------------------------------------#
# Spatial surface of the parameters estimations for t-Student errors with all observation #
est.mapce0.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$est[,1]), match.ID = F) 
spplot(est.mapce0.lt, zcol = 1, main = 'The intercept coeficients of t-Student GWER model', at = cort.b[1,], colorkey = list(space = "right", height = 1, at = cort.b[1,]), as.table=T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  est.mapce1.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$est[,i]), match.ID = F)
  print(spplot(est.mapce1.lt, zcol = 1, main = paste('The ', varnames.f[i], ' coeficients of GWER model', sep = ''), at = cort.b[i,], colorkey = list(space = "right", height = 1, at = cort.b[i,]), as.table = T, scales = list(draw = TRUE)))
}

est.mapcep.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$est[,length(varf)+1]), match.ID = F)
spplot(est.mapcep.lt, zcol = 1, main = 'The dispersion estimates of t-Student GWER model', at = cort.b[length(varf)+1,], colorkey = list(space = "right", height = 1, at = cort.b[length(varf)+1,]), as.table = T, scales = list(draw = TRUE)) 

# Spatial surface for the p-value of the parameters for t-Student errors with all observation #
pv.mapce0.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$pvalue[,1]), match.ID = F) 
spplot(pv.mapce0.lt, zcol = 1, main = 'p-Value of the intercept coeficients for normal GWER model', col.regions = colpv(4), at = cortpv, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  pv.mapce1.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$pvalue[,i]), match.ID = F)
  print(spplot(pv.mapce1.lt, zcol = 1, main = paste('p-Value of the ', varnames.f[i], ' coeficients of GWER model', sep = ''), col.regions = colpv(4), at = cortpv, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)))
}


# Intercept - t-Student GWER model #
est.mapc0.lt <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwer.c1$coef$est[,1]), match.ID = F) 
spplot(est.mapc0.lt, zcol = 1, main = paste('The intercept coeficients of t-Student GWER model (-', exc, ')', sep = ''), at = cort.b[1,], sp.layout = obsc1, colorkey = list(space = "right", height = 1, at = cort.b[1,]), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  est.mapc1.lt <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwer.c1$coef$est[,i]), match.ID = F) 
  print(spplot(est.mapc1.lt, zcol = 1,  main = paste('The ', varnames.f[i], ' coeficients of GWER model (-', exc, ')', sep = ''), at = cort.b[i,], sp.layout = obsc1, colorkey = list(space = "right", height = 1, at = cort.b[i,]), as.table = T, scales = list(draw = TRUE)))
}

est.mapcp.lt <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwer.c1$coef$est[,length(varf)+1]), match.ID = F)
spplot(est.mapcp.lt, zcol = 1, main = paste('The dispersion estimates of t-Student GWER model (-', exc, ')', sep = ''), at = cort.b[length(varf)+1,], sp.layout = obsc1, colorkey = list(space = "right", height = 1, at = cort.b[length(varf)+1,]), as.table = T, scales = list(draw = TRUE)) 

# Spatial surface for the p-value of the parameters for t-Student errors without observation #108 #
pv.mapc0.lt <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwer.c1$coef$pvalue[,1]), match.ID = F) 
spplot(pv.mapc0.lt, zcol = 1, main = paste('p-Value of the intercept coeficients for t-Student GWER model (-', exc, ')', sep = ''), col.regions = colpv(4), at = cortpv, sp.layout = obsc1, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  pv.mapc1.lt <- SpatialPolygonsDataFrame(spdata.aux[-exc,], data.frame(ajuste.gwer.c1$coef$pvalue[,i]), match.ID = F) 
  print(spplot(pv.mapc1.lt, zcol = 1,  main = paste('The ', varnames.f[i], ' coeficients of GWER model (-', exc, ')', sep = ''),col.regions = colpv(4), at = cortpv, sp.layout = obsc1, colorkey = list(space = "right", height = 1, col.regions = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)))
}



#-----------------------------------------------------#
# Interpretation of GWER model (with t-Student error) #
#-----------------------------------------------------#
# Estimate of the coeficients #
est.map0.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$est[,1]), match.ID = F) # Estimates of beta_0 parameters
spplot(est.map0.lt, zcol = 1, main = 'Estimates of intercept coefficients', colorkey = list(space = "right", height = 1), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  est.map1.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$est[,i]), match.ID = F) # Estimates of beta_1 parameters
  print(spplot(est.map1.lt, zcol = 1, main = paste('Estimates of ', varnames.f[i], ' coefficients', sep = ''), colorkey = list(space = "right", height = 1), as.table=T, scales = list(draw = TRUE)))
}

est.mapp.lt <- SpatialPolygonsDataFrame(spdata.aux, data.frame(ajuste.gwer$coef$est[,6]), match.ID = F) # Estimates of phi parameters
spplot(est.mapp.lt, zcol = 1, main = 'Estimates of dispersion parameters', colorkey = list(space = "right", height = 1), as.table=T, scales = list(draw = TRUE)) 

# p-Values of the signifance for the coefficients #
pv.map0.lt <- SpatialPolygonsDataFrame(spdata.aux,data.frame(ajuste.gwer$coef$pvalue[,1]), match.ID = F) # p-Values of beta_0 parameters
spplot(pv.map0.lt, zcol = 1, main = 'p-Values of intercept coefficients', col.regions = colpv(4), at = cortpv, colorkey = list(space = "right", height = 1, col = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)) 

for(i in 2:length(varf)){
  pv.map1.lt <- SpatialPolygonsDataFrame(spdata.aux,data.frame(ajuste.gwer$coef$pvalue[,i]), match.ID = F) # p-Values of beta_1 parameters
  print(spplot(pv.map1.lt, zcol = 1, main = paste('p-Values of ', varnames.f[i], ' coefficients', sep = ''), col.regions = colpv(4), at = cortpv, colorkey = list(space = "right", height = 1, col = colpv(4), at = cortpv), as.table = T, scales = list(draw = TRUE)))
}



#----------------------------------------------#
# Simulated Envelope of Standardized residuals #
#----------------------------------------------#
# Setting the seed of simulations #
set.seed(1111)

# Loading the functions of simulated envelope #
source("envelope.R")

# Simulated envelope of normal GWER #
gwer.envelope(ajuste.gwr, B = 100, ident = 2)

# Simulated envelope of t-Student GWER #
gwer.envelope(ajuste.gwer, arg = 3, B = 100)

