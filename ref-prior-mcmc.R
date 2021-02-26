library(rgdal); library(sp); library(maps); library(maptools); library(dplyr); library(geoR)
library(usmap); library(ggplot2); library(reshape2); library(mvtnorm)
source("lonlat_smap_data_script.R")
# mydata contains the data
names(mydata) = c("longitude", "latitude", "aug6",
                  "aug7", "aug8", "aug9",
                  "aug10", "aug11", "aug12")
# The following code to convert coordinates to state names is from 
# https://stackoverflow.com/questions/8751497/latitude-longitude-coordinates-to-state-code-in-r
lonlat_to_state_sp <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

mydata$state = lonlat_to_state_sp(mydata[,1:2])

dak.min = mydata %>%
  dplyr::select(longitude, latitude, state, aug6, aug7, aug8) %>%
  filter(state %in% c("north dakota", "south dakota", "minnesota"))

transformed_data = usmap_transform(dak.min)
melted = melt(
  transformed_data,
  id.vars = c("longitude", "latitude", "state", "longitude.1", "latitude.1"),
  variable.name = "date"
)

soil = melted[complete.cases(melted),] 

# Raw data
soil_raw = soil[,c(1,2,7)]

# Averaging replicates
avgd.data = data.frame(dak.min[,1:2], value = rowMeans(dak.min[,4:6], na.rm = T))
soil_avg = avgd.data[complete.cases(avgd.data),]

# Jittering replicates
soil_jitter = data.frame(jitterDupCoords(soil[,1:2], max=1e-7), value=soil$value)

# Lists for geoR
in_geor = list(coords=soil_raw[,1:2], data=soil_raw$value)
in_geor_avg = list(coords=soil_avg[,1:2], data=soil_avg$value)
in_geor_jitter = list(coords=soil_jitter[,1:2], data=soil_jitter$value)

# Creating design matrix D
Y = soil_jitter$value
n = length(Y)
lon = soil_jitter$longitude
lat = soil_jitter$latitude
D = cbind(rep(1,n), lon, lat, lon^2, lat^2, lon*lat)

n = dim(D)[1]
k = dim(D)[2]

# Hyperparameters for sigma^2
a = 2
b = 1

thetaphi.post = function(t, p) {
  theta = 1/(1+exp(-t))
  phi = exp(p)
  # theta = sigma^2/(tau^2 + sigma^2)
  I = diag(n)
  rho.mat = varcov.spatial(in_geor_jitter$coords, cov.pars=c(1, phi), nugget=0, kappa=0.5)
  rho = rho.mat$varcov
  V = theta*rho + (1-theta)*I
  
  # Z = L^-1 Y
  # F = L^-1 D
  L = t(chol(V))
  L.inv = forwardsolve(L, I)
  Z = forwardsolve(L, log(Y))
  FF = forwardsolve(L, D)
  DVD = t(FF) %*% FF
  
  # R_1 beta.hat = (Q^T Y)_1
  # S^2 = ||(Q^TY)_2||
  QR = qr(FF)
  Q = qr.Q(QR, complete=T)
  R = qr.R(QR, complete=T)
  
  QZ = t(Q) %*% Z
  
  R1 = R[1:k,]
  QZ1 = QZ[1:k,]
  QZ2 = QZ[(k+1):n,]
  
  beta.hat = backsolve(R1, QZ1)
  S2 = t(QZ2) %*% QZ2
  
  # Find inverse of (D V^-1 D)^-1 for beta covariance
  # F^-1 = R^-1 Q^T
  F.inv = backsolve(R, t(Q))
  DVD.inv = F.inv %*% t(F.inv)
  
  # Determinants
  E = t(chol(DVD))
  logdet.V = 2*sum(log(diag(L)))
  logdet.DVD = 2*sum(log(diag(E)))
  
  # Prior
  P = I - D %*% (F.inv %*% t(F.inv)) %*% (t(FF) %*% L.inv)
  K = rho - I
  r = as.matrix(dist(D[,1:2], upper=T, diag=T))
  J = theta/phi^2 * r * rho
  V.inv = t(L.inv) %*% L.inv
  W1 = K %*% V.inv %*% P
  W2 = J %*% V.inv %*% P
  trW1 = sum(diag(W1))
  trW2 = sum(diag(W2))
  trW1sq = sum(diag(W1 %*% W1))
  trW2sq = sum(diag(W2 %*% W2))
  trW1W2 = sum(diag(W1 %*% W2))
  ref.mat = matrix(c(n-k, trW1, trW2,
                     trW1, trW1sq, trW1W2,
                     trW2, trW1W2, trW2sq), ncol=3, nrow=3, byrow=T)
  prior = determinant(ref.mat)
  
  logdensity = - 0.5*logdet.V - 0.5*logdet.DVD - ((n-k)/2 + a)*log(S2 + 2*b) + 0.5*prior$modulus
  out = list(logdensity=logdensity, beta.hat=beta.hat, S2=S2, DVD.inv=DVD.inv)
  return(out)
}

n.iter = 11000
beta.mcmc = matrix(c( -102, -0.392, 3.58, 0.00141, -0.0212, 0.0158), ncol=6, nrow=n.iter+1, byrow=T)
sigma.mcmc = rep(0.13, n.iter+1)
phi.mcmc = rep(2, n.iter+1)
theta.mcmc = rep(2, n.iter+1)

sigma.tune = 1

for (i in 1:n.iter) {
  prop = rmvnorm(1, c(1,1), sigma=sigma.tune*diag(2))
  phi.prop = prop[1]
  theta.prop = prop[2]
  
  pi.p = thetaphi.post(theta.prop, phi.prop)
  pi.c = thetaphi.post(theta.mcmc[i], phi.mcmc[i])
  
  log.ratio = pi.p$logdensity - pi.c$logdensity + dmvnorm(c(theta.mcmc[i], phi.mcmc[i]), sigma=sigma.tune*diag(2)) - dmvnorm(prop, sigma=sigma.tune*diag(2))
  
  if (log(runif(1)) < log.ratio) {
    phi.mcmc[i+1] = exp(phi.prop)
    theta.mcmc[i+1] = 1/(1+exp(-theta.prop))
    sigma.mcmc[i+1] = 1/rgamma(1, (n+k)/2+a, (pi.p$S2 + 2*b)/2)
    beta.mcmc[i+1,] = rmvnorm(1, pi.p$beta.hat, sigma.mcmc[i+1]*pi.p$DVD.inv)
  } else {
    phi.mcmc[i+1] = phi.mcmc[i]
    theta.mcmc[i+1] = theta.mcmc[i]
    sigma.mcmc[i+1] = sigma.mcmc[i]
    beta.mcmc[i+1,] = beta.mcmc[i,]
  }
  if (i %% 100 == 0) print(i)
}

write.csv(cbind(beta.mcmc, sigma.mcmc, theta.mcmc, phi.mcmc), "refprior_samples2.csv")

