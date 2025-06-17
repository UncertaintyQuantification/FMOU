# This file is used to plot the Fig 7 
library(geometry)
library(plot3D)
library(MASS)
library(FastGaSP)

# functions used for NIF and FMOU
source('../functions/FMOU.R')
source('../functions/functions_geophysics_gppca.R')

## Read data and results summarized from Noel's MATLAB source
# load Greens functions
GF = read.csv("../real_data/Cascadia_data/Green_fun_large_mesh.csv", header=FALSE)
# load observations in 3 directions over 88 days (in meters)
obs = read.csv("../real_data/Cascadia_data/obs_KF.csv", header=FALSE)
# reads date (epochs) of observations
doy = read.csv("../real_data/Cascadia_data/DOY.csv", header=FALSE)
# get the name of stations we used in analysis:
names = read.csv("../real_data/Cascadia_data/station_names.csv", header=FALSE)
obs_uncertainty = read.csv("../real_data/Cascadia_data/obs_uncertainty.csv", header=FALSE)

GF = GF_ori = as.matrix(GF)
obs = as.matrix(obs)
names = as.vector(as.matrix(names))
doy = as.numeric(doy)
obs_uncertainty_ori = obs_uncertainty = as.matrix(obs_uncertainty)

epoch = dim(obs)[2]
N_site = dim(obs)[1] / 3 # number of GPS stations
N_xi = dim(GF)[2] # number of fault patches

# Stations in the index list: (5  16  17  24  30  80  83  88 104 105 107 109 110 111 112 113 114) are missing
missing_location = c(5,  16,  17,  24,  30,  80,  83,  88, 104, 105, 107, 109, 110, 111, 112, 113, 114)
# We only model the measurments in East and North
N_site_used = N_site - length(missing_location)
GF = GF[-c(3*missing_location, 3*missing_location-1, 3*missing_location-2), ]
GF= GF[-seq(3, 3*N_site_used, by=3),]
observed_location = (1:N_site)[-missing_location]
obs_uncertainty = obs_uncertainty[c(mapply(c, 3*observed_location-2, 3*observed_location-1)), ]

obs_used = obs
# fill the missing parts at each station:
t = c(1: epoch)
for(i in 1:length(observed_location)){
  index_here = c(3*observed_location[i]-2, 3*observed_location[i]-1, 3*observed_location[i])
  for(j in 1:3){
    if(sum(obs[index_here[j], ] == "NaN")>0){
      NA_index = as.numeric(which(obs[index_here[j], ] == "NaN"))
      input_here = t[-NA_index]
      output_here = as.numeric(obs[index_here[j], - NA_index])
      m_rgasp = rgasp(input_here, output_here, zero.mean="Yes", nugget.est=T, kernel_type = "pow_exp", alpha=c(1))
      
      ## fill the missing values:
      if (length(NA_index)>1){
        pred.model = predict(m_rgasp, testing_input=as.numeric(t[NA_index]))
        set.seed(i*j)
        obs_used[index_here[j], NA_index] = rnorm(length(pred.model$mean), mean=pred.model$mean, sd=pred.model$sd)
      }
      else{ 
        pred.model=predict(m_rgasp, testing_input=as.numeric(c(t[1], t[NA_index])))
        obs_used[index_here[j], NA_index] = rnorm(1, mean=pred.model$mean[2], sd=pred.model$sd[2])
      }
    }
  }
}

obs_cleaned = obs_used[c(mapply(c, 3*observed_location-2, 3*observed_location-1)), ] 
sum(is.na(obs_cleaned))
# Make the initial value of observations as zero (then cummulative displacements will be zero at t_1)
output = obs_cleaned  - obs_cleaned[, 1] 

# Try to add frame motions, where the displacements at all stations in each direction 
# is the same. It's equivalinet to model Y = G*s + (f_1, f_2)_{300 x 88} + \epsilon
G = cbind(GF, rep(c(1,0), 100), rep(c(0,1), 100)) # add 2 columns for East/North frame motions
dim(G)


eigen_Phi=eigen(G%*%t(G))
Lambda = eigen_Phi$values

# We can use the observation uncertainty provided to approximate the errors approximately
obs_uncertainty_E = mean(obs_uncertainty[seq(1,200,2), ], na.rm=TRUE) # 4.151058e-06
obs_uncertainty_N = mean(obs_uncertainty[seq(2,200,2), ], na.rm=TRUE) # 4.205995e-06
avg_obs_uncertainty = 0.5*mean(obs_uncertainty[seq(1,200,2), ], na.rm=TRUE) + 0.5*mean(obs_uncertainty[seq(2,200,2), ], na.rm=TRUE) # 4.178526e-06

U = eigen_Phi$vectors
Psi = t(U)%*%G    # compute basis functions

# And we select "d" based on observations' uncertainties
d_list = seq(10, 100)
sd_slip = rep(NA, length(d_list))
noise_est = rep(NA, length(d_list))

time1 = Sys.time()
left = d_list[1]
right = rev(d_list)[1]
while(left <= right){
  mid = left + floor((right-left)/2)
  
  em_fit = EM_alg_FMOU(output, d=mid, est_U0=F, U0=U[,1:mid], est_sigma0_2=T, threshold=10^{-6}) 
  est_noise = em_fit$record_sigma2_0[length(em_fit$record_sigma2_0)]
  if(est_noise < avg_obs_uncertainty){
    # overfit, reduce d
    right = mid-1
  }
  else{
    # underfit, increase d
    left = mid + 1
  }
}
Nbasis_EM  = mid
time2 = Sys.time()
time_FMOU_select_d = time2 - time1
time_FMOU_select_d # 3.564413 secs


tilde_output = t(U[,1:Nbasis_EM])%*%output
system.time({
  em_alg = EM_alg_FMOU(output, d=Nbasis_EM, est_U0=F, est_sigma0_2=F, threshold = 10^{-6},
                       U0=U[,1:Nbasis_EM], M=100, sigma0_2 = avg_obs_uncertainty)
  Sp_latent_state_EM = t(Psi[1:Nbasis_EM,])%*%diag(1/Lambda[1:Nbasis_EM])%*%em_alg$Z
  fit_latent_state_EM = U[,1:Nbasis_EM]%*%em_alg$Z 
})[3]



################################### plot the real data ###################################
names_site_idx = rep(NA, 6)
# names_list = c("P374", "P396","P407","P408","P411","P414","P415","P425","PCOL")
names_list = c("P374", "P407","P408","P411","P415","P425")
for(l in 1:6){
  names_site_idx[l] = which(names == names_list[l])
}

east_north_data = obs[c(mapply(c, 3*observed_location-2, 3*observed_location-1)), ] 
cum_data = 100*(east_north_data - east_north_data[,1]) # in cm
#cum_data = 100*apply(output, 1, function(x){mean(x,na.rm=T)})
cum_pred_data = 100*(fit_latent_state_EM)

FMOU_pred_mean_uq = em_alg$pred_mean_var

# pdf('real_data_sites.pdf',height=8.55, width=5.25)
plot(as.numeric(coastline[,1]), as.numeric(coastline[,2]), type="l", lwd=0.5, lty=1, xlab="", ylab="", 
     xlim=c(-124.5, -122.5), ylim=c(43.5, 48.6), xaxt="n", asp=0.65)
mtext("Longitude", side = 1, outer = F, line = 2.2, cex=1.3)
mtext("Latitude", side = 2, outer = F, line = 2.2, , cex=1.3)
mtext("(a) Locations of GPS stations", side = 3, outer = F, line = 1, cex=1.45)
axis(1, at=c(-125, -124, -123, -122), labels=c(-125, -124, -123, -122))
#lines(lon, lat, pch=1, type="p", col="grey")
for(ix in 1:100){
  x_end <- lon_cleaned[ix] + cum_data[2 * ix - 1, 88] * 0.65
  y_end <- lat_cleaned[ix] + cum_data[2 * ix, 88 ] * 0.65
  x_pred_end <- lon_cleaned[ix] + cum_pred_data[2 * ix - 1, 88] * 0.65
  y_pred_end <- lat_cleaned[ix] + cum_pred_data[2 * ix, 88 ] * 0.65
  if(lon_cleaned[ix]> -124.6 & lon_cleaned[ix] < -121.5 & lat_cleaned[ix]>43.6 & lat_cleaned[ix]<48.5){
    #segments(x0=lon_cleaned[ix], y0=lat_cleaned[ix], x1=x_end, y1=y_end, col="black", lwd=1.5)
    #segments(x0=lon_cleaned[ix], y0=lat_cleaned[ix], x1=x_pred_end, y1=y_pred_end, col="blue", lwd=1.5)
    arrows(x0=lon_cleaned[ix], y0=lat_cleaned[ix], x1=x_end, y1=y_end, col="black", lwd=1.1, length=0.08, angle=10)
    #arrows(x0=lon_cleaned[ix], y0=lat_cleaned[ix], x1=x_pred_end, y1=y_pred_end, col="blue", lwd=1.1, length=0.1, angle=10)
    # draw.ellipse(x_end, y_end, 100*sqrt(obs_uncertainty[2*ix-1, 88])* 0.2, 100*sqrt(obs_uncertainty[2*ix, 88])* 0.2, 
    #              border="blue", lwd=0.5, lty=1)
    # draw.ellipse(x_end, y_end, 100*sqrt(FMOU_pred_mean_uq[2*ix-1, 88]+ obs_uncertainty_E)* 0.25, 
    #              100*sqrt(FMOU_pred_mean_uq[2*ix, 88]+ obs_uncertainty_N)* 0.25, 
    #              border="blue", lwd=0.5, lty=1)
  }
  if(lon_cleaned[ix]> -124.6 & lon_cleaned[ix] < -121.5 & lat_cleaned[ix]>43.6 & lat_cleaned[ix]<48.5){
    points(lon_cleaned[ix], lat_cleaned[ix], cex=0.85, pch=20)
  }
  # for(ix in names_site_idx){
  #   points(lon_cleaned[ix], lat_cleaned[ix], col="darkred", cex=1.1, pch =20)
  # }
}
#text(-124.35, 43.55, cex=0.9)
legend(x=-123, y=43.68, legend = "0.2 cm ", cex=0.9)
segments(x0=-122.85, y0=43.62, x1=-122.72, y1=43.62, col="black", lwd=3)
text(lon_cleaned[names_site_idx], lat_cleaned[names_site_idx], names_list, pos=4, cex=0.8, col="darkred")
# dev.off()


### Use this if we want to split the regions fro E and N by offset:
# pdf('real_data_east_and_north_with_FMOU_uncertainty.pdf',height=6.3, width=7.35)
# old_par = par()
# par(mfrow=c(2,3))
# par(mar=c(4.5, 3.8, 3.8, 2) + 0.1)
# for(j in 1: length(names_site_idx)){
#   plot_id = names_site_idx[j]
#   yrange = 100*c(min(c(east_north_data_diff[2*plot_id-1,], east_north_data_diff[2*plot_id,]), na.rm=T)-0.0025,
#                  max(c(east_north_data_diff[2*plot_id-1,], east_north_data_diff[2*plot_id,]),na.rm=T)+0.015)
#   # East
#   east_ub = 100*(em_alg$pred_mean_95upper[2*plot_id-1,] + sqrt(obs_uncertainty_E))
#   east_lb = 100*(em_alg$pred_mean_95lower[2*plot_id-1,] - sqrt(obs_uncertainty_E))
#   
#   plot(100*east_north_data_diff[2*plot_id-1,], type="p", pch=21, col=rgb(0, 0, 0.6, alpha = 0.8), cex=0.55,
#        xlab="", ylab="", main=names_list[j], ylim=yrange)
#   polygon(c(1:88, 88:1), c(east_ub, rev(east_lb)), col = rgb(0, 0, 0.6, alpha = 0.35), border = NA)
#   lines(100*fit_latent_state_EM[2*plot_id-1,], type="l", lwd=2, col=rgb(0, 0, 0.6, alpha = 1))
#   # points(1:88, 100*east_north_data_diff[2*plot_id-1,], col="cyan", pch=15, cex=0.55)
#   
#   # North (offset by 1 for clarity)
#   north_ub = 100*(em_alg$pred_mean_95upper[2*plot_id,] + sqrt(obs_uncertainty_N)) + 1
#   north_lb = 100*(em_alg$pred_mean_95lower[2*plot_id,] - sqrt(obs_uncertainty_N)) + 1
#   lines(100*east_north_data_diff[2*plot_id,]+1, type="p", pch=21, cex=0.55, col=rgb(1, 0.4, 0.7, alpha = 0.8))
#   polygon(c(1:88, 88:1), c(north_ub, rev(north_lb)), col = rgb(1, 0.4, 0.7, alpha = 0.35), border = NA, fillOddEven=T)
#   lines(100*fit_latent_state_EM[2*plot_id,]+1, type="l", lwd=2, col=rgb(1, 0.4, 0.7, alpha =1))
#   if(j==1){
#     legend("topleft", pch=c(21, 21), col=c(col=rgb(0, 0, 0.6, alpha = 0.8), rgb(1, 0.4, 0.7, alpha = 0.8)),
#            legend=c("East", "North"), cex=0.9, horiz=T)
#   }
# }
# par(old_par)
# par(mfrow=c(1,1))
# mtext("Days since 06/03/2011", side = 1, outer = F, line = 3.3, cex=1.05)
# mtext("Displacements (cm)", side = 2, outer = F, line = 3.3, cex=1.05)
# mtext("(b) Displacements at selected stations", side = 3, outer = F, line = 3, cex=1.235)
# # dev.off()



# pdf('real_data_east_and_north_with_FMOU_uncertainty.pdf',height=6.35, width=7.38)
old_par = par()
par(mfrow=c(2,3))
# par(mar=c(4.3, 3.9, 4.25, 1.8) + 0.1)
par(mar=c(3.6, 3.9, 3, 1.8) + 0.1, oma=c(1, 0, 1, 0))
for(j in 1: length(names_site_idx)){
  plot_id = names_site_idx[j]
  yrange = 100*c(min(c(east_north_data_diff[2*plot_id-1,], east_north_data_diff[2*plot_id,]), na.rm=T)-0.0025,
                 max(c(east_north_data_diff[2*plot_id-1,], east_north_data_diff[2*plot_id,]),na.rm=T)+0.0025)
  # East
  east_ub = 100*(em_alg$pred_mean_95upper[2*plot_id-1,] + sqrt(obs_uncertainty_E))
  east_lb = 100*(em_alg$pred_mean_95lower[2*plot_id-1,] - sqrt(obs_uncertainty_E))

  plot(100*east_north_data_diff[2*plot_id-1,], type="p", pch=21, col=rgb(0, 0, 0.6, alpha = 0.8), cex=0.55,
       xlab="", ylab="", main=names_list[j], ylim=yrange, xaxt="n")
  axis(1, at=c(0, 26, 52, 78), labels=c("Jun 3", "Jun 29", "Jul 4", "Aug 20"))
  polygon(c(1:88, 88:1), c(east_ub, rev(east_lb)), col = rgb(0, 0, 0.6, alpha = 0.35), border = NA)
  lines(100*fit_latent_state_EM[2*plot_id-1,], type="l", lwd=2, col=rgb(0, 0, 0.6, alpha = 1))
  # points(1:88, 100*east_north_data_diff[2*plot_id-1,], col="cyan", pch=15, cex=0.55)

  # North (offset by 1 for clarity)
  north_ub = 100*(em_alg$pred_mean_95upper[2*plot_id,] + sqrt(obs_uncertainty_N))
  north_lb = 100*(em_alg$pred_mean_95lower[2*plot_id,] - sqrt(obs_uncertainty_N))
  lines(100*east_north_data_diff[2*plot_id,], type="p", pch=21, cex=0.55, col=rgb(1, 0.4, 0.7, alpha = 0.8))
  polygon(c(1:88, 88:1), c(north_ub, rev(north_lb)), col = rgb(1, 0.4, 0.7, alpha = 0.35), border = NA, fillOddEven=T)
  lines(100*fit_latent_state_EM[2*plot_id,], type="l", lwd=2, col=rgb(1, 0.4, 0.7, alpha =1))
  if(j==1){
    legend("topleft", pch=c(21, 21), col=c(col=rgb(0, 0, 0.6, alpha = 0.8), rgb(1, 0.4, 0.7, alpha = 0.8)),
           legend=c("East", "North"), cex=0.9, horiz=T)
  }
}
par(old_par)
par(mfrow=c(1,1))
# mtext("Days since 06/03/2011", side = 1, outer = F, line = 3.3, cex=1.05)
mtext("Displacements (cm)", side = 2, outer = F, line = 3.3, cex=1.05)
mtext("(b) Displacements at selected stations", side = 3, outer = F, line = 3, cex=1.235)
# dev.off()
