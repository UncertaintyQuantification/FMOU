## this file is used to generate figure 8 for the real data

library(geometry)
library(plot3D)
library(MASS)
library(FastGaSP)
library(dplyr)
library(ggplot2)
library(gridExtra)

# 1. load estimated slip rate data
slip_rate_NIF = read.csv('saved_data/slip_rate_NIF_large_mesh.csv')*100/365
slip_rate_NIF_selected_d = read.csv('saved_data/slip_rate_NIF_large_mesh_select_d.csv')*100/365
slip_rate_NIF_modified = read.csv('saved_data/slip_rate_modified_NIF_large_mesh.csv')*100/365
slip_rate_FMOU = read.csv('saved_data/slip_rate_FMOU_large_mesh.csv')*100/365


# 2. load slip's lat and lon slips
nd_ll = read.csv("../real_data/Cascadia_data/nd_ll_large_mesh.csv", header=FALSE)
el = read.csv("../real_data/Cascadia_data/el_large_mesh.csv", header=FALSE)
lon = read.csv("../real_data/Cascadia_data/lon_G.csv", header=FALSE)
lat = read.csv("../real_data/Cascadia_data/lat_G.csv", header=FALSE)

nd_ll = as.matrix(nd_ll)
el = as.matrix(el)
lon = as.matrix(lon)
lat = as.matrix(lat)

slip_lat_lon = matrix(0, nrow=dim(slip_rate_FMOU)[1] , 2) # longitude and latitude
for(i in 1:dim(slip_rate_FMOU)[1]){
  slip_lat_lon[i,1]= mean(nd_ll[el[i, ], 1])
  slip_lat_lon[i,2]= mean(nd_ll[el[i, ], 2])
}
slip_lat_lon = as.data.frame(slip_lat_lon)
names(slip_lat_lon) = c("lon", "lat")

# 3. load tremor's lat and lon, remove rows where the region are outside slips boundary
tremor = read.csv("../Cascadia_data/tremor_2011.csv", header=T)
tremor = as.data.frame(tremor)

tremor_cleaned = tremor %>%
  mutate(lat = as.numeric(lat),
         lon = as.numeric(lon),) %>%
  filter(lat>min(slip_lat_lon$lat), lat<max(slip_lat_lon$lat),
         lon>min(slip_lat_lon$lon), lon<max(slip_lat_lon$lon))

# 4. define boundary
lon_min = min(slip_lat_lon$lon)
lon_max = max(slip_lat_lon$lon)
lat_min = min(slip_lat_lon$lat)
lat_max = max(slip_lat_lon$lat)

# 5. create grid
ratio_lat_lon=(lat_max-lat_min)/(lon_max-lon_min)
lon_breaks <- seq(lon_min, lon_max, length.out = 40)
lat_breaks <- seq(lat_min, lat_max, length.out =floor(45*ratio_lat_lon))

tremor_cleaned = tremor_cleaned %>%
  mutate(lat_index = cut(lat, breaks = lat_breaks, labels = FALSE, include.lowest = TRUE),
         lon_index = cut(lon, breaks = lon_breaks, labels = FALSE, include.lowest = TRUE),
         index = paste(lat_index, lon_index, sep="-"))
slip_lat_lon = slip_lat_lon %>%
  mutate(lat_index = cut(lat, breaks = lat_breaks, labels = FALSE, include.lowest = TRUE),
         lon_index = cut(lon, breaks = lon_breaks, labels = FALSE, include.lowest = TRUE),
         index = paste(lat_index, lon_index, sep="-"))

# 6. cutoff for high slip rate
slip_cutoff_set = seq(0.0025, 0.02,by=0.0025)

# 7. create dates data - every 5 days
data_start_date = as.Date("2011-06-03")
five_days_start_date = seq(from = data_start_date+1, by = 5, length.out = ceiling(88 / 5))[1:16]
five_days_end_date = seq(from = data_start_date + (5), by = 5, length.out = ceiling(88 / 5))[1:16]
tt_start = as.numeric(five_days_start_date-data_start_date)
tt_end = as.numeric(five_days_end_date-data_start_date)

period_name = paste(five_days_start_date,five_days_end_date,sep="to" )
method_name = c("FMOU","NIF","NIF_modified","NIF_selected_d")

# 8. define metrics
count_grid_high_slip_list = as.list(1:length(slip_cutoff_set)) # how many grids contain high slip rates
count_grid_tremor = rep(NA, length(period_name)) # how many grids contain tremor
count_grid_identified_tremor_list = as.list(1:length(slip_cutoff_set)) # how many tremors were identify

names(count_grid_tremor) = period_name

# 9. compute metric under different cutoff
for(index_cutoff in 1:length(slip_cutoff_set)){
  print(index_cutoff)
  slip_cutoff = slip_cutoff_set[index_cutoff]
  
  # create matrix to save results for current cutoff
  count_grid_high_slip <- matrix(NA, length(period_name), 4)
  count_grid_identified_tremor = matrix(NA, length(period_name), 4)
  rownames(count_grid_high_slip) = period_name
  colnames(count_grid_high_slip) = method_name
  rownames(count_grid_identified_tremor) = period_name
  colnames(count_grid_identified_tremor) = method_name
  
  for(ts in 1:length(period_name)){
    tremor_five_days = tremor_cleaned %>% filter(date>=five_days_start_date[ts] & date<=five_days_end_date[ts])
    count_grid_tremor[ts] = length(unique(tremor_five_days$index))
    slip_rate_five_days_FMOU = apply(slip_rate_FMOU[,tt_start[ts]:tt_end[ts]],1,mean)
    slip_rate_five_days_NIF = apply(slip_rate_NIF[,tt_start[ts]:tt_end[ts]],1,mean)
    slip_rate_five_days_NIF_modified = c(apply(slip_rate_NIF_modified[,tt_start[ts]:tt_end[ts]],1,mean),NA)
    slip_rate_five_days_NIF_selected_d = apply(slip_rate_NIF_selected_d[,tt_start[ts]:tt_end[ts]],1,mean)
    
    high_slip_rate_index_FMOU = which(slip_rate_five_days_FMOU>slip_cutoff)
    high_slip_rate_index_NIF = which(slip_rate_five_days_NIF>slip_cutoff)
    high_slip_rate_index_NIF_modified = which(slip_rate_five_days_NIF_modified>slip_cutoff)
    high_slip_rate_index_NIF_selected_d = which(slip_rate_five_days_NIF_selected_d>slip_cutoff)
    
    # compute metrics
    count_grid_high_slip[ts,'FMOU'] = length(unique(slip_lat_lon[high_slip_rate_index_FMOU, "index"]))
    count_grid_high_slip[ts, 'NIF'] = length(unique(slip_lat_lon[high_slip_rate_index_NIF, "index"]))
    count_grid_high_slip[ts, 'NIF_modified'] = length(unique(slip_lat_lon[high_slip_rate_index_NIF_modified, "index"]))
    count_grid_high_slip[ts, 'NIF_selected_d'] = length(unique(slip_lat_lon[high_slip_rate_index_NIF_selected_d, "index"]))
    
    count_grid_identified_tremor[ts, 'FMOU'] = sum(unique(tremor_five_days$index) %in% unique(slip_lat_lon[high_slip_rate_index_FMOU, "index"]))
    count_grid_identified_tremor[ts, 'NIF'] = sum(unique(tremor_five_days$index) %in% unique(slip_lat_lon[high_slip_rate_index_NIF, "index"]))
    count_grid_identified_tremor[ts, 'NIF_modified'] = sum(unique(tremor_five_days$index) %in% unique(slip_lat_lon[high_slip_rate_index_NIF_modified, "index"]))
    count_grid_identified_tremor[ts, 'NIF_selected_d'] = sum(unique(tremor_five_days$index) %in% unique(slip_lat_lon[high_slip_rate_index_NIF_selected_d, "index"]))
  }
  # save the results in a list
  count_grid_high_slip_list[[index_cutoff]] = count_grid_high_slip
  count_grid_identified_tremor_list[[index_cutoff]] = count_grid_identified_tremor
}


# 10. summary results

# number of grids that has tremor
count_grid_tremor

# number of grids that has high slip rate
count_grid_high_slip_summary = matrix(NA, length(method_name), length(slip_cutoff_set))
for(index_cutoff in 1:length(slip_cutoff_set)){
  count_grid_high_slip_summary[,index_cutoff] = apply(count_grid_high_slip_list[[index_cutoff]], 2, mean)
}

rownames(count_grid_high_slip_summary) = method_name
colnames(count_grid_high_slip_summary) = slip_cutoff_set

prop_tremor_identified_summary = matrix(NA, length(method_name), length(slip_cutoff_set))
for(index_cutoff in 1:length(slip_cutoff_set)){
  prop_tremor_identified_summary[,index_cutoff] = apply((count_grid_identified_tremor_list[[index_cutoff]]/count_grid_tremor)[-3,], 2, function(x) mean(x, na.rm=T))
}
rownames(prop_tremor_identified_summary) = method_name
colnames(prop_tremor_identified_summary) = slip_cutoff_set


## 11. Visualization
library(tidyr)
prop_tremor_plot_df = as.data.frame(prop_tremor_identified_summary) %>%
  mutate(method=method_name) %>%
  pivot_longer(cols=-method,
               names_to="threshold",
               values_to = "prop") %>%
  filter(method!="NIF_selected_d") %>%
  mutate(method=ifelse(method=="NIF_modified", "Modified NIF", method),
         method=factor(method, levels=c("FMOU","NIF","Modified NIF")),
         threshold=as.numeric(threshold))

real_data_p1 <- ggplot(prop_tremor_plot_df) + 
  geom_line(aes(x=threshold, y=prop*100, group=method,col=method)) +
  geom_point(aes(x=threshold, y=prop*100, group=method,col=method)) + 
  scale_x_continuous(breaks = c(0.0025,  0.0125,0.0200)) + 
  xlab("threshold of slip rates (cm/day)") + 
  ylab("percentage") + 
  scale_color_manual(
    values = c(
      "FMOU" = "#00a6fb",
      "NIF" = "#EE4E4E",
      "Modified NIF" = "#FF9C73")) +
  theme_bw() +
  theme(legend.position = c(0.20, 0.25),
        legend.title= element_blank(),
        legend.text = element_text(size = 8, angle=0),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size=12, hjust=0.5),
        legend.margin = margin(0.01, 0.01, 0.01, 0.01)) +
  ggtitle("(a) % of detected tremor events")
real_data_p1

count_grid_high_slip_plot_df = as.data.frame(count_grid_high_slip_summary) %>%
  mutate(method=method_name) %>%
  pivot_longer(cols=-method,
               names_to="threshold",
               values_to = "count_grid") %>%
  filter(method!="NIF_selected_d") %>%
  mutate(method=ifelse(method=="NIF_modified", "Modified NIF", method),
         method=factor(method, levels=c("FMOU","NIF","Modified NIF")),
         threshold=as.numeric(threshold))

real_data_p2 <- ggplot(count_grid_high_slip_plot_df) + 
  geom_line(aes(x=threshold, y=count_grid, group=method,col=method)) +
  geom_point(aes(x=threshold, y=count_grid, group=method,col=method)) +
  scale_x_continuous(breaks = c(0.0025, 0.0125,0.0200)) + 
  xlab("threshold of slip rates (cm/day)") + 
  ylab("number of grids") + 
  scale_color_manual(
    values = c(
      "FMOU" = "#00a6fb",
      "NIF" = "#EE4E4E",
      "Modified NIF" = "#FF9C73"
    )) +
  theme_bw() + 
  theme(legend.position="none",
        #legend.position = c(0.79, 0.74),
        #legend.title= element_blank(),
        legend.text = element_text(size = 8, angle=0),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size=12, hjust=0.5)) + 
  #legend.margin = margin(0.1, 0.1, 0.1, 0.1))  +
  ggtitle("(b) # grids having high slip rates")
real_data_p2

computation_time <- data.frame(method=c("FMOU","Modified NIF","NIF"),
                               computational_time = c(0.92+3.564, 13898.872,1560.47)) %>%
  mutate(method=factor(method, levels=c("FMOU","NIF","Modified NIF")))

real_data_p3 <- ggplot(computation_time,aes(x=method, y=computational_time,fill=method)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(computational_time, 2)),  # Add values on top
            vjust = -0.3, size = 4)  +
  ylab("running time (sec)") +
  ylim(0, 14800) + 
  scale_fill_manual(
    values = c(
      "FMOU" = "#00a6fb",
      "NIF" = "#EE4E4E",
      "Modified NIF" = "#FF9C73"
    )) +
  theme_bw() + 
  ggtitle("(c) running time (seconds)") + 
  theme(legend.position="none",
        plot.title = element_text(size=12, hjust=0.5))
real_data_p3

#pdf("../figures/real_data_identified_alarm.pdf", height=2.3,width=9)
grid.arrange(real_data_p1,real_data_p2, real_data_p3 , nrow=1)
#dev.off()
