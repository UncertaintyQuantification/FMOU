grid_slip_update = function(t0, current_time, semi_minor, semi_major_ini, center=c(0, 0)){
  slip_grid = matrix(0, nrow=length(y_grid), ncol=length(x_grid))
  delta_t = current_time - t0
  semi_major = delta_t * major_grow_rate + semi_major_ini
  for(ix in 1:length(x_grid)){
    for(iy in 1:length(y_grid)){
      dis_from_center = (x_grid[ix]-center[1])^2/semi_minor^2 + (y_grid[iy]-center[2])^2/semi_major^2
      if(dis_from_center >1){
        slip_grid[length(y_grid)-iy+1, ix] = NA
      }
      else{
        slip_grid[length(y_grid)-iy+1, ix] = c*(1-(dis_from_center)^{alpha}/r)
      }
    }
  }
  return(slip_grid)
}

grid_slip_rate_update = function(t0, current_time, semi_minor, semi_major_ini, center=c(0, 0)){
  slip_rate_grid = matrix(0, nrow=length(y_grid), ncol=length(x_grid))
  delta_t = current_time - t0
  semi_major = delta_t * major_grow_rate + semi_major_ini
  for(ix in 1:length(x_grid)){
    for(iy in 1:length(y_grid)){
      dis_from_center = (x_grid[ix]-center[1])^2/semi_minor^2 + (y_grid[iy]-center[2])^2/semi_major^2
      if(dis_from_center >1){
        slip_rate_grid[length(y_grid)-iy+1, ix] = NA
      }
      else{
        slip_rate_grid[length(y_grid)-iy+1, ix] = 2*alpha*major_grow_rate* (y_grid[iy]-center[2])^2 * dis_from_center^{alpha-1} /(semi_major^3)
        slip_rate_grid[length(y_grid)-iy+1, ix] = slip_rate_grid[length(y_grid)-iy+1, ix] *c
      }
    }
  }
  return(slip_rate_grid)
}


point_slip_update = function(t0, current_time, semi_minor, semi_major_ini, center=c(0, 0)){
  slip_grid = rep(0, dim(el)[1])
  delta_t = current_time - t0
  semi_major = delta_t * major_grow_rate + semi_major_ini
  for(ix in 1:dim(el)[1]){
    dis_from_center = (tricenter[ix,1]-center[1])^2/semi_minor^2 + (tricenter[ix,2]-center[2])^2/semi_major^2
    if(dis_from_center > 1){
      slip_grid[ix] =0 
    }
    else{
      slip_grid[ix] = c*(1-(dis_from_center)^{alpha}/r)
    }
  }
  return(slip_grid)
}


point_slip_rate_update = function(t0, current_time, semi_minor, semi_major_ini, center=c(0, 0)){
  slip_rate_grid = rep(0, dim(el)[1])
  delta_t = current_time - t0
  semi_major = delta_t * major_grow_rate + semi_major_ini
  for(ix in 1:dim(el)[1]){
    dis_from_center = (tricenter[ix,1]-center[1])^2/semi_minor^2 + (tricenter[ix,2]-center[2])^2/semi_major^2
    if(dis_from_center > 1){
      slip_rate_grid[ix] =0
    }
    else{
      slip_rate_grid[ix] = 2*alpha*major_grow_rate* (tricenter[ix,2]-center[2])^2 * dis_from_center^{alpha-1} /(semi_major^3)
      slip_rate_grid[ix] = slip_rate_grid[ix]*c
    }
  }
  return(slip_rate_grid)
}