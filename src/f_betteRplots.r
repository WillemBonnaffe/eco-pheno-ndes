######
##  ##
######

## Goal:

## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

###############
## FUNCTIONS ##
###############

## add_axis_and_grid
## Goal: Add custom axis and grid to plot given the vector of x and y data.
## Arguments:
## * x - vector - vector of x coordinates of the data to plot
## * y - vector - vector of y coordinates of the data to plot
add_axes_and_grid = function(x, y, alpha, label)
{
  ## format data x
  alpha_x = alpha[1]
  lb_x = floor(min(x)/alpha_x)*alpha_x
  rb_x = ceiling(max(x)/alpha_x)*alpha_x
  dx = 2.5
  x = seq(lb_x, rb_x, dx * alpha_x)
  
  ## format data y
  alpha_y = alpha[2]
  lb_y = floor(min(y)/alpha_y)*alpha_y
  rb_y = ceiling(max(y)/alpha_y)*alpha_y
  dy = 2.5
  y = seq(lb_y, rb_y, dy * alpha_y)
  
  ## background
  coords = par("usr")
  coords_x = coords[1:2]
  coords_y = coords[3:4]
  polygon(x=c(coords_x, rev(coords_x)), y=c(c(coords_y[1],coords_y[1]), c(coords_y[2],coords_y[2])), col=adjustcolor("lightgrey",alpha=0.2), border=NA)
  
  ## grid guides
  for (l in 1:length(y)) lines(c(x[1]-10,x[length(x)]+10), c(y[l], y[l]), col="white")
  for (l in 1:length(x)) lines(c(x[l], x[l]), c(y[1]-10,y[length(y)]+10), col="white")
  
  ## x axis
  axis(1, label=x, at=x, lwd=0, lwd.ticks=1)
  axis(2, label=y, at=y, lwd=0, lwd.ticks=1)
  
  ## add the label above the plot
  x <- par("usr")[2] - 0.2  # Adjust the x-coordinate as needed
  y <- par("usr")[4] + 0.05  # Adjust the y-coordinate to position it above the plot
  mtext(text = label, side = 3, line = 1, at = x, cex = 1.25)
  
}

## add_axis_and_grid_at
## Goal: Add custom axis and grid to plot given the vector of x and y data.
## Arguments:
## * x - vector - vector of x coordinates of the data to plot
## * y - vector - vector of y coordinates of the data to plot
add_axes_and_grid_at = function(x, y, at_x, at_y, label)
{
 
  ## background
  coords = par("usr")
  coords_x = coords[1:2]
  coords_y = coords[3:4]
  polygon(x=c(coords_x, rev(coords_x)), y=c(c(coords_y[1],coords_y[1]), c(coords_y[2],coords_y[2])), col=adjustcolor("lightgrey",alpha=0.2), border=NA)
  
  ## grid guides
  for (l in 1:length(at_y)) lines(c(at_x[1]-10,at_x[length(at_x)]+10), c(at_y[l], at_y[l]), col="white")
  for (l in 1:length(at_x)) lines(c(at_x[l], at_x[l]), c(at_y[1]-10,at_y[length(at_y)]+10), col="white")
  
  ## x axis
  axis(1, label=at_x, at=at_x, lwd=0, lwd.ticks=1)
  axis(2, label=at_y, at=at_y, lwd=0, lwd.ticks=1)
  
  ## add the label above the plot
  x <- par("usr")[2] - 0.2  # Adjust the x-coordinate as needed
  y <- par("usr")[4] + 0.05  # Adjust the y-coordinate to position it above the plot
  mtext(text = label, side = 3, line = 1, at = x, cex = 1.25)
  
}

## add_axis_and_grid_short
add_grid = function()
{
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray95", border=NA) # Add gray background
  abline(v = axTicks(1), h = axTicks(2), col="white", lwd=1.5) # Add white gridlines at tick positions
}


#
###