#
#      A stochastic cellular automaton model
#          of vegetation cover
#
# --------------------------------------------------
#
# The MIT License (MIT)
# 
# Copyright (c) 2014 Florian D. Schneider
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

source("simfunctions.r")

# specify lattice
width = 50
height = 50

mapping(width, height)

# initial cell states
states = c("+","0","-")
prob = c(1/10,1/10,8/10)

# time and resolution of simulation
timesteps = 400
delta = .1

# defining parameter set
parameters = list(
  m = 0.1,   	# intrinsic mortality
  b = 0.4, 		# beta*eps 
  d = 0.2,		# degradation
  c_ = 0.3, 		# beta*g  
  del = 0.0,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
  r = 0.0, 	# regeneration rate
  f = 0.8 		# local fascilitation
)

parms_temp <- parameters # copying parms for the simulation

# sampling the initial random grid into a list object
parms_temp$rho_plus <- 0.8

# how many initial plants
init_plant <- as.integer(width*height*parms_temp$rho_plus)

# vector of empty cells in state "0"
cells <- factor(rep("0", times = width*height), levels = states) # 1:length(states))
# replace init_plant cells, randomly drawn, with "+"
cells[sample(1:(width*height), init_plant, replace = FALSE)] <- "+"

# which cells are still empty?
empty <- which(cells != "+")
# select which will be degraded? fixed to 50% of non occupied cells. 
init_degraded <- sample(empty, length(empty) * (parms_temp$r + parms_temp$f/2) )  
# replace cell state. 
cells[init_degraded] <- "-"

initial <- list(  
  dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
  cells = cells#contains a random row-wise, factorial vector to fill the grid 
)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

color <- c("black","grey80", "white") # define colors for the cell state levels

#plot the initial state
#par(mar = c(0,0,0,0), mfrow = c(1,1))
#plot(initial, cols = color)


# initialising result list 
result <- list()   # create an output list and
result$time <- seq(0, timesteps, delta) # write simulated timesteps

#save all level global densities: rho_*
result$rho_plus  <- vector("numeric", length = length(result$time))
result$rho_plus[1] <- parms_temp$rho_plus

result$timeseries <- list()
result$timeseries[1:timesteps] <- rep(list(initial), times = n_snaps) 

#  initialising first iteration 
x_old <- initial    #old landscape for first timestep   

for(i in 1:(timesteps/delta)+1) {    #calculation loop
  
  x_new <- x_old 		# copy x_old into an object x_new to allocate memory
  
  # model specific part:

# 2 - drawing random numbers
rnum <- runif(width*height) # one random number between 0 and 1 for each cell

# 4 - setting probabilities for transitions in this timestep

  # count local density of occupied fields for each cell: 
  parms_temp$Q_plus <- count(x_old, "+")/4
  

  # calculate recolonisation rates of all cells
  recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*delta)
  
  # calculate death rates
  death <- with(parms_temp, m*delta)
 
  # correct for overshooting death prob
  #death[death > 1] <- 1
  
  regeneration <- with(parms_temp, (r + f*Q_plus)*delta)

degradation <- with(parms_temp, (d *delta))

# check for sum of probabilities to be inferior 1 and superior 0
if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 

  # 5 - applying the rules to fill the cells in x_new

  x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
  x_new$cells[which(x_old$cells == "+"  & rnum <= death)] <- "0"

x_new$cells[which(x_old$cells == "0"  & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"  
x_new$cells[which(x_old$cells == "-"   & rnum <= regeneration)] <- "0"

  
  # 6 saving state of the new grid		

parms_temp$rho_plus <- sum(x_new$cells == "+")/(width*height) 
result$rho_plus[i] <- parms_temp$rho_plus

  # activate to save each single timeseries step
  if(result$time[i] %in% 1:timesteps) result$timeseries[[as.integer(result$time[i])]] <- x_new  #the whole grid is saved to timeseries
  # activate to save each single iteration
  #result$timeseries[[i]] <- x_new  #the whole grid is saved to timeseries

  # activate for plotting during simulation (very slow !!!)
  #if(i %in% seq(1,(timesteps/delta)+1, 10)) plot(x_new, col = color, add = TRUE)
  
  x_old <- x_new 
  gc()  #garbage collection
} # end of simulation.


plot(result$rho_plus ~ result$time, type = "l", ylim = c(0,1))

# FIGURE 3 -- animated gif (in Windows: ImageMagick is required)
library(animation)
if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 

saveGIF( 
  for(i in seq(1, length(result$timeseries), 1)) {
    par(mar = c(0,0,0,0))
    plot(result$timeseries[[i]], grid = FALSE, cols = color, ani = TRUE)
  }
  , movie.name = "facilitation.gif", img.name = "grid", convert = "convert", interval = 0.01/1,
  cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())

