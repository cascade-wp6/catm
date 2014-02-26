#################################################
#
#		model: facilitation model (Kefi et al 2007)
#
#
#author: Flo
#date: 02.04.2013
#
#
#################################################


rm(list=ls())

# plotting function for objects of class "landscape".
plot.landscape <- function(x, grid = FALSE, axis = FALSE, cols = "auto", add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") cols = greyscale(nlev) # default color value
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, x$dim[1]+0.5+adj), ylim = c( x$dim[2]+0.5+adj, 0+0.5+adj), bty = c("n", "o")[grid+1], xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", ... ) 
  
  if(axis && !add) axis(3) 
  if(axis && !add) axis(2)
  
  if(grid) border = "grey80" else border = cols[as.numeric(x$cells)]
  
  rect(rep(1:x$dim[1], times = x$dim[2])-.5, rep(1:x$dim[2], each = x$dim[1])-.5, rep(1:x$dim[1], times = x$dim[2])+.5, rep(1:x$dim[2], each = x$dim[1])+.5, col = cols[as.numeric(x$cells)], border = border)
  
  if(grid) box()
 }

count  <- function(x, neighbor) {
			neighbors <- numeric(length = prod(x$dim))
			x_logical_with_border <- (x$cells == neighbor)[x_with_border]
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
			}
			return(neighbors)	
			}
	
# function which counts the density of a state within neighboring cells of a particular state
# requires a previous global definition of the transformation vectors "x_to_evaluate" and "x_with_border" and the neighborhood vector "interact"
localdens <- function(x, neighbor, state) {
			neighbors <- numeric(length = length(which(x$cells == state))) # preallocate memory
			
			# create a logical vector x(which cells are in state "state"?) 
			# and transfer the vector x to a vector x_with_border: x[x_with_border]
			x_logical_with_border <- (x$cells == neighbor)[x_with_border]
			
			#of those values which are original values (x_to_evaluate) and which are "state": x_to_evaluate[initial$cells == "state"]
			x_state <- x_to_evaluate[x$cells == state]
			
			for(k in interact) {
				neighbors <- neighbors + x_logical_with_border[x_state+k]
			}
			# derive a vector which contains the value to count (also "state") in one of the neighbor cells j
			# cumulate the result for each of the four neighboring cells (for(k in interact))

			mean(neighbors/4)	
			# devide by 4 to get density
			# average over all cells
			}
			

# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
prob = c(1/10,1/10,8/10)

# time and resolution of simulation
timesteps = 400
delta = 1/10



# derive helper vectors for counting: 
# transformation vector for evaluation at the border of the grid
	# set evaluation matrix 
	X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
	# setting the border of the evaluation matrix X
	X <- cbind(X[,width], X, X[,1] )  
	X <- rbind(X[height,], X, X[1,] ) 
	# transformation vector which adds the border to the lattice:
	x_with_border <- as.integer(t(X))
	
	# from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
	x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]	)		

# defining the neighborhood which is to be evaluated	
	# set interaction matrix
	I <- matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)	
	# coordinates of neighbours in Interaction matrix I: 
	neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
	# coordinates relative to the evaluated cell (=  which(is.na(I) ) 
	relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
	relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
	
	# relative position of the four direct neighbours of a cell
	interact <- (relrow * dim(X)[2] + relcol)

	

# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

color <- c("black","grey80", "white") # define colors for the cell state levels

#plot the initial state
par(mar = c(0,0,0,0), mfrow = c(1,1))
plot(initial, cols = color)

# defining parameter set
parameters = list(
	m = 0.2, 		# intrinsic mortality
	b = 0.62, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.0, 	# regeneration rate
	f = 0.9 		# local fascilitation
)


# initialising result list 
result <- list()   # create an output list and
			result$time <- seq(0, timesteps, delta) # write simulated timesteps
			
			#save all level global densities: rho_*
			result$rho  <- list()			
			for(i in 1:length(states)) {
				result$rho[[i]] <- numeric(length = timesteps/delta)
				result$rho[[i]][1]  <- sum(initial$cells == states[i])/(width*height) # write initial rho 
			}
			
			result$q_  <- list()			
			for(j in 1:length(states)) {
			result$q_[[j]] <- numeric(length = timesteps/delta)
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, 	function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == states[j]]+k]) )) /4# write initial rho 
			}
			
			result$timeseries <- list() # create a subordinate list and
			result$timeseries[[1]] <- initial # write 'x' as the first entry of the list
			for(i in 1:(timesteps/delta)+1) result$timeseries[[i]] <- initial	# allocate memory for each timeseries object
			

#  initialising first iteration 
	x_old <- initial    #old landscape for first timestep   
	parms_temp <- parameters # copying parms for the simulation and multiplying
	
for(i in 1:(timesteps/delta)+1) {    #calculation loop

		x_new <- x_old 		# copy x_old into an object x_new to allocate memory

# model specific part:
# 1 - setting time-step parameters

		parms_temp$rho <- result$rho[[1]][i-1]   # get old rho plus for this timestep 
		
	 # count local density of occupied fields for each cell: 
 		parms_temp$Q_plus <- count(x_old, "+")/4   
		# This is a vectorised evaluation: first transform the original grid x_old, into a grid with border. Here check for all cells of the original grid (x_to_evaluate) the content of all four neighbouring cells, whose relative position is defined in interact. The number of cells with content "+" are count and divided by 4.
	# count local density of degraded fields for each cell
		#parms_temp$Q_minus <- rowSums(sapply(interact, function(j)  (x_old$cells == "-")[x_with_border][x_to_evaluate+j]))/4 
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		death <- with(parms_temp, m*delta)
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
		x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"
		x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
		x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"


# 5 saving state of the new grid		
		for(j in 1:length(states)) {
			result$rho[[j]][i]  <- sum(x_new$cells == states[j])/(width*height) # write rhovalues 
			}

		
		for(j in 1:length(states)) {
		result$q_[[j]][i]  <- localdens(x_new, states[j], "+")
		}
		
		# activate to save each single timeseries step
		if(i %in% seq(1,(timesteps/delta)+1, 10)) result$timeseries[[i]] <- x_new  #the whole grid is saved to timeseries
		
		# activate for plotting during simulation (very slow !!!)
		#if(i %in% seq(1,(timesteps/delta)+1, 10)) plot(x_new, col = color, add = TRUE)

		x_old <- x_new 
		gc()  #garbage collection
	} # end of simulation.


#str(result)

	
	
	
	
# the rest is graphical output

# FIGURE 1 	-- final states of the grid
par(mar = c(0,0,0,0))
plot(x_new, cols = color)


# FIGURE 2 	-- global densities and local densities around cells in state "+" 
par(mfrow = c(2,1))
par(mar = c(4,4,2,1)+.1, xaxt = "s", yaxt ="s")
plot(NA,NA, ylim = c(0, 1), xlim = range(result$time), type ="l", xlab = "time", ylab = expression(rho ['*']), xaxs ="i", yaxs = "i" )
lines(result$time, result$rho[[1]])
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]]+result$rho[[3]], 0), col = "white")
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]]+result$rho[[2]], 0), col = "grey50")
polygon(c(0,result$time,max(result$time)), c(0, result$rho[[1]], 0), col = "black")
	
par(mar = c(4,4,2,1)+.1, xaxt = "s", yaxt ="s")
plot(NA,NA, ylim = c(0, 1), xlim = range(result$time), type ="l", xlab = "time", ylab = "q*", xaxs ="i", yaxs = "i" )
lines(result$time, result$q_[[1]])
polygon(c(0,result$time,max(result$time)), c(0, result$q_[[1]]+result$q_[[2]]+result$q_[[3]], 0), col = "white")
polygon(c(0,result$time,max(result$time)), c(0, result$q_[[1]]+result$q_[[2]], 0), col = "grey50")
polygon(c(0,result$time,max(result$time)), c(0, result$q_[[1]], 0), col = "black")
par(mfrow = c(1,1))


# FIGURE 3 -- animated gif (in Windows: ImageMagick is required)
library(animation)
if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 

saveGIF( 
	for(i in seq(1, length(result$timeseries), 1/delta)) {
	    par(mar = c(0,0,0,0))
		plot(result$timeseries[[i]], grid = FALSE, cols = color, ani = TRUE)
	}
, movie.name = "facilitation_test.gif", img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
