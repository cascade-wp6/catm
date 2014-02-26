#################################################
#
#		model: iterator for cellular automata
#
#
#author: Flo
#date: 02.04.2013
#
#
#
#################################################



rm(list=ls())
source("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\automatafunctions.r")
setwd("C:\\Users\\SCHNEIDER\\Dropbox\\Projects\\CASCADE\\cellularautomata\\code v1.0\\thread herbivory\\")
library(foreach)
library(doSNOW)
winWorker <- list(host = "localhost", 
				rscript = "C:\\R\\R-2.15.3\\bin\\x64\\Rscript.exe", 
				snowlib = "C:\\R\\R-2.15.3\\library"
	)
winWorker <- list(host = "localhost", 
				rscript = "C:/R/R-2.15.3/bin/x64/Rscript.exe", 
				snowlib = "C:/R/R-2.15.3/library"
	)
ubuWorker <- list(host = "kefi118", user = "schneider",
				rscript = "/usr/lib/R/bin/Rscript", #R/bin/
				snowlib = "/usr/lib/R/site-library/"
	)
#	, rep(list(ubuWorker), times = 23)
workerlist <- rep(list(winWorker), times = 11)
cl <- makeSOCKcluster(workerlist)
registerDoSNOW(cl)
#foreach(i = 1:10) %dopar% {
#runif(10^7)
#i
#}
#stopCluster(cl)


# specify lattice
width = 50  
height = 50


# time and resolution of simulation
timesteps = 800
#addgrazing = 300
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


# defining parameter set
# Effect of mortality (m). d = 0.1, r = 0.01, c = 0.2, d = 0.1, f = 0.9. (Kefi et al 2007 TPB Fig 4c)
parameters = list(
	m = 0.1, 		# intrinsic mortality
	b = 0.5, 		# beta*eps 
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9 		# local fascilitation
)
	
# initial cell states
states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels
	
env <- 	seq(0.2,1,length = 200)
mort <- seq(0.1, 0.3, length = 3)

lmort <- length(mort)
lenv <- length(env)

header_output <- data.frame(
			ID = NA,
			mortality = NA, 
			environmt = NA,  
			rho_plus = NA, 
			rho_plus_ungrazed = NA, 
			rho_minus = NA,
			q_plus = NA)
			

write.table(header_output[-1,], "output_Kefi.csv", row.names = FALSE, sep = ",")

write.table(matrix(NA, ncol = 2+(width*height)), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",")

env <- rep(rep(env, times = lmort), times = 2)
mort <- rep(rep(mort, each = lenv), times = 2)



foreach(iteration = 1:(lenv*lmort)) %dopar% {

parameters$b <- env[iteration]
parameters$m <- mort[iteration]

flag = "low"

while(flag %in% c("low", "high")) {

set.seed(iteration)
if(flag == "low") prob = c(.1/10,.9/10,9/10)
if(flag == "high") prob = c(9/10,.9/10,0.1/10)

# sampling the initial random grid into a list object
initial <- list(  
	dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
	cells = sample(factor(1:length(states)), width*height, replace = T, prob = prob ) #second element contains a random row-wise, factorial vector to fill the grid 
	)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

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
			result$q_[[j]][1]  <- mean(rowSums( sapply(interact, 	function(k) (initial$cells == states[j])[x_with_border][x_to_evaluate[initial$cells == states[j]]+k]) )) /4# write initial local densities 
			}
			
			#result$timeseries <- list()
			#result$timeseries$initial <- initial
			
#write.table(matrix(c(iteration, flag, initial$cells), nrow = 1), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
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
		parms_temp$Q_unveg <- (count(x_old, "0")+count(x_old, "-"))/4 
		
# 2 - drawing random numbers
		rnum <- runif(width*height) # one random number between 0 and 1 for each cell
	
# 4 - applying the rules to fill the cells in x_new
		
		recolonisation <- with(parms_temp, (del*rho+(1-del)*Q_plus)*(b-c_*rho)*delta)
		degradation <- with(parms_temp, (d *delta))
		death <- with(parms_temp, m*delta)
		regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
		
		if(any(c(recolonisation+degradation, death, regeneration) > .95 )) warning(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
				
		if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
		
		
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
		
		# activate to save timeseries at activation of herbivory
		#if(i == addgrazing/delta) result$timeseries$stable_without_grazing <- x_new  #the whole grid is saved to timeseries
		
		# activate for plotting during simulation (very slow !!!)
		#if(i %in% seq(1,(timesteps/delta)+1, 10)) plot(x_new, col = color, add = TRUE)
		#if(i %in% seq(1,timesteps, 50)/delta) write.table(c(iteration, i*delta, flag, x_new$cells), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

		x_old <- x_new
		
		gc()  #garbage collection
	} # end of simulation.

	#result$timeseries$final <- x_new
	
#write.table(matrix(c(iteration, flag, x_new$cells), nrow = 1), "grids.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
	
gc()
	t_eval <- ((timesteps/delta)-50/delta):(timesteps/delta)+1
	
	out <- data.frame(ID = iteration,
			mortality = parameters$m, 
			environmt = parameters$b,  
			rho_plus = mean(result$rho[[1]][t_eval]), 
			rho_plus_ungrazed = 0, #mean(result$rho[[1]][(t_eval-1)-(timesteps-addgrazing)/delta]), 
			rho_minus = mean(result$rho[[2]][t_eval]),
			q_plus = mean(result$q_[[2]][t_eval]))

write.table(out, "output_Kefi.csv", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)

if(flag == "high") {flag = "ok"}
if(flag == "low" && out$rho_plus < 10e-3) {flag = "high"} else {flag = "ok"}


#}

}


} 

output <- read.csv("output_Kefi.csv", header = TRUE)   


pdf("herbivory_stability.pdf", height = 8, width = 4, paper = "special")
par(mfrow = c(3,1), mar = c(4,4,1,1))
with(data = output[output$mortality == 0.1,], plot(environmt, rho_plus, xlim = c(0,1), ylim = c(0,1), pch = 20))
with(data = output[output$mortality == 0.2,], plot(environmt, rho_plus, xlim = c(0,1), ylim = c(0,1), pch = 20))
with(data = output[output$mortality == 0.3,], plot(environmt, rho_plus, xlim = c(0,1), ylim = c(0,1), pch = 20))
dev.off()

grids <- read.csv("grids.csv", header = FALSE, sep = ",")  
grids[,c(1:2)]

class(grids[,c(3:5)])
show <- list(dim = c(width, height),
			cells = as.factor(as.integer(grids[323,]))
			)
levels(show$cells) <- states
class(show) <- c("list", "landscape")
plot(show, col = color)
par(mfrow = c(1,3), mar = c(1,1,1,1))
plot(saveall[[n]]$initial, col = color)
plot(saveall[[n]]$stable_without_grazing, col = color)
plot(saveall[[n]]$final_with_grazing, col = color)


stopCluster(cl)
