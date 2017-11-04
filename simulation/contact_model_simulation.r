
rm(list=ls(all=TRUE))

require(ggplot2)
require(reshape2)
require(boot)
require(geosphere)
require(ggmap)
require(stpp)
require(gridExtra)
require(animation)
require(fields)
require(rgdal)
require(raster)
require(maptools)
require(rgeos)
require(SDMTools)
require(spatstat)
require(ecospat)

set.seed(1)


## Key data structures ##
#########################

# 1) pop_grid = a matrix whose elements are the population densities of corresponding grids
# 2) grid_lines = a dataframe stores the straight lines that form the grids;  the numbering (k_line): bottom to top, then left to right; orientation of line, 1 for horizontal, 2 for vertical; also, left endpoints (coor_x_1,coor_y_1) and right endpoints (coor_x_2,coor_y_2) that define the lines
# 3) simulated_epi = a dataframe to store simulated epidemic, with column names k (kth case, counted from k=0), coor_x, coor_y, t_e, t_i, t_r, age, infected_source (i.e. source of infection)

path <- "xxx" # working dir

## !! need to read from files !! ##
## attached text files are taken from Western Area Sierra Leone ##
pop_grid <- as.matrix(read.csv(file=paste(path,"pop_grid.txt",sep=""),header=FALSE))
grid_lines <- read.csv(file=paste(path,"grid_lines.txt",sep=""),stringsAsFactors=FALSE)
###

source(file=paste(path,"contact_model_simulation_functions.r",sep=""))



grid_size <- 0.1 #in km width of the square grid

age_level <- c(0,1) # levels of age gp
age_dist  <- c(0.5,0.5) # proportion of each age gp

m_start  <- 1 # number of initial cases 

t_max <- 30 #  observation duration
t_intervention <- 10 # day of intervention

## parameter values ##
alpha <- 0.1
beta_1 <- 2
beta_2 <- 1.5
k_1 <- 0.4
mu_lat <- 5
var_lat <- 1
c <- 4
omega <- 0.2

##

min_coor_x <- min(c(grid_lines$coor_x_1,grid_lines$coor_x_2)) # bounds of the region 
max_coor_x <- max(c(grid_lines$coor_x_1,grid_lines$coor_x_2))
min_coor_y <- min(c(grid_lines$coor_y_1,grid_lines$coor_y_2))
max_coor_y <- max(c(grid_lines$coor_y_1,grid_lines$coor_y_2))

x_intervals <- seq(min_coor_x,max_coor_x,grid_size)
y_intervals <- seq(min_coor_y,max_coor_y,grid_size)

n_line <- max(grid_lines$k) + 1 # number of lines
n_row_grid <- nrow(pop_grid) # number of rows of grids
n_col_grid <- ncol(pop_grid)  # number of cols of grids

###
simulated_epi <- data.frame(k=numeric(0), coor_x=numeric(0), coor_y=numeric(0), t_e=numeric(0), t_i=numeric(0), t_r=numeric(0), age=numeric(0), infected_source=numeric(0))

## initialize the index cases ##
index_k <- 0:(m_start-1)
index_coor_x <- runif(m_start,min_coor_x,max_coor_x)
index_coor_y <- runif(m_start,min_coor_y,max_coor_y)
index_t_e <- rep(0,m_start) # all assumed to infected at time=0
index_t_i <- index_t_e + rgamma(m_start,shape=mu_lat*mu_lat/var_lat, scale=var_lat/mu_lat)
index_t_r <- index_t_i + rexp(m_start,rate=1/c)
index_age <- sample(age_level,size=m_start,prob=age_dist,replace=T)
index_source <- rep(9999,m_start) # all assumed to be infected by the background

simulated_epi[1:m_start,] <- c(index_k, index_coor_x,index_coor_y,index_t_e,index_t_i,index_t_r,index_age,index_source)
#####

t_now <- min(simulated_epi$t_i)

simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious

num_infection <- nrow(simulated_epi)


t_next <- t_now # to start the while loop
while(t_next<(t_max)){

	### simulate the timings, and the source, of next infection ###

	simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious

	if(nrow(simulated_epi_sub)>=1){
		beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))
		total_beta <- sum(beta_infectious)
		joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
		min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])
	}

	if(nrow(simulated_epi_sub)<1){
		total_beta <- 0
	}


	t_next <- simulate_NHPP_next_event (t_now=t_now, t_intervention=t_intervention, sum_beta=total_beta, alpha=alpha, omega=omega, t_max=t_max) # simulate the next infection time using thinning algorithm

	cat("t_next =",t_next, "\n")

	while(t_next>=min_I_R & t_next!=Inf){

		t_now <- min_I_R
		simulated_epi_sub <- subset(simulated_epi, simulated_epi$t_i<=t_now & simulated_epi$t_r>t_now) # those are currently infectious


		if(nrow(simulated_epi_sub)>=1){

			beta_infectious <- sapply(simulated_epi_sub$age, FUN=beta_by_age, c(beta_1,beta_2))

			total_beta <- sum(beta_infectious)

			joint_I_R <- c(simulated_epi_sub$t_i,simulated_epi_sub$t_r)
			min_I_R <- min(joint_I_R[which(joint_I_R>t_now)])
			t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, alpha=alpha, omega=omega, t_max=t_max)

		}

		if(nrow(simulated_epi_sub)<1){
			total_beta <- 0
			t_next <- simulate_NHPP_next_event (t_now=t_now,  t_intervention=t_intervention, sum_beta=total_beta, alpha=alpha, omega=omega, t_max=t_max)
			break # needed when consider background infection
		}

	} # end of while(t_next>=min_I_R)


	k <- num_infection + 1 - 1 # k=0,1,2...
	t_now <- t_next
	t_i_new <- t_now +  rgamma(1,shape=mu_lat*mu_lat/var_lat, scale=var_lat/mu_lat)
	t_r_new <- t_i_new + rexp(1,rate=1/c)


	if(nrow(simulated_epi_sub)>=1) source <- sample(c(9999,simulated_epi_sub$k),size=1, prob=c(alpha,beta_infectious)) # 9999 = background

	if(nrow(simulated_epi_sub)<1) source <- 9999 # 9999 = background

	age <- sample(age_level,size=1,prob=age_dist)

	### simulate the coordinates of the new infection (above) ###

	x_new = min_coor_x - 5 # to start the while loop
	y_new = min_coor_y - 5 # to start the while loop
	while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){

			if (source!=9999){

					r <- rexp(1,rate=k_1)		

					set_points <- circle_line_intersections (circle_x=simulated_epi$coor_x[source+1],circle_y=simulated_epi$coor_y[source+1], r,  n_line=n_line, grid_lines=grid_lines)

					n_set_points = nrow(set_points)

					if (n_set_points>=1) {
						
						arcs <- func_arcs_attributes(set_points, pop_grid, r, min_coor_x, min_coor_y, grid_size, n_row_grid, n_col_grid) 
						arcs$mass <- arcs$den*arcs$len
						arcs <- arcs[order(arcs$theta_abs),]

						sum_arcs_den <- sum(arcs$den)

						if (sum_arcs_den>0){
							k_segment <- sample(1:nrow(arcs),size=1,prob=arcs$mass) # decide which segment the new infection would lie
							theta_within_segment <- runif(1, min=0, max=arcs$theta_abs[k_segment]) # uniformly draws a point within the segment chosen above

							if (k_segment==1) theta_from_y_eq_0 <- theta_within_segment # the measure of theta from y=0
							if (k_segment!=1) theta_from_y_eq_0 <- sum(arcs$theta_abs[1:(k_segment-1)]) + theta_within_segment # the measure of theta from y=0

							x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
							y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection


						}

						###
						if (sum_arcs_den==0){
							theta_from_y_eq_0 <-  runif(1, min=0,max=(2*pi)) # draw theta uniformly between [0,2pi]
							x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
							y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection

				
						}
						###


					}


					if (n_set_points<1) {
						theta_from_y_eq_0 <-  runif(1, min=0,max=(2*pi)) # draw theta uniformly between [0,2pi]

						x_new <- simulated_epi$coor_x[source+1] + r*cos(theta_from_y_eq_0) # the coor_x of the new infection
						y_new <- simulated_epi$coor_y[source+1] + r*sin(theta_from_y_eq_0) # the coor_x of the new infection



					}


			}

			if(source==9999){

				k_grid <- sample(1:length(as.numeric(pop_grid)),size=1, prob=as.numeric(pop_grid))
				m_grid <- k_grid%%nrow(pop_grid) # the mth row of the grids
				if(m_grid==0) m_grid <- nrow(pop_grid)			
				n_grid <- ceiling(k_grid/nrow(pop_grid))  # nth column ..
				x_new <- runif(1,min=x_intervals[n_grid],max=x_intervals[n_grid+1]) # random lands at a point in the grid selected
				y_new <- runif(1,min=y_intervals[m_grid],max=y_intervals[m_grid+1])

							
			}



	} # end while(x_new<min_coor_x | x_new>max_coor_x | y_new<min_coor_y | y_new>max_coor_y){


	simulated_epi[k+1,] <- c(k, x_new, y_new, t_now, t_i_new, t_r_new, age, source) 
	####

	num_infection <- num_infection + 1

} # end of while(t_next<t_max)




simulated_epi


ggplot(data=simulated_epi) + geom_point(aes(x=coor_x,y=coor_y),size=0.4) 

