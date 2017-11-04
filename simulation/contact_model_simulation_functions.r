


### functional form of the intensity (beta) as a function of t ###

func_time_beta <- function(t, t_intervention, sum_beta, alpha, omega){

	if (t>=t_intervention) beta_t <- (sum_beta+alpha)*exp(-1*(omega)*(t- t_intervention))
	if(t<t_intervention) beta_t <- sum_beta + alpha

# beta_t <- sum_beta + alpha

beta_t


}

### simulate the next arrival time of a non-homogeneous Poisson process  ###
simulate_NHPP_next_event <- function(t_now, t_intervention, sum_beta, alpha, omega,t_max) {

u=1
acp_pr=0

total_rate_now <-  func_time_beta(t_now, t_intervention, sum_beta, alpha, omega)


if ((total_rate_now)>0) {
# if ((sum_beta+alpha)>0) {

	while(u>acp_pr){
		
		t_sim <- t_now + rexp(1,rate=total_rate_now)
		acp_pr <- func_time_beta(t_sim, t_intervention, sum_beta, alpha, omega)/(total_rate_now)

		# t_sim <- t_now + rexp(n=1,rate=(sum_beta+alpha))
		# acp_pr <- func_time_beta(t_sim, t_intervention, sum_beta, alpha, omega)/(sum_beta+alpha)

		u <- runif(1)

		if(u<=acp_pr) t_next <- t_sim

	}


    if (t_sim<=t_max){
        t_next = t_sim
    } else {t_next=Inf}


}

if ((total_rate_now)<=0) {
# if ((sum_beta+alpha)<=0) {
	t_next <- Inf
}

t_next

}


### assign the proper beta for an individual by age ###
beta_by_age <- function(age, beta_by_age_vector){

	beta_by_age_vector[age+1]

}

### pow(a,b) in C++ ###
pow <- function(a,b){
	a^b
}

### return a set of intersection points between the cirlce and the grid lines ###
circle_line_intersections <- function (circle_x, circle_y, r,  n_line, grid_lines){ 
## note that formula assume the circle sits at (0,0), so need a shift of coordinates before and after##

# set_points = data.frame(coor_x=numeric(0),coor_y=numeric(0),theta=numeric(0))
set_points = matrix(NA, ncol=3)


x_lb = -r;
x_ub = r;
y_lb = -r;
y_ub = r;


for (i in 1:n_line){

    if ( (grid_lines$orient_line[i]==1 & (grid_lines$coor_y_1[i]- circle_y)>=y_lb & (grid_lines$coor_y_1[i]- circle_y)<=y_ub) | (grid_lines$orient_line[i]==2 & (grid_lines$coor_x_1[i]-circle_x)>=x_lb & (grid_lines$coor_x_1[i]-circle_x)<=x_ub) ) { # only need to consider those grids line intersect with the smallest bounding box that contains the circle

       	dx = grid_lines$coor_x_2[i] - grid_lines$coor_x_1[i];
        dy = grid_lines$coor_y_2[i] - grid_lines$coor_y_1[i];
        dr = sqrt(pow(dx,2.0)+pow(dy,2.0));
        D = (grid_lines$coor_x_1[i]-circle_x)*(grid_lines$coor_y_2[i]-circle_y) - (grid_lines$coor_x_2[i]-circle_x)*(grid_lines$coor_y_1[i]-circle_y);

        sgn_dy=1;
        if (dy<0) sgn_dy = -1;

        Delta = pow(r,2.0)*pow(dr,2.0) - pow(D,2.0);

        # switch(Delta>=0){
        #     case 1:{ // tangent (one intersection) or two intersection



                if (Delta>0){

                    x_1 = (D*dy + sgn_dy*dx*sqrt(Delta))/pow(dr,2.0) + circle_x;
                    x_2 = (D*dy - sgn_dy*dx*sqrt(Delta))/pow(dr,2.0) + circle_x;

                    y_1 = (-D*dx + abs(dy)*sqrt(Delta))/pow(dr,2.0) + circle_y;
                    y_2 = (-D*dx - abs(dy)*sqrt(Delta))/pow(dr,2.0) + circle_y;

                    theta_1 = atan2((y_1-circle_y),(x_1-circle_x));
                    theta_2 = atan2((y_2-circle_y),(x_2-circle_x));

                    if (theta_1<0) theta_1 = 2*pi + theta_1;
                    if (theta_2<0) theta_2 = 2*pi + theta_2;

       
                    # set_points_struct tmp_1, tmp_2;
                    # tmp_1.coor_x= x_1;
                    # tmp_1.coor_y=y_1;
                    # tmp_2.coor_x=x_2;
                    # tmp_2.coor_y=y_2;
                    # tmp_1.theta=theta_1;
                    # tmp_2.theta=theta_2;
                    # set_points.push_back(tmp_1);
                    # set_points.push_back(tmp_2);


                    set_points = rbind( set_points, c(x_1,y_1,theta_1));
                    set_points = rbind( set_points, c(x_2,y_2,theta_2));

  
                }


                if (Delta==0){
                    x = D*dy/pow(dr,2.0)+ circle_x;

                  	y = -D*dx/pow(dr,2.0) + circle_y;

                   	theta = atan2((y-circle_y),(x-circle_x));


                    if (theta<0) theta = 2*pi + theta;

    
                    # set_points_struct tmp;           
                    # tmp.coor_x=(x);
                    # tmp.coor_y=(y);    
                    # tmp.theta=(theta);
                    # set_points.push_back(tmp);

                    set_points = rbind( set_points, c(x,y,theta));


                }

            # break;
            # }

            # case 0:{ // no intersection
            #     // do nothing
            # break;
            # }

        # }


    
    }



}

# std::sort(set_points.begin(), set_points.end(), by_theta());

set_points = set_points[-1,];

colnames(set_points) = c("coor_x", "coor_y", "theta");
set_points = data.frame(set_points);


set_points =  set_points[order(set_points$theta),];


set_points;

}



### return a set of segments (with length, density, and absolute angle) correspond to a set_points ##
func_arcs_attributes <- function (set_points, pop_grid,  r, x_min, y_min, grid_size, n_row_grid, n_col_grid){ 

    n_segments = nrow(set_points);

    segments = data.frame(theta_abs=rep(NA,n_segments),len=rep(NA,n_segments),den=rep(NA,n_segments),m=rep(NA,n_segments),n=rep(NA,n_segments));

    # segments.resize(n_segments);

    for (i in 1:n_segments){
        if(i==1) {

            segments$theta_abs[i] = set_points$theta[1] + (2*pi - set_points$theta[n_segments]);

            x_1= set_points$coor_x[n_segments];
            y_1= set_points$coor_y[n_segments];
           	x_2= set_points$coor_x[1];
            y_2= set_points$coor_y[1];
            midpoint_x = (x_1+x_2)/2.0;
            midpoint_y = (y_1+y_2)/2.0;
            m = ceiling((midpoint_y - y_min)/grid_size); #at mth row of the grid
            n = ceiling((midpoint_x - x_min)/grid_size); # at nth col..

            segments$m[i] = m;
            segments$n[i] = n;


            # switch(m>para_other.n_row_grid | n>para_other.n_col_grid | m<=0| n<=0 ){
            #     case 1:{
            #          segments[i].den = 0.0;
            #     break;
            #     }
            #     case 0:{

            #         segments[i].den = pop_grid[m-1][n-1];
            #     break;
            #     }
            # }

            if (m>n_row_grid | n>n_col_grid | m<=1| n<=1 ){
                     segments$den[i] = 0.0;
            } else {
                    segments$den[i] = pop_grid[m,n];   
            }


        }

        if(i!=1) {

            segments$theta_abs[i] = set_points$theta[i] - set_points$theta[i-1];

            x_1= set_points$coor_x[i-1];
            y_1= set_points$coor_y[i-1];
            x_2= set_points$coor_x[i];
            y_2= set_points$coor_y[i];
            midpoint_x = (x_1+x_2)/2.0;
            midpoint_y = (y_1+y_2)/2.0;
            m = ceiling((midpoint_y - y_min)/grid_size); # at mth row of the grid
            n = ceiling((midpoint_x - x_min)/grid_size); # at nth col..

            segments$m[i] = m;
            segments$n[i] = n;


            # switch(m>para_other.n_row_grid | n>para_other.n_col_grid | m<=0| n<=0){
            #     case 1:{
            #          segments[i].den = 0.0;
            #     break;
            #     }
            #     case 0:{
            #         segments[i].den = pop_grid[m-1][n-1];
            #     break;
            #     }
            # }


            if (m>n_row_grid | n>n_col_grid | m<=1| n<=1 ){
                     segments$den[i] = 0.0;
            } else {
                    segments$den[i] = pop_grid[m,n];   
            }

        }


        segments$len[i] = r*(segments$theta_abs[i]);  



    }



segments;

}

