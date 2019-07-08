
#include "contact_model_extension_sim_inference_header.h"

//---
double func_time_beta (const double& beta, const double& t, const double& T,  const double& omega){


    double beta_t;

    switch(t<T){
        case 1:{
            beta_t =  beta;
        break;
        }

        case 0:{
            beta_t = beta*exp(-omega*(t-T));
        break;
        }

    }



return(beta_t);

}

//---
double func_time_alpha (const double& alpha, const double& t, const double& T,  const double& omega){


    double alpha_t;

    switch(t<T){
        case 1:{
            alpha_t =  alpha;
        break;
        }

        case 0:{
            alpha_t = alpha*exp(-omega*(t-T));
        break;
        }

    }



return(alpha_t);

}

//------
vector<set_points_struct> circle_line_intersections (double circle_x, double circle_y, double r, int& n_line, grid_lines_struct& grid_lines){ // return a set of intersection points between the cirlce and the grid lines

// note that formula assume the circle sits at (0,0), so need a shift of coordinates before and after//

vector<set_points_struct>  set_points;
set_points.reserve(n_line*2);

// set_points_struct set_points;
// set_points.coor_x.reserve(n_line*2);// each line has max 2  intersections with the circle
// set_points.coor_y.reserve(n_line*2);
// set_points.theta.reserve(n_line*2);

double x_lb = -r;
double x_ub = r;
double y_lb = -r;
double y_ub = r;


for (int i=0;i<=(n_line-1);i++){

    if ( (grid_lines.orient_line[i]==1 & (grid_lines.coor_y_1[i]- circle_y)>=y_lb & (grid_lines.coor_y_1[i]- circle_y)<=y_ub) | (grid_lines.orient_line[i]==2 & (grid_lines.coor_x_1[i]-circle_x)>=x_lb & (grid_lines.coor_x_1[i]-circle_x)<=x_ub) ) { // only need to consider those grids line intersect with the smallest bounding box that contains the circle

        double dx = grid_lines.coor_x_2[i] - grid_lines.coor_x_1[i];
        double dy = grid_lines.coor_y_2[i] - grid_lines.coor_y_1[i];
        double dr = sqrt(pow(dx,2.0)+pow(dy,2.0));
        double D = (grid_lines.coor_x_1[i]-circle_x)*(grid_lines.coor_y_2[i]-circle_y) - (grid_lines.coor_x_2[i]-circle_x)*(grid_lines.coor_y_1[i]-circle_y);

        int sgn_dy=1;
        if (dy<0) sgn_dy = -1;

        double Delta = pow(r,2.0)*pow(dr,2.0) - pow(D,2.0);

        switch(Delta>=0){
            case 1:{ // tangent (one intersection) or two intersection



                if (Delta>0){

                    double x_1 = (D*dy + sgn_dy*dx*sqrt(Delta))/pow(dr,2.0) + circle_x;
                    double x_2 = (D*dy - sgn_dy*dx*sqrt(Delta))/pow(dr,2.0) + circle_x;

                    double y_1 = (-D*dx + abs(dy)*sqrt(Delta))/pow(dr,2.0) + circle_y;
                    double y_2 = (-D*dx - abs(dy)*sqrt(Delta))/pow(dr,2.0) + circle_y;

                    double theta_1 = atan2((y_1-circle_y),(x_1-circle_x));
                    double theta_2 = atan2((y_2-circle_y),(x_2-circle_x));

                    if (theta_1<0) theta_1 = 2*M_PI + theta_1;
                    if (theta_2<0) theta_2 = 2*M_PI + theta_2;

                    // set_points.coor_x.push_back(x_1);
                    // set_points.coor_y.push_back(y_1);

                    // set_points.coor_x.push_back(x_2);
                    // set_points.coor_y.push_back(y_2);

                    // set_points.theta.push_back(theta_1);
                    // set_points.theta.push_back(theta_2);


                    set_points_struct tmp_1, tmp_2;

                    tmp_1.coor_x= x_1;
                    tmp_1.coor_y=y_1;

                    tmp_2.coor_x=x_2;
                    tmp_2.coor_y=y_2;

                    tmp_1.theta=theta_1;
                    tmp_2.theta=theta_2;

                    set_points.push_back(tmp_1);
                    set_points.push_back(tmp_2);


  
                }


                if (Delta==0){
                    double x = D*dy/pow(dr,2.0)+ circle_x;

                    double y = -D*dx/pow(dr,2.0) + circle_y;

                    double theta = atan2((y-circle_y),(x-circle_x));


                    if (theta<0) theta = 2*M_PI + theta;

                    // set_points.coor_x.push_back(x);
                    // set_points.coor_y.push_back(y);    

                    // set_points.theta.push_back(theta);

    
                    set_points_struct tmp;
                    
                    tmp.coor_x=(x);
                    tmp.coor_y=(y);    

                    tmp.theta=(theta);

                    set_points.push_back(tmp);

                }

            break;
            }

            case 0:{ // no intersection
                // do nothing
            break;
            }

        }

        // ofstream myfile_out; 
        // myfile_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);
        // myfile_out << Delta <<endl;
        // myfile_out.close();
    
    }



}

std::sort(set_points.begin(), set_points.end(), by_theta());

return(set_points);
}


// //----
// bool intersects(double circle_r, double circle_x, double circle_y, double grid_x, double grid_y, double grid_size)
// {
//     double circleDistance_x = abs(circle_x - grid_x);
//     double circleDistance_y = abs(circle_y - grid_y);

//     if (circleDistance_x > (grid_size/2 + circle_r)) { return false; }
//     if (circleDistance_y > (grid_size/2 + circle_r)) { return false; }

//     if (circleDistance_x <= (grid_size/2)) { return true; } 
//     if (circleDistance_y <= (grid_size/2)) { return true; }

//     double cornerDistance_sq = pow((circleDistance_x - grid_size/2),2) +
//                          pow((circleDistance_y - grid_size/2),2);

//     return (cornerDistance_sq < pow(circle_r,2));
// }

//---


vector<segments_struct> func_segments_attributes (vector<set_points_struct>& set_points, vector < vector<double> >& pop_grid, double& r, para_aux& para_other){ // return a set of segments (with length, density, and absolute angle) correspond to a set_points

    vector<segments_struct> segments;

    int n_segments = set_points.size();

    segments.resize(n_segments);

    // ofstream myfile_out_temp; 
    // myfile_out_temp.open((string(path4)+string("000.txt")).c_str(),ios::app);
    // myfile_out_temp << "n_segments" << endl;
    // myfile_out_temp << n_segments;
    // myfile_out_temp.close();

    for (int i=0; i<=(n_segments-1); i++){
        if(i==0) {

            segments[i].theta_abs = set_points[0].theta + (2*M_PI - set_points[set_points.size()-1].theta);

            double x_1= set_points[set_points.size()-1].coor_x;
            double y_1= set_points[set_points.size()-1].coor_y;
            double x_2= set_points[0].coor_x;
            double y_2= set_points[0].coor_y;
            double midpoint_x = (x_1+x_2)/2.0;
            double midpoint_y = (y_1+y_2)/2.0;
            int m = ceil((midpoint_y - para_other.y_min)/para_other.grid_size); // at mth row of the grid
            int n = ceil((midpoint_x - para_other.x_min)/para_other.grid_size); // at nth col..



            // switch(midpoint_x>para_other.x_max|midpoint_x<para_other.x_min|midpoint_y>para_other.y_max|midpoint_y<para_other.y_min){
            switch(m>para_other.n_row_grid | n>para_other.n_col_grid | m<=0| n<=0 ){
                case 1:{
                     segments[i].den = 0.0;
                break;
                }
                case 0:{

                    // ofstream myfile_out_temp; 
                    // myfile_out_temp.open((string(path4)+string("000.txt")).c_str(),ios::app);
                    // myfile_out_temp << "m" << "," << "n" << endl;
                    // myfile_out_temp << m << "," << n <<endl ;
                    // myfile_out_temp.close();

                    segments[i].den = pop_grid[m-1][n-1];
                break;
                }
            }


        }

        if(i!=0) {

            segments[i].theta_abs = set_points[i].theta - set_points[i-1].theta;

            double x_1= set_points[i-1].coor_x;
            double y_1= set_points[i-1].coor_y;
            double x_2= set_points[i].coor_x;
            double y_2= set_points[i].coor_y;
            double midpoint_x = (x_1+x_2)/2.0;
            double midpoint_y = (y_1+y_2)/2.0;
            int m = ceil((midpoint_y - para_other.y_min)/para_other.grid_size); // at mth row of the grid
            int n = ceil((midpoint_x - para_other.x_min)/para_other.grid_size); // at nth col..



            switch(m>para_other.n_row_grid | n>para_other.n_col_grid | m<=0| n<=0){
                case 1:{
                     segments[i].den = 0.0;
                break;
                }
                case 0:{

                    // ofstream myfile_out_temp; 
                    // myfile_out_temp.open((string(path4)+string("000.txt")).c_str(),ios::app);
                    // myfile_out_temp << "m" << "," << "n" << endl;
                    // myfile_out_temp << m << "," << n <<endl ;
                    // myfile_out_temp.close();

                    segments[i].den = pop_grid[m-1][n-1];
                break;
                }
            }

        }

        segments[i].len = r*(segments[i].theta_abs);  



    }



return(segments);

}


//-----
double dtnorm(double x, double mean, double sd, double a){ // pdf of truncated normal with lower bound = a; upper bound =Inf

double d;

switch(x>=a){
	case 1:{
		double num = (1.0/sd)*gsl_ran_ugaussian_pdf((x-mean)/sd);
		double denom = 1.0 - gsl_cdf_ugaussian_P((a-mean)/sd);
		d = num/denom;
	break;
	}
	
	case 0:{
		d = 0.0;
    	break;
	}
}

return(d);

}

/*-------------------------------------------*/
// double func_distance (double x_1 , double y_1, double x_2, double y_2){


// // x is Lat; y is Lon//

// double pi = 3.1415926535897;
// double rad = pi/180.0;
// double a1 = x_1 * rad;
// double a2 = y_1 * rad;
// double b1 = x_2 * rad;
// double b2 = y_2 * rad;
// double dlon = (b2 - a2);
// double dlat = (b1 - a1);
// double a = pow((sin(dlat/2.0)),2.0) + cos(a1) * cos(b1) * pow((sin(dlon/2.0)),2.0);
// double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

// double R = 6371.145;
// long double eucli_dist  = R * c;


// return(eucli_dist);

// }

/*-------------------------------------------*/

inline long double func_kernel (double x_1 , double y_1, double x_2, double y_2, double k_1_arg, double k_2_arg){


// x is Lat; y is Lon//

// double pi = 3.1415926535897;
// double rad = pi/180.0;
// double a1 = x_1 * rad;
// double a2 = y_1 * rad;
// double b1 = x_2 * rad;
// double b2 = y_2 * rad;
// double dlon = (b2 - a2);
// double dlat = (b1 - a1);
// double a = pow((sin(dlat/2.0)),2.0) + cos(a1) * cos(b1) * pow((sin(dlon/2.0)),2.0);
// double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

// double R = 6371.145;
// long double eucli_dist  = R * c;

long double eucli_dist = sqrt( (x_1-x_2)*(x_1-x_2) + (y_1-y_2)*(y_1-y_2) );

//double func_ker=0.0;

// if (kernel_type_arg=="exp_kernel") func_ker = exp((-par_kernel_1_arg)*pow(eucli_dist,par_kernel_2_arg)) ;
//if (kernel_type_arg=="exp_kernel") func_ker = 20*k_1_arg*exp((-k_1_arg)*pow(eucli_dist,k_2_arg)) ;
//if (kernel_type_arg=="exp_kernel") func_ker = 0.4*exp((-k_1_arg)*pow(eucli_dist,k_2_arg)) ;
//if (kernel_type_arg=="exp_kernel") func_ker = k_2_arg*exp((-k_1_arg)*eucli_dist) ;

// if (kernel_type_arg=="cauchy_kernel") func_ker = (1/(1+par_kernel_1_arg*pow((eucli_dist),par_kernel_2_arg)));
//if (kernel_type_arg=="cauchy_kernel") func_ker = 20*(1/(3.141592*k_2_arg*(1+pow((eucli_dist-k_1_arg)/k_2_arg,2))));

//long double func_ker = 20*k_1_arg*exp((-k_1_arg)*eucli_dist);
//long double func_ker = 20*(1/(3.141592*k_2_arg*(1+((eucli_dist)/k_2_arg)*((eucli_dist)/k_2_arg))));

//long double func_ker = 1.0/(1.0+pow(eucli_dist,k_2_arg)); //power-law

long double func_ker = exp((-k_1_arg)*eucli_dist); //exp
// long double func_ker = 1.0/(1+pow(eucli_dist,k_2_arg)); //power-law


return(func_ker);
}


/*-------------------------------------------*/

inline long double func_contact (double distance, double k_1_arg, double k_2_arg){


long double func_cont= k_1_arg*exp((-k_1_arg)*distance); //exp distribution

// long double func_cont = (2.0/(M_PI*k_1_arg))*(1.0/(1.0+pow((distance/k_1_arg),2.0)));

return(func_cont);
}

/*-------------------------------------------*/


inline long double func_distance (double x_1 , double y_1, double x_2, double y_2){


// x is Lat; y is Lon//

// double pi = 3.1415926535897;
// double rad = pi/180.0;
// double a1 = x_1 * rad;
// double a2 = y_1 * rad;
// double b1 = x_2 * rad;
// double b2 = y_2 * rad;
// double dlon = (b2 - a2);
// double dlat = (b1 - a1);
// double a = pow((sin(dlat/2.0)),2.0) + cos(a1) * cos(b1) * pow((sin(dlon/2.0)),2.0);
// double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

// double R = 6371.145;
// double eucli_dist  = R * c;

long double eucli_dist = sqrt( (x_1-x_2)*(x_1-x_2) + (y_1-y_2)*(y_1-y_2) );


return(eucli_dist);
}

/*----------------*/

long double func_distance_for_main (double x_1 , double y_1, double x_2, double y_2){


// x is Lat; y is Lon//

// double pi = 3.1415926535897;
// double rad = pi/180.0;
// double a1 = x_1 * rad;
// double a2 = y_1 * rad;
// double b1 = x_2 * rad;
// double b2 = y_2 * rad;
// double dlon = (b2 - a2);
// double dlat = (b1 - a1);
// double a = pow((sin(dlat/2.0)),2.0) + cos(a1) * cos(b1) * pow((sin(dlon/2.0)),2.0);
// double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

// double R = 6371.145;
// double eucli_dist  = R * c;


long double eucli_dist = sqrt( (x_1-x_2)*(x_1-x_2) + (y_1-y_2)*(y_1-y_2) );


return(eucli_dist);
}
/*-------------------------------------------*/
inline long double func_latent_pdf(double t , double mu_lat, double var_lat){

long double func_lat_pdf;

double a_lat = mu_lat*mu_lat/var_lat; // when use GAMMA latent
double b_lat = var_lat/mu_lat;
func_lat_pdf = gsl_ran_gamma_pdf(t, a_lat, b_lat);

// func_lat_pdf = gsl_ran_gamma_pdf(t, mu_lat, var_lat);


// double a_lat = log(mu_lat*mu_lat/(sqrt(var_lat +mu_lat*mu_lat))); // when use lognormal latent
// double b_lat = sqrt(log(var_lat/(mu_lat*mu_lat)+1));
// func_lat_pdf = gsl_ran_lognormal_pdf(t, a_lat, b_lat);

// double a_lat = mu_lat; // when use exponential latent
// func_lat_pdf = gsl_ran_exponential_pdf(t, a_lat);

return(func_lat_pdf);
}
/*-------------------------------------------*/

inline long double func_latent_cdf(double t , double mu_lat, double var_lat){

long double func_lat_cdf;

double a_lat = mu_lat*mu_lat/var_lat; // when use GAMMA latent
double b_lat = var_lat/mu_lat;
func_lat_cdf = gsl_cdf_gamma_P(t, a_lat, b_lat);

// func_lat_cdf = gsl_cdf_gamma_P(t, mu_lat, var_lat);

// double a_lat = log(mu_lat*mu_lat/(sqrt(var_lat +mu_lat*mu_lat))); // when use lognormal latent
// double b_lat = sqrt(log(var_lat/(mu_lat*mu_lat)+1));
// func_lat_cdf = gsl_cdf_lognormal_P(t, a_lat, b_lat);

// double a_lat = mu_lat; // when use exponential latent
// func_lat_cdf = gsl_cdf_exponential_P(t, a_lat);

return(func_lat_cdf);
}
/*-------------------------------------------*/

double log_lh_func (lh_SQUARE lh_square_arg, int n_arg) {

double log_lh_value =0.0;

for (int i=0; i<=(n_arg-1);i++){
// log_lh_value = log_lh_value + log(lh_square_arg.f_U.at(i)*lh_square_arg.f_E.at(i)*lh_square_arg.f_I.at(i)*lh_square_arg.f_R.at(i)*lh_square_arg.f_EnI.at(i)*lh_square_arg.f_InR.at(i));
log_lh_value = log_lh_value + log(lh_square_arg.f_U.at(i)*lh_square_arg.f_E.at(i)*lh_square_arg.f_I.at(i)*lh_square_arg.f_R.at(i)*lh_square_arg.f_EnI.at(i)*lh_square_arg.f_InR.at(i)*lh_square_arg.f_Grid.at(i)*lh_square_arg.f_Arc.at(i));

}

log_lh_value = log_lh_value + log(lh_square_arg.f_Surv);

return(log_lh_value);


}


/*------------------------------------------------*/
void FUNC::initialize_kernel_mat (vector< vector<double> >& kernel_mat_arg, vector<double>& norm_const_arg) {


for (int i=0;i<=(n_Clh-1);i++) {
 for (int j=0;j<=(n_Clh-1);j++) {
 if (i==j) kernel_mat_arg[i][j]=0.0;
 if (i<j) kernel_mat_arg[i][j] = func_kernel (coordinate_Clh[i][0],coordinate_Clh[i][1],coordinate_Clh[j][0],coordinate_Clh[j][1],k_1_Clh,k_2_Clh);
 if (i>j) kernel_mat_arg[i][j]=kernel_mat_arg[j][i];
 }
}

for (int j=0;j<=(n_Clh-1);j++) {
norm_const_arg.at(j) = 0.0;
 for (int i=0;(i<=(n_Clh-1)); i++) {
norm_const_arg.at(j) = norm_const_arg.at(j) +  kernel_mat_arg[i][j];
}
}

}

/*------------------------------------------------*/
void FUNC::initialize_delta_mat (vector< vector<double> >& delta_mat_arg){
 

if (xi_U_Clh.size()>=1){
for (int i=0;i<= (int)(xi_U_Clh.size()-1);i++){

    for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
    
    switch (t_r_Clh.at(xi_I_Clh.at(j))>t_max_Clh) {
    case 1:{ // not yet recovered
    delta_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)] = t_max_Clh - t_i_Clh.at(xi_I_Clh.at(j));
    break;
    }
    case 0:{ // recovered
    delta_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)] = t_r_Clh.at(xi_I_Clh.at(j)) - t_i_Clh.at(xi_I_Clh.at(j));
    break;
    }    
    }

    }
}
}


//----------//

if (xi_E_minus_Clh.size()>=1){
for (int i=0;i<= (int)(xi_E_minus_Clh.size()-1);i++){
   
    for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
   
        if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_minus_Clh.at(i))) { 

        switch (t_r_Clh.at(xi_I_Clh.at(j))>=t_e_Clh.at(xi_E_minus_Clh.at(i))) {
        case 1:{ // not yet recovered at e_i
        delta_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)] = t_e_Clh.at(xi_E_minus_Clh.at(i)) - t_i_Clh.at(xi_I_Clh.at(j));
        break;
        }
        case 0:{ // recovered before e_i
        delta_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)] = t_r_Clh.at(xi_I_Clh.at(j)) - t_i_Clh.at(xi_I_Clh.at(j));
        break;
        }    
        }

        } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)

}
}

for (int i=0;i<= (int)(index_Clh.size()-1);i++){
for (int j=0;j<= (int) (n_Clh-1);j++){
delta_mat_arg[index_Clh.at(i)][j] = 0.0;
}
}

}
/*------------------------------------------------*/

void FUNC::initialize_lh_square (lh_SQUARE& lh_square_arg, vector< vector<double> > kernel_mat_arg, vector< vector<double> > delta_mat_arg, vector<double>& norm_const_arg){



// if (xi_U_Clh.size()>=1){
// for (int i=0;i<= (int)(xi_U_Clh.size()-1);i++){

//     for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
 
//        double delta_t = delta_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)]; 
// //     double delta_t = 0.0;
// //     switch (t_r_Clh.at(xi_I_Clh.at(j))>t_max_Clh) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_Clh - t_i_Clh.at(xi_I_Clh.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_Clh.at(xi_I_Clh.at(j)) - t_i_Clh.at(xi_I_Clh.at(j));
// //     break;
// //     }    
// //     }

//     lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)) = lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)) + beta_age_Clh.at(age_gp_Clh.at(xi_I_Clh.at(j)))*delta_t*kernel_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)]/norm_const_arg.at(xi_I_Clh.at(j));
//     }

// //lh_square_arg.q_T.at(xi_U_Clh.at(i)) = alpha_Clh*t_max_Clh + beta_Clh*stb_Clh.at(xi_U_Clh.at(i))*lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)); // up to t_max


// lh_square_arg.q_T.at(xi_U_Clh.at(i)) = alpha_Clh*t_max_Clh*eta_Clh*stb_Clh.at(xi_U_Clh.at(i)) + eta_Clh*stb_Clh.at(xi_U_Clh.at(i))*lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)); // up to t_max
// // lh_square_arg.q_T.at(xi_U_Clh.at(i)) = alpha_Clh*t_r_Clh.at(xi_U_Clh.at(i))*eta_Clh*stb_Clh.at(xi_U_Clh.at(i)) + eta_Clh*stb_Clh.at(xi_U_Clh.at(i))*lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)); // up to the removal date


// lh_square_arg.f_U.at(xi_U_Clh.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_arg.q_T.at(xi_U_Clh.at(i)),1.0);

// }
// }

//----------//

// ofstream myfile_out; 
// myfile_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);
// for (int i=0;i<=(infected_source_Clh.size()-1);i++) {
// myfile_out << infected_source_Clh.at(i) <<endl;
// }
// myfile_out.close();



if (xi_E_minus_Clh.size()>=1){

for (int i=0;i<= (int)(xi_E_minus_Clh.size()-1);i++){
   
    // double total_beta =0.0;

    // for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){

    //     if (t_i_Clh.at(xi_I_Clh.at(j))<=t_e_Clh.at(xi_E_minus_Clh.at(i)) & t_r_Clh.at(xi_I_Clh.at(j))>t_e_Clh.at(xi_E_minus_Clh.at(i))) { 


    //         total_beta = total_beta + beta_age_Clh.at(age_gp_Clh.at(xi_I_Clh.at(j)));
    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)



    switch(infected_source_Clh.at(xi_E_minus_Clh.at(i))){
    case 9999:{ // by background (this part is redundant now but retained for future developments)
        // lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i)) = 0.0; // update k_sum_E

        double alpha_t= func_time_alpha( alpha_Clh, t_e_Clh.at(xi_E_minus_Clh.at(i)), para_other_Clh.T, omega_Clh);

        lh_square_arg.g_E.at(xi_E_minus_Clh.at(i)) = alpha_t;

        lh_square_arg.f_Grid.at(xi_E_minus_Clh.at(i)) = (pop_Clh[xi_E_minus_Clh.at(i)]/para_other_Clh.total_pop_of_grid)/pow(para_other_Clh.grid_size,2.0);
        lh_square_arg.f_Arc.at(xi_E_minus_Clh.at(i)) = 1.0;

        // lh_square_arg.f_Grid.at(xi_E_minus_Clh.at(i)) = 1.0/((para_other_Clh.x_max-para_other_Clh.x_min)*(para_other_Clh.y_max-para_other_Clh.y_min));
        // lh_square_arg.f_Arc.at(xi_E_minus_Clh.at(i)) = 1.0;


    break;
    }

    default :{ // not by background
        // lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i)) = kernel_mat_arg[xi_E_minus_Clh.at(i)][infected_source_Clh.at(xi_E_minus_Clh.at(i))]/norm_const_arg.at(infected_source_Clh.at(xi_E_minus_Clh.at(i))); // update k_sum_E
        double distance= func_distance (coordinate_Clh[infected_source_Clh.at(xi_E_minus_Clh.at(i))][0],coordinate_Clh[infected_source_Clh.at(xi_E_minus_Clh.at(i))][1],coordinate_Clh[xi_E_minus_Clh.at(i)][0],coordinate_Clh[xi_E_minus_Clh.at(i)][1]);
        
        double beta_t= func_time_beta( beta_age_Clh.at(age_gp_Clh.at(infected_source_Clh.at(xi_E_minus_Clh.at(i)))), t_e_Clh.at(xi_E_minus_Clh.at(i)), para_other_Clh.T, omega_Clh);
        lh_square_arg.g_E.at(xi_E_minus_Clh.at(i)) = beta_t*func_contact(distance, k_1_Clh, k_2_Clh);


        // ofstream myfile_out; 
        // myfile_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);
        // if (distance==0) myfile_out << xi_E_minus_Clh.at(i) << "," << infected_source_Clh.at(xi_E_minus_Clh.at(i)) << endl;
        // myfile_out.close();

        //---
     



        vector<set_points_struct> set_points= circle_line_intersections (coordinate_Clh[infected_source_Clh.at(xi_E_minus_Clh.at(i))][0],coordinate_Clh[infected_source_Clh.at(xi_E_minus_Clh.at(i))][1], distance, para_other_Clh.n_line, grid_lines_Clh);

            switch(set_points.size()>=1){
                case 1:{// have at least on intersection point
                    vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_Clh, distance, para_other_Clh);
                 
                    double mass = 0.0;
                  
                    for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
                        mass = mass + segments[i_segment].den*segments[i_segment].len;
                    }

                    lh_square_arg.f_Grid.at(xi_E_minus_Clh.at(i)) = pop_Clh[xi_E_minus_Clh.at(i)]/mass;
                    lh_square_arg.f_Arc.at(xi_E_minus_Clh.at(i)) = 1.0;

                break;
                }    

                case 0:{

                lh_square_arg.f_Grid.at(xi_E_minus_Clh.at(i)) = 1.0/(2.0*M_PI*distance);
                lh_square_arg.f_Arc.at(xi_E_minus_Clh.at(i)) = 1.0;


                break;
                }

            }



    //--

    break;
    }

    }


    // lh_square_arg.q_E.at(xi_E_minus_Clh.at(i)) = alpha_Clh*t_e_Clh.at(xi_E_minus_Clh.at(i))*eta_Clh*stb_Clh.at(xi_E_minus_Clh.at(i)) + eta_Clh*stb_Clh.at(xi_E_minus_Clh.at(i))*lh_square_arg.kt_sum_E.at(xi_E_minus_Clh.at(i));
    // lh_square_arg.q_E.at(xi_E_minus_Clh.at(i)) = alpha_Clh + lh_square_arg.kt_sum_E.at(xi_E_minus_Clh.at(i));

    // lh_square_arg.h_E.at(xi_E_minus_Clh.at(i)) = gsl_ran_exponential_pdf(lh_square_arg.q_E.at(xi_E_minus_Clh.at(i)),1.0);

    // lh_square_arg.f_E.at(xi_E_minus_Clh.at(i)) = lh_square_arg.g_E.at(xi_E_minus_Clh.at(i))*lh_square_arg.h_E.at(xi_E_minus_Clh.at(i));
    lh_square_arg.f_E.at(xi_E_minus_Clh.at(i)) = lh_square_arg.g_E.at(xi_E_minus_Clh.at(i));

}

}

//----------//


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_Clh;

t_change.insert( t_change.end(), t_r_Clh.begin(), t_r_Clh.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_Clh), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_Clh); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());


// ofstream myfile_out; 
// myfile_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);
// for (int i=0;i<= (int)(t_change.size()-1);i++){
//         if (i<(t_change.size()-1)) myfile_out  << t_change.at(i) << ",";
//         if (i==(t_change.size()-1)) myfile_out  << t_change.at(i) << endl;
// }
// myfile_out.close();

double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_Clh.at(xi_I_Clh.at(j))<=lb & t_r_Clh.at(xi_I_Clh.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            beta_j = beta_age_Clh.at(age_gp_Clh.at(xi_I_Clh.at(j)));

            if(para_other_Clh.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_Clh.T<lb){
                risk_surv_j = (beta_j/(-omega_Clh))*(exp(-omega_Clh*(ub-para_other_Clh.T ))-exp(-omega_Clh*(lb-para_other_Clh.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_Clh.T>=lb & para_other_Clh.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_Clh.T-lb);
                double risk_surv_j_2 = (beta_j/(-omega_Clh))*(exp(-omega_Clh*(ub-para_other_Clh.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }

            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// ofstream myfile_out; 
// myfile_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);
// myfile_out  << risk_surv << endl;

// myfile_out.close();


// risk_surv = risk_surv + alpha_Clh*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + alpha_Clh*( para_other_Clh.T - t_change[0] ) + (alpha_Clh/(-omega_Clh))*(exp(-omega_Clh*(t_change[t_change.size()-1]-para_other_Clh.T ))-1 ) ;


lh_square_arg.f_Surv = exp(-risk_surv);




//----------//

if (xi_I_Clh.size()>=1){
for (int i=0;i<= (int)(xi_I_Clh.size()-1);i++){
//lh_square_arg.f_I.at(xi_I_Clh.at(i)) = gsl_ran_gamma_pdf(t_i_Clh.at(xi_I_Clh.at(i)) - t_e_Clh.at(xi_I_Clh.at(i)), a_Clh, b_Clh);
lh_square_arg.f_I.at(xi_I_Clh.at(i)) = func_latent_pdf(t_i_Clh.at(xi_I_Clh.at(i)) - t_e_Clh.at(xi_I_Clh.at(i)), mu_lat_Clh,var_lat_Clh);

}
}

//--------//

if (xi_R_Clh.size()>=1){
for (int i=0;i<= (int)(xi_R_Clh.size()-1);i++){
// lh_square_arg.f_R.at(xi_R_Clh.at(i)) = gsl_ran_weibull_pdf(t_r_Clh.at(xi_R_Clh.at(i)) - t_i_Clh.at(xi_R_Clh.at(i)), c_Clh, d_Clh);
lh_square_arg.f_R.at(xi_R_Clh.at(i)) = gsl_ran_exponential_pdf(t_r_Clh.at(xi_R_Clh.at(i)) - t_i_Clh.at(xi_R_Clh.at(i)), c_Clh); // c_Clh is the mean

}
}

//-------//

if (xi_EnI_Clh.size()>=1){
for (int i=0;i<= (int)(xi_EnI_Clh.size()-1);i++){
//lh_square_arg.f_EnI.at(xi_EnI_Clh.at(i)) = 1.0 -  gsl_cdf_gamma_P(t_max_Clh - t_e_Clh.at(xi_EnI_Clh.at(i)), a_Clh, b_Clh);
lh_square_arg.f_EnI.at(xi_EnI_Clh.at(i)) = 1.0 -  func_latent_cdf(t_max_Clh - t_e_Clh.at(xi_EnI_Clh.at(i)),mu_lat_Clh,var_lat_Clh);
}
}
//-------//

if (xi_InR_Clh.size()>=1){
for (int i=0;i<= (int)(xi_InR_Clh.size()-1);i++){
// lh_square_arg.f_InR.at(xi_InR_Clh.at(i)) = 1.0 -  gsl_cdf_weibull_P(t_max_Clh - t_i_Clh.at(xi_InR_Clh.at(i)), c_Clh, d_Clh);
lh_square_arg.f_InR.at(xi_InR_Clh.at(i)) = 1.0 -  gsl_cdf_exponential_P(t_max_Clh - t_i_Clh.at(xi_InR_Clh.at(i)), c_Clh);

}
}

}

/*------------------------------------------------*/




void mcmc_UPDATE::alpha_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed

double alpha_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

alpha_proposed = para_current_arg.alpha + 0.1*gsl_ran_gaussian(r_c,1.0);

// switch (alpha_proposed<=0) {
// case 1: {
// alpha_proposed = -alpha_proposed; //reflection
// break;
// }
// case 0: {
// alpha_proposed = alpha_proposed;
// break;
// }
// }

switch (alpha_proposed<=0 | alpha_proposed>=10) {
case 1: {
alpha_proposed = para_current_arg.alpha; // keep the current value
}
case 0: {
alpha_proposed = alpha_proposed;
break;
}
}


// if (xi_U_arg.empty()==0){
// for (int i=0; i<=(int)(xi_U_arg.size()-1);i++){

// log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below

// // lh_square_modified.q_T.at(xi_U_arg.at(i)) = alpha_proposed*t_r_arg.at(xi_U_arg.at(i)) + stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
// lh_square_modified.q_T.at(xi_U_arg.at(i)) = alpha_proposed*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))+ para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));

// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }

//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   
 
    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 


    //         total_beta = total_beta + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));


    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background 

            double alpha_t = func_time_alpha (alpha_proposed, t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);

            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_t;

            // lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = (pop_CUPDATE[xi_E_minus_arg.at(i)]/para_other_CUPDATE.total_pop_of_grid)/pow(para_other_CUPDATE.grid_size,2.0);
            // lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;       
        break;
        }

        default :{ // not by background

            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);

            double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);

            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //     lh_square_modified.f_Grid.at(xi_E_minus_CUPDATE.at(i)) = 1.0/(2.0*M_PI*distance);
            //     lh_square_modified.f_Arc.at(xi_E_minus_CUPDATE.at(i)) = 1.0;


            //     break;
            //     }

            // }



    //--        
        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = alpha_proposed + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}

//----------//



// log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


// vector<double> t_change; // the time points where number of infectious changes

// t_change = t_i_arg;

// t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

// t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// // if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
// t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


// std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
// t_change.assign( s.begin(), s.end() );

// sort(t_change.begin(),t_change.end());

// double total_risk_secondary = (-log(lh_square_current_arg.f_Surv)) - (para_current_arg.alpha)*(para_other_CUPDATE.T -t_change[0]) - (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ; // the original total risk due to secondary cases
// double risk_surv_modified = total_risk_secondary +  (alpha_proposed)*(para_other_CUPDATE.T - t_change[0]) + (alpha_proposed/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;
// lh_square_modified.f_Surv = exp(-risk_surv_modified);


// log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//-----

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 

// double total_risk_secondary = (-log(lh_square_current_arg.f_Surv)) - (para_current_arg.alpha)*(t_max_CUPDATE- t_i_arg[index_arg[0]]); // the original total risk due to secondary cases
// double risk_surv_modified = total_risk_secondary + (alpha_proposed)*(t_max_CUPDATE- t_i_arg[index_arg[0]]);
// lh_square_modified.f_Surv = exp(-risk_surv_modified);

double total_risk_secondary = (-log(lh_square_current_arg.f_Surv)) - (para_current_arg.alpha)*(para_other_CUPDATE.T - t_i_arg[index_arg[0]]) - (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_max_CUPDATE-para_other_CUPDATE.T ))-1 ) ; // the original total risk due to secondary cases
double risk_surv_modified = total_risk_secondary +  (alpha_proposed)*(para_other_CUPDATE.T - t_i_arg[index_arg[0]]) + (alpha_proposed/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_max_CUPDATE-para_other_CUPDATE.T ))-1 ) ;
lh_square_modified.f_Surv = exp(-risk_surv_modified);


log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//-----
acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(para_current_arg.alpha/alpha_proposed));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("0.txt")).c_str(),ios::app);
// myfile_mcmc_out << t_change[t_change.size()-1] << ","<<  t_max_CUPDATE<< endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000_.txt")).c_str(),ios::app);
// myfile_mcmc_out << acp_pr<< "," << alpha_proposed << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.alpha = alpha_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.alpha = para_current_arg.alpha;
break;

}

gsl_rng_free(r_c);

}


//*------------------------------------------------*/

void mcmc_UPDATE::omega_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

double omega_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;


const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed

// omega_proposed = para_current_arg.omega + 0.1*gsl_ran_gaussian(r_c,1.0);
// // gsl_rng_free(r_c);
// switch (omega_proposed<=0) {
// case 1: {
// omega_proposed = -omega_proposed; //reflection
// break;
// }
// case 0: {
// omega_proposed = omega_proposed;
// break;
// }
// }

while(omega_proposed<=0){
omega_proposed = para_current_arg.omega + 0.1*gsl_ran_gaussian(r_c,1.0);
}


//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

      
//     }

// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }

//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background   

             double alpha_t = func_time_alpha (para_current_arg.alpha, t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, omega_proposed);

            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_t;

        break;
        }

        default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);


            double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, omega_proposed);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);


            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = 1.0/(2.0*M_PI*distance);
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;


            //     break;
            //     }

            // }

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_arg;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());


double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_arg.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            double beta_j = para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-omega_proposed))*(exp(-omega_proposed*(ub-para_other_CUPDATE.T ))-exp(-omega_proposed*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-omega_proposed))*(exp(-omega_proposed*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }
                                   

            // risk_surv_j = (beta_j/(-omega_proposed))*(exp(-omega_proposed*(ub-para_other_CUPDATE.T ))-exp(-omega_proposed*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval


            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );

risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-omega_proposed))*(exp(-omega_proposed*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;


lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//---

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.omega = omega_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.omega = para_current_arg.omega;
break;

}

gsl_rng_free(r_c);

}



//*------------------------------------------------*/

void mcmc_UPDATE::beta_1_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

double beta_1_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector<double> beta_age_modified = para_current_arg.beta_age;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
beta_1_proposed = para_current_arg.beta_age.at(0) + 0.3*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (beta_1_proposed<=0) {
case 1: {
beta_1_proposed = -beta_1_proposed; //reflection
break;
}
case 0: {
beta_1_proposed = beta_1_proposed;
break;
}
}

beta_age_modified.at(0) = beta_1_proposed;


//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

      
//     }

// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }

//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background   
        break;
        }

        default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);


            double beta_t = func_time_beta (beta_age_modified.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);


            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = 1.0/(2.0*M_PI*distance);
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;


            //     break;
            //     }

            // }

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_arg;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());


double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_arg.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            beta_j = beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }
                 

            // risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval


            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;


lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//---

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.beta_age = beta_age_modified;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.beta_age = para_current_arg.beta_age;
break;

}

gsl_rng_free(r_c);

}


/*------------------------------------------------*/


void mcmc_UPDATE::beta_2_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

double beta_2_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector<double> beta_age_modified = para_current_arg.beta_age;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
beta_2_proposed = para_current_arg.beta_age.at(1) + 0.3*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (beta_2_proposed<=0) {
case 1: {
beta_2_proposed = -beta_2_proposed; //reflection
break;
}
case 0: {
beta_2_proposed = beta_2_proposed;
break;
}
}

beta_age_modified.at(1) = beta_2_proposed;


//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

      
//     }

// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }
//----------
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background   
        break;
        }

        default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);


            double beta_t = func_time_beta (beta_age_modified.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);


            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = 1.0/(2.0*M_PI*distance);
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;


            //     break;
            //     }

            // }

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_arg;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());

double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_arg.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            beta_j = beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }
                 

            // risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval


            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;

lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//---


acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.beta_age = beta_age_modified;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.beta_age = para_current_arg.beta_age;
break;

}

gsl_rng_free(r_c);

}


/*------------------------------------------------*/

void mcmc_UPDATE::beta_3_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

double beta_3_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector<double> beta_age_modified = para_current_arg.beta_age;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
beta_3_proposed = para_current_arg.beta_age.at(2) + 0.3*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (beta_3_proposed<=0) {
case 1: {
beta_3_proposed = -beta_3_proposed; //reflection
break;
}
case 0: {
beta_3_proposed = beta_3_proposed;
break;
}
}

beta_age_modified.at(2) = beta_3_proposed;


//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

      
//     }

// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }
//----------
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background   
        break;
        }

        default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);


            double beta_t = func_time_beta (beta_age_modified.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);


            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = 1.0/(2.0*M_PI*distance);
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;


            //     break;
            //     }

            // }

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_arg;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());

double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_arg.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            beta_j = beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }
                 

            // risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval


            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;

lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//---

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.beta_age = beta_age_modified;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.beta_age = para_current_arg.beta_age;
break;

}

gsl_rng_free(r_c);

}


/*------------------------------------------------*/

void mcmc_UPDATE::beta_4_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

double beta_4_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector<double> beta_age_modified = para_current_arg.beta_age;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
beta_4_proposed = para_current_arg.beta_age.at(3) + 0.3*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (beta_4_proposed<=0) {
case 1: {
beta_4_proposed = -beta_4_proposed; //reflection
break;
}
case 0: {
beta_4_proposed = beta_4_proposed;
break;
}
}

beta_age_modified.at(3) = beta_4_proposed;


//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

      
//     }

// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }
//----------
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background   
        break;
        }

        default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);


            double beta_t = func_time_beta (beta_age_modified.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);


            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = 1.0/(2.0*M_PI*distance);
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;


            //     break;
            //     }

            // }

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_arg;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());

double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_arg.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            beta_j = beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }
                 

            // risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval


            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;

lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//---

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.beta_age = beta_age_modified;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.beta_age = para_current_arg.beta_age;
break;

}

gsl_rng_free(r_c);

}



/*------------------------------------------------*/


void mcmc_UPDATE::beta_5_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

double beta_5_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector<double> beta_age_modified = para_current_arg.beta_age;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
beta_5_proposed = para_current_arg.beta_age.at(4) + 0.3*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (beta_5_proposed<=0) {
case 1: {
beta_5_proposed = -beta_5_proposed; //reflection
break;
}
case 0: {
beta_5_proposed = beta_5_proposed;
break;
}
}

beta_age_modified.at(4) = beta_5_proposed;


//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

      
//     }

// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }
//----------
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    // double total_beta = 0.0;

    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background   
        break;
        }

        default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);


            double beta_t = func_time_beta (beta_age_modified.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);


            // vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1], distance, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

            // switch(set_points.size()>=1){
            //     case 1:{// have at least on intersection point
            //         vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance, para_other_CUPDATE);
                 
            //         double mass = 0.0;
                  
            //         for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
            //             mass = mass + segments[i_segment].den*segments[i_segment].len;
            //         }

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = pop_CUPDATE[xi_E_minus_arg.at(i)]/mass;
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;

            //     break;
            //     }    

            //     case 0:{

            //         lh_square_modified.f_Grid.at(xi_E_minus_arg.at(i)) = 1.0/(2.0*M_PI*distance);
            //         lh_square_modified.f_Arc.at(xi_E_minus_arg.at(i)) = 1.0;


            //     break;
            //     }

            // }

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_arg;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());

double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_arg.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
            // n_dt = n_dt + 1.0;

            beta_j = beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }
                 

            // risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval


            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;

lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//---

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.beta_age = beta_age_modified;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.beta_age = para_current_arg.beta_age;
break;

}

gsl_rng_free(r_c);

}



/*------------------------------------------------*/

void mcmc_UPDATE::mu_lat_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg,const vector<double>& t_r_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector<int>& infected_source_arg,int iter){

double mu_lat_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed

mu_lat_proposed = para_current_arg.mu_lat + 0.5*gsl_ran_gaussian(r_c,1.0);
//mu_lat_proposed = para_current_arg.mu_lat + gsl_ran_flat(r_c,-0.05,0.05);

// double log_mu_lat_proposed;
// log_mu_lat_proposed = log(para_current_arg.mu_lat) + 0.1*gsl_ran_gaussian(r_c,1.0);
// mu_lat_proposed = exp(log_mu_lat_proposed );

// gsl_rng_free(r_c);

switch (mu_lat_proposed<=0 | mu_lat_proposed>=20) {
case 1: {
mu_lat_proposed = para_current_arg.mu_lat; 
break;
}
case 0: {
mu_lat_proposed = mu_lat_proposed;
break;
}
}

if (xi_I_arg.empty()==0){
for (int i=0; i<=(int)(xi_I_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_I.at(xi_I_arg.at(i)) = func_latent_pdf(t_i_arg.at(xi_I_arg.at(i)) - t_e_arg.at(xi_I_arg.at(i)), mu_lat_proposed, para_current_arg.var_lat);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_EnI_arg.empty()==0){
for (int i=0; i<=(int)(xi_EnI_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_EnI.at(xi_EnI_arg.at(i)) =  1.0 - func_latent_cdf(t_max_CUPDATE - t_e_arg.at(xi_EnI_arg.at(i)), mu_lat_proposed, para_current_arg.var_lat);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //add back part of likelihood that updated above
}
}


acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(para_current_arg.mu_lat/mu_lat_proposed)); // when sampling on log of a


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_mu_lat.txt")).c_str(),ios::app);
// myfile_mcmc_out << mu_lat_proposed << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.mu_lat = mu_lat_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.mu_lat = para_current_arg.mu_lat;
break;

}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::var_lat_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<double>& t_r_arg,const vector<int>& index_arg, para_key& para_current_arg,const vector<int>& infected_source_arg, int iter){

double var_lat_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed

var_lat_proposed = para_current_arg.var_lat + 1*gsl_ran_gaussian(r_c,1.0);
//var_lat_proposed = para_current_arg.var_lat + gsl_ran_flat(r_c,-0.03,0.03);

// double log_var_lat_proposed;
// log_var_lat_proposed = log(para_current_arg.var_lat) + 0.1*gsl_ran_gaussian(r_c,1.0);
// var_lat_proposed = exp(log_var_lat_proposed);

switch (var_lat_proposed<=0) {
case 1: {
var_lat_proposed = -var_lat_proposed; //reflection
break;
}
case 0: {
var_lat_proposed = var_lat_proposed;
break;
}
}

if (xi_I_arg.empty()==0){
for (int i=0; i<=(int)(xi_I_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_I.at(xi_I_arg.at(i)) = func_latent_pdf(t_i_arg.at(xi_I_arg.at(i)) - t_e_arg.at(xi_I_arg.at(i)), para_current_arg.mu_lat, var_lat_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_EnI_arg.empty()==0){
for (int i=0; i<=(int)(xi_EnI_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_EnI.at(xi_EnI_arg.at(i)) =  1.0 - func_latent_cdf(t_max_CUPDATE - t_e_arg.at(xi_EnI_arg.at(i)), para_current_arg.mu_lat, var_lat_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //add back part of likelihood that updated above
}
}

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(para_current_arg.var_lat/var_lat_proposed));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_var_lat.txt")).c_str(),ios::app);
// myfile_mcmc_out << var_lat_proposed<< "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.var_lat = var_lat_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.var_lat = para_current_arg.var_lat;
break;

}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::c_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg,const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector<int>& infected_source_arg,int iter){

double c_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
c_proposed = para_current_arg.c + 0.03*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (c_proposed<=0) {
case 1: {
c_proposed = -c_proposed; //reflection
break;
}
case 0: {
c_proposed = c_proposed;
break;
}
}

if (xi_R_arg.empty()==0){
for (int i=0; i<=(int)(xi_R_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //subtract part of likelihood that would be updated below

// lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_weibull_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), c_proposed, para_current_arg.d);
lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_exponential_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), c_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_InR_arg.empty()==0){
for (int i=0; i<=(int)(xi_InR_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //subtract part of likelihood that would be updated below

// lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_weibull_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)), c_proposed, para_current_arg.d);
lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)), c_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //add back part of likelihood that updated above
}
}

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.c = c_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.c = para_current_arg.c;
break;

}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::d_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg,  const vector<int>& index_arg, para_key& para_current_arg,const vector<int>& infected_source_arg, int iter){

double d_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
d_proposed = para_current_arg.d + 0.01*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (d_proposed<=0) {
case 1: {
d_proposed = -d_proposed; //reflection
break;
}
case 0: {
d_proposed = d_proposed;
break;
}
}

if (xi_R_arg.empty()==0){
for (int i=0; i<=(int)(xi_R_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //subtract part of likelihood that would be updated below

// lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_weibull_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), para_current_arg.c, d_proposed );
lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_exponential_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), para_current_arg.c );

log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_InR_arg.empty()==0){
for (int i=0; i<=(int)(xi_InR_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //subtract part of likelihood that would be updated below

// lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_weibull_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)),  para_current_arg.c, d_proposed);
lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)),  para_current_arg.c);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //add back part of likelihood that updated above
}
}

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.d = d_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.d = para_current_arg.d;
break;
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::k_1_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg,const vector<int>& infected_source_arg, int iter){

double k_1_proposed = 0.0;
double acp_pr = 0.0;

double log_prior_backward = 0.0;
double log_prior_forward =0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
// vector< vector<double> > kernel_mat_modified = kernel_mat_current_arg;
// vector<double> norm_const_modified = norm_const_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
//k_1_proposed = para_current_arg.k_1 + 1.0*gsl_ran_flat(r_c,0.0,0.01);
k_1_proposed = para_current_arg.k_1 + 1*gsl_ran_gaussian(r_c,1.0);


// gsl_rng_free(r_c);

switch (k_1_proposed<=0.0) {
case 1: {
k_1_proposed = -k_1_proposed;//reflection
break;
}
case 0: {
k_1_proposed = k_1_proposed;
break;
}
}

// switch (k_1_proposed<=0.03) {
// case 1: {
// k_1_proposed = 0.03+abs(0.03-k_1_proposed);//reflection
// break;
// }
// case 0: {
// k_1_proposed = k_1_proposed;
// break;
// }
// }

// log_prior_backward = log(dtnorm(k_1_proposed, 0.2, 0.02, 0.0));
// log_prior_forward = log(dtnorm(para_current_arg.k_1, 0.2, 0.02, 0.0));

//----------
// for (int i=0;i<=(n_CUPDATE-1);i++) {
//  for (int j=0;j<=(n_CUPDATE-1);j++) {
//  if (i==j) kernel_mat_modified[i][j]=0.0;
//  if (i<j) kernel_mat_modified[i][j] = func_kernel (coordinate_CUPDATE[i][0],coordinate_CUPDATE[i][1],coordinate_CUPDATE[j][0],coordinate_CUPDATE[j][1],k_1_proposed,para_current_arg.k_2);
//  if (i>j) kernel_mat_modified[i][j]=kernel_mat_modified[j][i];
//  }
// }

// for (int j=0;j<=(n_CUPDATE-1);j++) {
// norm_const_modified.at(j)=0.0;
//  for (int i=0;(i<=(n_CUPDATE-1));i++) {
// norm_const_modified.at(j)= norm_const_modified.at(j) + kernel_mat_modified[i][j];
// }
// }

//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_modified[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j));

      
//     }

// // lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }

//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
 

    // double total_beta = 0.0;

   
    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   
    

    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background (this part is redundant now but retained for future developments)


        double alpha_t = func_time_alpha (para_current_arg.alpha, t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);

        lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_t;

        // lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha;
        
        break;
        }

        default :{ // not by background
        double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);
  

        double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega );
        lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = beta_t*stb_arg.at(xi_E_minus_arg.at(i))*func_contact(distance, k_1_proposed,para_current_arg.k_2);
  

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}
//----------

acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) +(log_prior_forward - log_prior_backward)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("k_1_proposed_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  k_1_proposed << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();



switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
// kernel_mat_current_arg = kernel_mat_modified;
// norm_const_current_arg = norm_const_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.k_1 = k_1_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.k_1 = para_current_arg.k_1;
break;
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


/*------------------------------------------------*/

void mcmc_UPDATE::k_2_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , vector<double>& norm_const_current_arg,const vector<int>& infected_source_arg, int iter){

double k_2_proposed = 0.0;
double acp_pr = 0.0;

double log_prior_backward = 0.0;
double log_prior_forward =0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
// vector< vector<double> > kernel_mat_modified = kernel_mat_current_arg;
// vector<double> norm_const_modified = norm_const_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
//k_1_proposed = para_current_arg.k_1 + 1.0*gsl_ran_flat(r_c,0.0,0.01);
k_2_proposed = para_current_arg.k_2 + 0.001*gsl_ran_gaussian(r_c,1.0);


// gsl_rng_free(r_c);

switch (k_2_proposed<=0.0) {
case 1: {
k_2_proposed = -k_2_proposed;//reflection
break;
}
case 0: {
k_2_proposed = k_2_proposed;
break;
}
}

// switch (k_1_proposed<=0.03) {
// case 1: {
// k_1_proposed = 0.03+abs(0.03-k_1_proposed);//reflection
// break;
// }
// case 0: {
// k_1_proposed = k_1_proposed;
// break;
// }
// }

// log_prior_backward = log(dtnorm(k_1_proposed, 0.2, 0.02, 0.0));
// log_prior_forward = log(dtnorm(para_current_arg.k_1, 0.2, 0.02, 0.0));

//----------
// for (int i=0;i<=(n_CUPDATE-1);i++) {
//  for (int j=0;j<=(n_CUPDATE-1);j++) {
//  if (i==j) kernel_mat_modified[i][j]=0.0;
//  if (i<j) kernel_mat_modified[i][j] = func_kernel (coordinate_CUPDATE[i][0],coordinate_CUPDATE[i][1],coordinate_CUPDATE[j][0],coordinate_CUPDATE[j][1],k_1_proposed,para_current_arg.k_2);
//  if (i>j) kernel_mat_modified[i][j]=kernel_mat_modified[j][i];
//  }
// }

// for (int j=0;j<=(n_CUPDATE-1);j++) {
// norm_const_modified.at(j)=0.0;
//  for (int i=0;(i<=(n_CUPDATE-1));i++) {
// norm_const_modified.at(j)= norm_const_modified.at(j) + kernel_mat_modified[i][j];
// }
// }

//----------

// if (xi_U_arg.empty()==0){
// for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

//     log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
//     for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

// //     double delta_t = 0.0;
// //     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
// //     case 1:{ // not yet recovered
// //     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }
// //     case 0:{ // recovered
// //     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// //     break;
// //     }    
// //     }

//     double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

//     lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t*kernel_mat_modified[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j));

      
//     }

// // lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
// lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));


// lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
// }
// }
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){

    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;


    // double total_beta = 0.0;
    
    // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   


    //     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

    //         total_beta = total_beta + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));



    //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


    switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
        case 9999:{ // by background (this part is redundant now but retained for future developments)
        double alpha_t = func_time_alpha (para_current_arg.alpha, t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);

        lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_t;
        break;
        }

        default :{ // not by background
        double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);
   
        double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega );
        lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = beta_t*stb_arg.at(xi_E_minus_arg.at(i))*func_contact(distance, para_current_arg.k_1,k_2_proposed);
  

        break;
        }

    }

    // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
    // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

    lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));

    log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}
//----------
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) +(log_prior_forward - log_prior_backward)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("k_1_proposed_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  k_1_proposed << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();



switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
// kernel_mat_current_arg = kernel_mat_modified;
// norm_const_current_arg = norm_const_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.k_2 = k_2_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.k_2 = para_current_arg.k_2;
break;
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::t_i_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int subject_proposed){

double t_proposed; // new t_i to be proposed
double t_low, t_up;

double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
// vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_i_modified = t_i_arg;


const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,subject_proposed); // set a seed

// int subject_proposed = xi_I_arg.at(gsl_rng_uniform_int (r_c, xi_I_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn

//--
vector<int> infecting_list;
vector<double> t_e_list;
infecting_list.reserve(xi_E_minus_arg.size());
t_e_list.reserve(xi_E_minus_arg.size());

for (int j=0;j<=(int) (xi_E_minus_arg.size()-1);j++){ 
    if (infected_source_arg.at(xi_E_minus_arg.at(j))==subject_proposed) {
        infecting_list.push_back(xi_E_minus_arg.at(j));
        t_e_list.push_back(t_e_arg.at(xi_E_minus_arg.at(j)));
    }
}


//--

switch(infecting_list.size()>=1){
    case 1:{
        double min_t_e = *min_element(t_e_list.begin(),t_e_list.end());
        t_up =  min_t_e;

        double dt = t_r_arg.at(subject_proposed) - min_t_e;
        t_low = std::max(t_e_arg.at(subject_proposed), t_up-(20-dt));  //make sure the biggest difference between t_i and t_r is  20


// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  dt << endl;
// myfile_mcmc_out.close();

      
    break;
    }
    case 0:{
        t_up =  min(t_r_arg.at(subject_proposed), t_max_CUPDATE);
        t_low = std::max(t_e_arg.at(subject_proposed), t_up-20);

    break;
    }
}


// t_low = std::max(t_e_arg.at(subject_proposed), t_up-20);


// t_proposed= gsl_ran_flat(r_c,t_low, t_up);

t_proposed = 9e+10;
while(t_proposed<t_low | t_proposed>t_up){
t_proposed= t_i_arg.at(subject_proposed) + 0.1*gsl_ran_gaussian(r_c,1.0);
}


t_i_modified.at(subject_proposed) = t_proposed;


double log_prior_y =0.0;
double log_prior_x=0.0;


// log_prior_y = log(gsl_ran_gamma_pdf(t_proposed, t_o/pow(0.5,2.0), pow(0.5,2.0)));
// log_prior_x = log(gsl_ran_gamma_pdf(t_i_arg.at(subject_proposed), t_o/pow(0.5,2.0), pow(0.5,2.0)));


// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  min_t_e << endl;
// myfile_mcmc_out.close();




//----------------------------------------------------------------------------------//

// for (int j=0;j<=(int) (xi_U_arg.size()-1);j++){ 

// log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(j))); //subtract part of likelihood that would be updated below


// switch (t_r_arg.at(subject_proposed)>=t_max_CUPDATE) {
//     case 1:{
//         delta_mat_modified[xi_U_arg.at(j)][subject_proposed] = t_max_CUPDATE  - t_proposed;
//     break;
//     }
//     case 0:{
//         delta_mat_modified[xi_U_arg.at(j)][subject_proposed] = t_r_arg.at(subject_proposed) - t_proposed;
//     break;
//     }
// }   


// lh_square_modified.kt_sum_U.at(xi_U_arg.at(j)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(j)) - para_current_arg.beta_age.at(age_gp_CUPDATE.at(subject_proposed))*delta_mat_current_arg[xi_U_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_U_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(subject_proposed))*delta_mat_modified[xi_U_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_U_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed)  ; //subtract the infectious challenge  due to the infectious subject chosen THEN add the updated one back

// // lh_square_modified.q_T.at(xi_U_arg.at(j)) = para_current_arg.alpha*stb_arg.at(xi_U_arg.at(j))* t_max_CUPDATE  + para_current_arg.beta*stb_arg.at(xi_U_arg.at(j))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(j));
// lh_square_modified.q_T.at(xi_U_arg.at(j)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(j))  + para_current_arg.eta*stb_arg.at(xi_U_arg.at(j))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(j));

// lh_square_modified.f_U.at(xi_U_arg.at(j)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(j)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(j))); // add back the updated part of likelihood

// }
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    // switch ((t_e_arg.at(xi_E_minus_arg.at(i))>=t_proposed) ) {

    // case 1:{



        log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

        // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
        // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;


        // double total_beta = 0.0;

        // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
       
 
        //     if (t_i_modified.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

        //         total_beta = total_beta + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

        //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


        // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


        switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
            case 9999:{ // by background
            double alpha_t = func_time_alpha (para_current_arg.alpha, t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);

            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_t;
            break;
            }

            default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);

            double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);

            break;
            }

        }

        // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
        // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

        // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

        // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
        lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

        log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // break;
    // }

    // case 0:{
    // break;
    // }
    // }
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_modified;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());


double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_modified.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
  

            beta_j = para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }

            // risk_surv_j = (para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval

            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;

lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//----------

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_proposed - t_e_arg.at(subject_proposed), para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 


//----------

switch ( find(xi_R_arg.begin(), xi_R_arg.end(),subject_proposed) != (xi_R_arg.end()) ) { //return 1 when the subject is also in xi_R
case 1:{
    log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(subject_proposed)); //subtract part of likelihood that would be updated below
    // lh_square_modified.f_R.at(subject_proposed) = gsl_ran_weibull_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c, para_current_arg.d);
    lh_square_modified.f_R.at(subject_proposed) = gsl_ran_exponential_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c);
    log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(subject_proposed)); 
    break;
    }
    case 0:{
    break;
    }
}

//----------
switch ( find(xi_InR_arg.begin(), xi_InR_arg.end(),subject_proposed) != (xi_InR_arg.end()) ) { //return 1 when the subject is also in xi_InR
    case 1:{
    log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(subject_proposed)); //subtract part of likelihood that would be updated below
    // lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_weibull_P(t_max_CUPDATE - t_proposed, para_current_arg.c, para_current_arg.d);
    lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_exponential_P(t_max_CUPDATE - t_proposed, para_current_arg.c);
    log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(subject_proposed)); 
    
    break;
    }
    case 0:{
    break;
    }
}

//---------------
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) +(log_prior_y-log_prior_x)));

if ( isnan(log_lh_modified)==1 |  isinf(log_lh_modified)==1 ) acp_pr=0;


//      ofstream myfile_mcmc_out; 
//      myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
//      myfile_mcmc_out << acp_pr  << "," << t_i_arg.at(subject_proposed) << "," << t_proposed << endl;
//      myfile_mcmc_out.close();

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
    case 1: {
        lh_square_current_arg = lh_square_modified;
        // delta_mat_current_arg = delta_mat_modified;
        log_lh_current_arg = log_lh_modified;
        t_i_arg= t_i_modified;
    break;
    }
    
    case 0: {
    break;
    }
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


/*------------------------------------------------*/


void mcmc_UPDATE::t_i_update_2(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int subject_proposed){

double t_proposed; // new t_i to be proposed
double t_low, t_up;

double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
// vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_i_modified = t_i_arg;


const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,subject_proposed); // set a seed

// int subject_proposed = xi_I_arg.at(gsl_rng_uniform_int (r_c, xi_I_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn

//--
vector<int> infecting_list;
vector<double> t_e_list;
infecting_list.reserve(xi_E_minus_arg.size());
t_e_list.reserve(xi_E_minus_arg.size());

for (int j=0;j<=(int) (xi_E_minus_arg.size()-1);j++){ 
    if (infected_source_arg.at(xi_E_minus_arg.at(j))==subject_proposed) {
        infecting_list.push_back(xi_E_minus_arg.at(j));
        t_e_list.push_back(t_e_arg.at(xi_E_minus_arg.at(j)));
    }
}


//--

switch(infecting_list.size()>=1){
    case 1:{
        double min_t_e = *min_element(t_e_list.begin(),t_e_list.end());
        t_up =  min_t_e;

        double dt = t_r_arg.at(subject_proposed) - min_t_e;
        t_low = std::max(t_e_arg.at(subject_proposed), t_up-(1.0-dt));  //make sure the biggest difference between t_i and t_r is  20


// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  dt << endl;
// myfile_mcmc_out.close();

      
    break;
    }
    case 0:{
        t_up =  min(t_r_arg.at(subject_proposed), t_max_CUPDATE);
        t_low = std::max(t_e_arg.at(subject_proposed), t_up-1.0);

    break;
    }
}


// t_low = std::max(t_e_arg.at(subject_proposed), t_up-20);


// t_proposed= gsl_ran_flat(r_c,t_low, t_up);

t_proposed = 9e+10;
while(t_proposed<t_low | t_proposed>t_up){
t_proposed= t_i_arg.at(subject_proposed) + 0.1*gsl_ran_gaussian(r_c,1.0);
}


t_i_modified.at(subject_proposed) = t_proposed;


double log_prior_y =0.0;
double log_prior_x=0.0;


// log_prior_y = log(gsl_ran_gamma_pdf(t_proposed, t_o/pow(0.5,2.0), pow(0.5,2.0)));
// log_prior_x = log(gsl_ran_gamma_pdf(t_i_arg.at(subject_proposed), t_o/pow(0.5,2.0), pow(0.5,2.0)));


// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  min_t_e << endl;
// myfile_mcmc_out.close();




//----------------------------------------------------------------------------------//

// for (int j=0;j<=(int) (xi_U_arg.size()-1);j++){ 

// log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(j))); //subtract part of likelihood that would be updated below


// switch (t_r_arg.at(subject_proposed)>=t_max_CUPDATE) {
//     case 1:{
//         delta_mat_modified[xi_U_arg.at(j)][subject_proposed] = t_max_CUPDATE  - t_proposed;
//     break;
//     }
//     case 0:{
//         delta_mat_modified[xi_U_arg.at(j)][subject_proposed] = t_r_arg.at(subject_proposed) - t_proposed;
//     break;
//     }
// }   


// lh_square_modified.kt_sum_U.at(xi_U_arg.at(j)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(j)) - para_current_arg.beta_age.at(age_gp_CUPDATE.at(subject_proposed))*delta_mat_current_arg[xi_U_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_U_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(subject_proposed))*delta_mat_modified[xi_U_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_U_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed)  ; //subtract the infectious challenge  due to the infectious subject chosen THEN add the updated one back

// // lh_square_modified.q_T.at(xi_U_arg.at(j)) = para_current_arg.alpha*stb_arg.at(xi_U_arg.at(j))* t_max_CUPDATE  + para_current_arg.beta*stb_arg.at(xi_U_arg.at(j))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(j));
// lh_square_modified.q_T.at(xi_U_arg.at(j)) = para_current_arg.alpha*t_max_CUPDATE*para_current_arg.eta*stb_arg.at(xi_U_arg.at(j))  + para_current_arg.eta*stb_arg.at(xi_U_arg.at(j))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(j));

// lh_square_modified.f_U.at(xi_U_arg.at(j)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(j)),1.0);

// log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(j))); // add back the updated part of likelihood

// }
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    // switch ((t_e_arg.at(xi_E_minus_arg.at(i))>=t_proposed) ) {

    // case 1:{



        log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

        // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
        // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;


        // double total_beta = 0.0;

        // for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
       
 
        //     if (t_i_modified.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

        //         total_beta = total_beta + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

        //     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


        // } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


        switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
            case 9999:{ // by background
            double alpha_t = func_time_alpha (para_current_arg.alpha, t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);

            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_t;
            break;
            }

            default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);

            double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i)))), t_e_arg.at(xi_E_minus_arg.at(i)), para_other_CUPDATE.T, para_current_arg.omega);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = beta_t*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2);

            break;
            }

        }

        // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
        // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

        // lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

        // lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
        lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));

        log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // break;
    // }

    // case 0:{
    // break;
    // }
    // }
}
}


//----------//

log_lh_modified = log_lh_modified - log(lh_square_current_arg.f_Surv); 


vector<double> t_change; // the time points where number of infectious changes

t_change = t_i_modified;

t_change.insert( t_change.end(), t_r_arg.begin(), t_r_arg.end() ); // join t_i and t_r

t_change.erase(std::remove(t_change.begin(), t_change.end(), unassigned_time_CUPDATE), t_change.end()); // remove unassigned_time


// if (t_change[t_change.size()-1]!=t_max_Clh) t_change.push_back(t_max_Clh); // the last time point should be t_max
t_change.push_back(t_max_CUPDATE); // the last time point should be t_max


std::set<double> s( t_change.begin(), t_change.end() ); // return a unique set without the duplicated elements
t_change.assign( s.begin(), s.end() );

sort(t_change.begin(),t_change.end());


double risk_surv = 0.0;

for (int i=0;i<= (int)(t_change.size()-2);i++){

    double lb=t_change[i]; // sub-interval to be considered [lb,ub)
    double ub=t_change[i+1];

    // int n_dt=0;

    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

        double beta_j, risk_surv_j;

        if (t_i_modified.at(xi_I_arg.at(j))<=lb & t_r_arg.at(xi_I_arg.at(j))>=ub) { 
  

            beta_j = para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

            if(para_other_CUPDATE.T>ub){
                risk_surv_j = beta_j*(ub-lb);
            }

            if(para_other_CUPDATE.T<lb){
                risk_surv_j = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval
            }

            if(para_other_CUPDATE.T>=lb & para_other_CUPDATE.T<=ub){
                double risk_surv_j_1 = beta_j*(para_other_CUPDATE.T-lb);
                double risk_surv_j_2 = (beta_j/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-1 ); 
                risk_surv_j = risk_surv_j_1 + risk_surv_j_2;
            }

            // risk_surv_j = (para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(ub-para_other_CUPDATE.T ))-exp(-para_current_arg.omega*(lb-para_other_CUPDATE.T )) ); // the integration of intensity function in the interval

            risk_surv = risk_surv + risk_surv_j;
        }


    } 

}


// risk_surv = risk_surv + para_current_arg.alpha*( t_change[t_change.size()-1] - t_change[0] );
risk_surv = risk_surv + para_current_arg.alpha*( para_other_CUPDATE.T - t_change[0] ) + (para_current_arg.alpha/(-para_current_arg.omega))*(exp(-para_current_arg.omega*(t_change[t_change.size()-1]-para_other_CUPDATE.T ))-1 ) ;

lh_square_modified.f_Surv = exp(-risk_surv);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_Surv); 

//----------

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_proposed - t_e_arg.at(subject_proposed), para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 


//----------

switch ( find(xi_R_arg.begin(), xi_R_arg.end(),subject_proposed) != (xi_R_arg.end()) ) { //return 1 when the subject is also in xi_R
case 1:{
    log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(subject_proposed)); //subtract part of likelihood that would be updated below
    // lh_square_modified.f_R.at(subject_proposed) = gsl_ran_weibull_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c, para_current_arg.d);
    lh_square_modified.f_R.at(subject_proposed) = gsl_ran_exponential_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c);
    log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(subject_proposed)); 
    break;
    }
    case 0:{
    break;
    }
}

//----------
switch ( find(xi_InR_arg.begin(), xi_InR_arg.end(),subject_proposed) != (xi_InR_arg.end()) ) { //return 1 when the subject is also in xi_InR
    case 1:{
    log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(subject_proposed)); //subtract part of likelihood that would be updated below
    // lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_weibull_P(t_max_CUPDATE - t_proposed, para_current_arg.c, para_current_arg.d);
    lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_exponential_P(t_max_CUPDATE - t_proposed, para_current_arg.c);
    log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(subject_proposed)); 
    
    break;
    }
    case 0:{
    break;
    }
}

//---------------
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) +(log_prior_y-log_prior_x)));

if ( isnan(log_lh_modified)==1 |  isinf(log_lh_modified)==1 ) acp_pr=0;


//      ofstream myfile_mcmc_out; 
//      myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
//      myfile_mcmc_out << acp_pr  << "," << t_i_arg.at(subject_proposed) << "," << t_proposed << endl;
//      myfile_mcmc_out.close();

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
    case 1: {
        lh_square_current_arg = lh_square_modified;
        // delta_mat_current_arg = delta_mat_modified;
        log_lh_current_arg = log_lh_modified;
        t_i_arg= t_i_modified;
    break;
    }
    
    case 0: {
    break;
    }
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::t_i_update_with_obs(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, vector<double>& t_i_arg, vector<double>& t_obs_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int subject_proposed){ // propose t_i using t_obs 

double t_proposed; // new t_i to be proposed
double t_low, t_up;

double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
// vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_i_modified = t_i_arg;


const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,subject_proposed); // set a seed

// int subject_proposed = xi_I_arg.at(gsl_rng_uniform_int (r_c, xi_I_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn



//--
vector<int> infecting_list;
vector<double> t_e_list;
infecting_list.reserve(xi_E_minus_arg.size());
t_e_list.reserve(xi_E_minus_arg.size());

for (int j=0;j<=(int) (xi_E_minus_arg.size()-1);j++){ 
    if (infected_source_arg.at(xi_E_minus_arg.at(j))==subject_proposed) {
        infecting_list.push_back(xi_E_minus_arg.at(j));
        t_e_list.push_back(t_e_arg.at(xi_E_minus_arg.at(j)));
    }
}


//--

switch(infecting_list.size()>=1){
    case 1:{
        double min_t_e = *min_element(t_e_list.begin(),t_e_list.end());
        t_up =  min_t_e;
        double dt = t_r_arg.at(subject_proposed) - min_t_e;
        t_low = std::max(t_e_arg.at(subject_proposed), t_up-(20-dt));  //make sure the biggest difference between t_i and t_r is  20


// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  dt << endl;
// myfile_mcmc_out.close();

      
    break;
    }
    case 0:{
        t_up =  min(t_r_arg.at(subject_proposed), t_max_CUPDATE);
        t_low = std::max(t_e_arg.at(subject_proposed), t_up-20);

    break;
    }
}
//---


double sigma=0;
double t_obs = t_obs_arg.at(subject_proposed);

if((t_up-t_obs)<=7) sigma=1;
if((t_up-t_obs)>7 & (t_up-t_obs)<=14) sigma=1;
if((t_up-t_obs)>14) sigma=1;

t_proposed = t_obs + gsl_ran_gaussian (r_c, sigma);

switch (t_proposed<=t_low | t_proposed>=t_up ) {
case 1: {
t_proposed = t_i_arg.at(subject_proposed); // remain as current value
break;
}
case 0: {
t_proposed = t_proposed;
break;
}
}


//---


t_i_modified.at(subject_proposed) = t_proposed;


double log_q_x, log_q_y;
log_q_x = log(gsl_ran_gaussian_pdf (t_i_arg.at(subject_proposed)-t_obs, sigma));
log_q_y = log(gsl_ran_gaussian_pdf (t_proposed -t_obs, sigma));


double log_prior_y =0.0;
double log_prior_x=0.0;

// log_prior_y = log(gsl_ran_gamma_pdf(t_proposed, t_o/pow(0.5,2.0), pow(0.5,2.0)));
// log_prior_x = log(gsl_ran_gamma_pdf(t_i_arg.at(subject_proposed), t_o/pow(0.5,2.0), pow(0.5,2.0)));
//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){


    // switch ((t_e_arg.at(xi_E_minus_arg.at(i))>=t_proposed) ) {

    // case 1:{


        double num_infectious = 0.0;

        log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

        // lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
        lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
        
        for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
       
            // double distance = func_distance (coordinate_CUPDATE[xi_I_arg.at(j)][0],coordinate_CUPDATE[xi_I_arg.at(j)][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);

            if (t_i_modified.at(xi_I_arg.at(j))<=t_e_arg.at(xi_E_minus_arg.at(i)) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i))) { 

                // if (t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(xi_E_minus_arg.at(i)))  {
                num_infectious = num_infectious + 1.0;
                // }

                // double delta_t =0;
                // switch (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_minus_arg.at(i))) {
                // case 1:{ // not yet recovered at e_i
                // delta_t = t_e_arg.at(xi_E_minus_arg.at(i)) - t_i_arg.at(xi_I_arg.at(j));
                // break;
                // }
                // case 0:{ // recovered before e_i
                // delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
                // break;
                // }    
                // }

         

                // lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) + beta_age_modified.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_t; // update kt_sum_E
                lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j))); // update kt_sum_E

            } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


        } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


        switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
            case 9999:{ // by background (this part is redundant now but retained for future developments)
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = (1.0/(num_infectious+1))*para_current_arg.alpha/557;
            break;
            }

            default :{ // not by background
            double distance= func_distance (coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][0],coordinate_CUPDATE[infected_source_arg.at(xi_E_minus_arg.at(i))][1],coordinate_CUPDATE[xi_E_minus_arg.at(i)][0],coordinate_CUPDATE[xi_E_minus_arg.at(i)][1]);
            lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = (1.0/(num_infectious+1))*para_current_arg.beta_age.at(age_gp_CUPDATE.at(infected_source_arg.at(xi_E_minus_arg.at(i))))*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2)/distance;
            break;
            }

        }

        // lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
        lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));

        lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

        lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));

        log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    // break;
    // }

    // case 0:{
    // break;
    // }
    // }
}
}


//----------

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_proposed - t_e_arg.at(subject_proposed), para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 


//----------

switch ( find(xi_R_arg.begin(), xi_R_arg.end(),subject_proposed) != (xi_R_arg.end()) ) { //return 1 when the subject is also in xi_R
case 1:{
    log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(subject_proposed)); //subtract part of likelihood that would be updated below
    // lh_square_modified.f_R.at(subject_proposed) = gsl_ran_weibull_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c, para_current_arg.d);
    lh_square_modified.f_R.at(subject_proposed) = gsl_ran_exponential_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c);
    log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(subject_proposed)); 
    break;
    }
    case 0:{
    break;
    }
}

//----------
switch ( find(xi_InR_arg.begin(), xi_InR_arg.end(),subject_proposed) != (xi_InR_arg.end()) ) { //return 1 when the subject is also in xi_InR
    case 1:{
    log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(subject_proposed)); //subtract part of likelihood that would be updated below
    // lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_weibull_P(t_max_CUPDATE - t_proposed, para_current_arg.c, para_current_arg.d);
    lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_exponential_P(t_max_CUPDATE - t_proposed, para_current_arg.c);
    log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(subject_proposed)); 
    
    break;
    }
    case 0:{
    break;
    }
}

//---------------
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) +(log_q_y-log_q_x)));

//      ofstream myfile_mcmc_out; 
//      myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
//      myfile_mcmc_out << acp_pr  << "," << t_i_arg.at(subject_proposed) << "," << t_proposed << endl;
//      myfile_mcmc_out.close();

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
    case 1: {
        lh_square_current_arg = lh_square_modified;
        // delta_mat_current_arg = delta_mat_modified;
        log_lh_current_arg = log_lh_modified;
        t_i_arg= t_i_modified;
    break;
    }
    
    case 0: {
    break;
    }
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::source_t_e_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, vector<int>& infected_source_arg, vector<double>& infected_distance_arg, int iter){

int subject_proposed ;
double t_proposed = 0.0;
double acp_pr = 0.0;

double distance_modified = -99;

double log_pr_forward=0.0; 
double log_pr_backward=0.0;

double log_pr_t_e_forward=0.0; 
double log_pr_t_e_backward=0.0;

double log_pr_source_y = 0.0;
double log_pr_source_x = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
// vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
// vector <int> index_modified = index_arg;
// vector <int> xi_E_minus_modified = xi_E_minus_arg;

// vector <int> xi_U_modified = xi_U_arg;
// vector <int> xi_E_modified = xi_E_arg;
// vector <int> xi_EnI_modified = xi_EnI_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed


    
subject_proposed = xi_E_minus_arg.at(gsl_rng_uniform_int (r_c, xi_E_minus_arg.size()));


// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("subject_proposed_replct.txt")).c_str(),ios::app);         
// myfile_mcmc_out << subject_proposed  << endl;
// myfile_mcmc_out.close();      

int source_x = infected_source_arg.at(subject_proposed);



//-------------- propose t_e before proposing a source  --------------//



double t_up =  std::min(t_i_arg.at(subject_proposed), t_max_CUPDATE);            

double t_low = std::max(t_i_arg.at(index_arg.at(0)), t_up -30);



int source_y;

vector<int> source_pool; // vector contains the indices of possible source for the subject


int num_pool =0;

while(num_pool<1){

    source_pool.clear();

    t_proposed= gsl_ran_flat(r_c, t_low, t_up );

    t_e_modified.at(subject_proposed) = t_proposed;


    for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

    switch( (t_i_arg.at(xi_I_arg.at(i))<=t_proposed) & (t_r_arg.at(xi_I_arg.at(i))>t_proposed)){
            case 1:{
                source_pool.push_back(xi_I_arg.at(i));  
            break;
            }
            case 0:{
            break;
            }
        }
    }
        

    source_pool.insert(source_pool.begin(),9999);

    num_pool = (int)source_pool.size();


    if (num_pool>=1) source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_pool)); // uniformly choose a new source (including bg)
   
}





// if (num_pool>=1) source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_pool)); // uniformly choose a new source (including bg)
// if (num_pool<1) source_y = source_x;

//--

vector<int> source_pool_before; 

for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

switch( (t_i_arg.at(xi_I_arg.at(i))<=t_e_arg[subject_proposed]) & (t_r_arg.at(xi_I_arg.at(i))>t_e_arg[subject_proposed])){
        case 1:{
            source_pool_before.push_back(xi_I_arg.at(i));  
        break;
        }
        case 0:{
        break;
        }
    }
}
    

// source_pool_before.insert(source_pool_before.begin(),9999);

int num_pool_before = (int)source_pool_before.size();


log_pr_source_y = log(1.0/num_pool);
log_pr_source_x = log(1.0/num_pool_before);




/*
//-------------- propose a new source before proposing a t_e --------------//

int source_y;

vector<int> source_pool; // vector contains the indices of possible source for the subject

double t_bound = min(t_i_arg.at(subject_proposed), t_max_CUPDATE);

for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

//  switch( (t_i_arg.at(xi_I_arg.at(i))<t_e_subject) & (t_r_arg.at(xi_I_arg.at(i))>=t_e_subject)){
    switch( t_i_arg.at(xi_I_arg.at(i))<t_bound){
    // switch( t_i_arg.at(xi_I_arg.at(i))<t_bound & (t_bound-t_i_arg.at(xi_I_arg.at(i)))<=20){

        case 1:{
            source_pool.push_back(xi_I_arg.at(i));  
        break;
        }
        case 0:{
        break;
        }
    }
}
    

// source_pool.insert(source_pool.begin(),9999); // when include 9999

int num_pool = (int)source_pool.size();

//-----------------------------propose uniformly-------------------------------------------//

// source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_pool)); // uniformly choose a new source (including bg)


//-----propose according to infectious challenge--------------------//


vector<double> ic(num_pool);

// ic.at(0) = para_current_arg.alpha; // when include 9999

// switch(num_pool>=2){ // when include 9999
switch(num_pool>=1){ 

    case 1:{

        // for (int j=1;j<=(num_pool-1);j++){ // when include 9999
        for (int j=0;j<=(num_pool-1);j++){


        // double distance= func_distance (coordinate_CUPDATE[source_pool.at(j)][0],coordinate_CUPDATE[source_pool.at(j)][1],coordinate_CUPDATE[subject_proposed][0],coordinate_CUPDATE[subject_proposed][1]);

        // ic.at(j) = (1.0/(num_pool))*para_current_arg.beta_age.at(age_gp_CUPDATE.at(source_pool.at(j)))*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2)/distance;
        ic.at(j) = para_current_arg.beta_age.at(age_gp_CUPDATE.at(source_pool.at(j)));

        
        }

        double *P=&ic.at(0); // convert vector to array
        gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)ic.size(),P);
        int link= gsl_ran_discrete (r_c, g);
        gsl_ran_discrete_free (g);

        source_y = source_pool.at(link); // a new source
        log_pr_forward = log(ic.at(link));

        switch(source_x==9999){
            case 0:{

                // double distance= func_distance (coordinate_CUPDATE[source_x][0],coordinate_CUPDATE[source_x][1],coordinate_CUPDATE[subject_proposed][0],coordinate_CUPDATE[subject_proposed][1]);

                // double ic_source_x = (1.0/(num_pool))*para_current_arg.beta_age.at(age_gp_CUPDATE.at(source_x))*func_contact(distance, para_current_arg.k_1,para_current_arg.k_2)/distance;
                double ic_source_x = para_current_arg.beta_age.at(age_gp_CUPDATE.at(source_x));

                log_pr_backward =log(ic_source_x);
            break;
            }
            case 1:{
                // double ic_source_x =  (1.0/(num_pool))*para_current_arg.alpha;
                double ic_source_x =  para_current_arg.alpha;

                log_pr_backward =log(ic_source_x);
            break;
            }
        }

        break;
    }

    case 0:{ // only primary source from pool

        source_y = 9999;
        log_pr_forward = log(para_current_arg.alpha);

        double ic_source_x =  para_current_arg.alpha;
        log_pr_backward =log(ic_source_x);

    break;
    }
}



//--------------//

//-------- propose a new t_e---------------//
switch(source_y){

case 9999:{ // by background

    double t_up = std::min(t_i_arg.at(subject_proposed), t_max_CUPDATE);

    double t_low = std::max(t_i_arg.at(index_arg[0]), t_up-20);
    
    t_proposed= gsl_ran_flat(r_c,t_low, t_up );
    
    log_pr_t_e_forward = log(1.0/(t_up-t_low));

break;
}

default :{ // not by background

    double t_up = std::min( t_r_arg.at(source_y), std::min(t_i_arg.at(subject_proposed), t_max_CUPDATE) );            

    // double t_low = std::max(t_i_arg.at(source_y), t_up -20 );
    double t_low =t_i_arg.at(source_y);
    
    t_proposed= gsl_ran_flat(r_c, t_low, t_up );

    log_pr_t_e_forward = log(1.0/(t_up-t_low));

break;
}

}

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);         
// myfile_mcmc_out << t_proposed <<"," << t_e_arg.at(subject_proposed)  << endl;
// myfile_mcmc_out.close(); 


//------------//
t_e_modified.at(subject_proposed) = t_proposed;
//------------//

switch(source_x){

case 9999:{ // by background

    double t_up = std::min(t_i_arg.at(subject_proposed), t_max_CUPDATE);

    double t_low = std::max(t_i_arg.at(index_arg[0]), t_up-20);
                
    log_pr_t_e_backward = log(1.0/(t_up-t_low));

break;
}

default :{ // not by background

    double t_up = std::min( t_r_arg.at(source_x), min(t_i_arg.at(subject_proposed), t_max_CUPDATE) );            

    double t_low = std::max(t_i_arg.at(source_x), t_up -20 );

    log_pr_t_e_backward = log(1.0/(t_up-t_low));

break;
}

}


//--------------------------------------//

*/



// double total_beta_current = 0.0; // this changes when update t_e

// for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

//     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_arg.at(subject_proposed) & t_r_arg.at(xi_I_arg.at(j))>t_e_arg.at(subject_proposed)) { 
//         total_beta_current = total_beta_current + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));

//     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 

// } 

//---


log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below

log_lh_modified = log_lh_modified - log(lh_square_modified.f_Grid.at(subject_proposed)); //subtract part of likelihood that would be updated below
log_lh_modified = log_lh_modified - log(lh_square_modified.f_Arc.at(subject_proposed)); //subtract part of likelihood that would be updated below

// lh_square_modified.kt_sum_E.at(subject_proposed) =0.0;

// double total_beta_modified = 0.0; // this changes when update t_e

// for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){



//     if (t_i_arg.at(xi_I_arg.at(j))<=t_e_modified.at(subject_proposed) & t_r_arg.at(xi_I_arg.at(j))>t_e_modified.at(subject_proposed)) { 

//         total_beta_modified = total_beta_modified + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)));


//     } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


// } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)




switch(source_y){
    case 9999:{ // by background 

    // lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha;

    double alpha_t = func_time_alpha (para_current_arg.alpha, t_e_modified.at(subject_proposed), para_other_CUPDATE.T, para_current_arg.omega);

    lh_square_modified.g_E.at(subject_proposed) = alpha_t;


    lh_square_modified.f_Grid.at(subject_proposed) = (pop_CUPDATE[subject_proposed]/para_other_CUPDATE.total_pop_of_grid)/pow(para_other_CUPDATE.grid_size,2.0);
    lh_square_modified.f_Arc.at(subject_proposed) = 1.0;   

    // lh_square_modified.f_Grid.at(subject_proposed) =1.0/((para_other_CUPDATE.x_max-para_other_CUPDATE.x_min)*(para_other_CUPDATE.y_max-para_other_CUPDATE.y_min));
    // lh_square_modified.f_Arc.at(subject_proposed) = 1.0;   


    break;
    }

    default :{ // not by background
    distance_modified= func_distance (coordinate_CUPDATE[source_y][0],coordinate_CUPDATE[source_y][1],coordinate_CUPDATE[subject_proposed][0],coordinate_CUPDATE[subject_proposed][1]);

    double beta_t = func_time_beta (para_current_arg.beta_age.at(age_gp_CUPDATE.at(source_y)), t_e_modified.at(subject_proposed), para_other_CUPDATE.T, para_current_arg.omega);
    lh_square_modified.g_E.at(subject_proposed) = beta_t*func_contact(distance_modified, para_current_arg.k_1, para_current_arg.k_2);

    //--
    vector<set_points_struct> set_points= circle_line_intersections (coordinate_CUPDATE[source_y][0],coordinate_CUPDATE[source_y][1], distance_modified, para_other_CUPDATE.n_line, grid_lines_CUPDATE);

        switch(set_points.size()>=1){
            case 1:{// have at least on intersection point
                vector<segments_struct> segments= func_segments_attributes (set_points, pop_grid_CUPDATE, distance_modified, para_other_CUPDATE);
             
                double mass = 0.0;
              
                for (int i_segment=0; i_segment<=(segments.size()-1);i_segment++){
                    mass = mass + segments[i_segment].den*segments[i_segment].len;
                }

                lh_square_modified.f_Grid.at(subject_proposed) = pop_CUPDATE[subject_proposed]/mass;
                lh_square_modified.f_Arc.at(subject_proposed) = 1.0;

            break;
            }    

            case 0:{

            lh_square_modified.f_Grid.at(subject_proposed) = 1.0/(2.0*M_PI*distance_modified);
            lh_square_modified.f_Arc.at(subject_proposed) = 1.0;


            break;
            }

        }


    break;
    }
    //---

}



// lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i))*para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.eta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
// lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha + lh_square_modified.kt_sum_E.at(subject_proposed);

// lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);

// lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed);


log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
log_lh_modified = log_lh_modified + log(lh_square_modified.f_Grid.at(subject_proposed)); //subtract part of likelihood that would be updated below
log_lh_modified = log_lh_modified + log(lh_square_modified.f_Arc.at(subject_proposed)); //subtract part of likelihood that would be updated below


        
//--------------------//

switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_I.txt")).c_str(),ios::app);
// myfile_mcmc_out   <<  subject_proposed << "," << t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_I.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------

switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_EnI.txt")).c_str(),ios::app);
// myfile_mcmc_out  <<    subject_proposed << "," <<t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_EnI.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------
        
// acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_backward-log_pr_forward) ));      
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_backward-log_pr_forward) + (log_pr_source_x - log_pr_source_y) ));      


if (isnan(exp( (log_lh_modified-log_lh_current_arg) + (log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_backward-log_pr_forward) + (log_pr_source_x - log_pr_source_y) ))==1 | isinf(exp( (log_lh_modified-log_lh_current_arg) + (log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_backward-log_pr_forward) + (log_pr_source_x - log_pr_source_y) ))==1) acp_pr=0;

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
// delta_mat_current_arg = delta_mat_modified;
log_lh_current_arg = log_lh_modified;
t_e_arg= t_e_modified;
infected_source_arg.at(subject_proposed) =  source_y;
// index_arg = index_modified;
// xi_E_minus_arg = xi_E_minus_modified;

infected_distance_arg.at(subject_proposed) = distance_modified;

break;
}

case 0: {
break;
}
}




// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("0000.txt")).c_str(),ios::app);  
// if (func_distance (coordinate_CUPDATE[source_y][0],coordinate_CUPDATE[source_y][1],coordinate_CUPDATE[subject_proposed][0],coordinate_CUPDATE[subject_proposed][1])==0){       
// myfile_mcmc_out << coordinate_CUPDATE[source_y][0] << ","<< coordinate_CUPDATE[source_y][1]  <<"," << coordinate_CUPDATE[subject_proposed][0]<< "," << coordinate_CUPDATE[subject_proposed][1] << endl;
// }
// myfile_mcmc_out.close(); 



gsl_rng_free(r_c);

}


/*------------------------------------------------*/

void mcmc_UPDATE::t_e_replct(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg,int iter){

int subject_proposed ;
double t_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

vector <int> xi_U_modified = xi_U_arg;
vector <int> xi_E_modified = xi_E_arg;
vector <int> xi_EnI_modified = xi_EnI_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed




	
subject_proposed = xi_E_arg.at(gsl_rng_uniform_int (r_c, xi_E_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn



//ofstream myfile_mcmc_out; 
//myfile_mcmc_out.open((string(path4)+string("subject_proposed_replct.txt")).c_str(),ios::app); 		//myfile_mcmc_out << subject_proposed  << endl;
//myfile_mcmc_out.close();		



// t_proposed= gsl_ran_flat(r_c, 0.0, min( t_i_arg.at(subject_proposed), t_max_CUPDATE) );
t_proposed= gsl_ran_flat(r_c, std::max(0.0,t_i_arg.at(subject_proposed)-20) , min( t_i_arg.at(subject_proposed), t_max_CUPDATE) );



//t_proposed= t_e_arg.at(subject_proposed) + 1.0*gsl_ran_gaussian(r_c,1.0); 
// double ub = min(t_max_CUPDATE,t_i_arg.at(subject_proposed) ) ; // upper bound, note: t_i is an extreme value when it has not gone through class I
// switch((t_proposed>ub)|(t_proposed<0.0)) { 
// case 1: {
// if (t_proposed>ub) t_proposed = ub - ( (t_proposed - ub) - (ub-0.0)*floor((t_proposed-ub)/(ub-0.0)) ); //reflection at ub with period ub-0
// if (t_proposed<0.0)  t_proposed =  0.0 + ( (0.0 - t_proposed) - (ub-0.0)*floor((0.0-t_proposed)/(ub-0.0)) ); //reflection at 0 with period ub-0
// break;
// }
// case 0: {
// break;
// }
// }


//t_proposed = t_e_arg.at(index_arg.at(0)) +10.0;


t_e_modified.at(subject_proposed) = t_proposed;
// ofstream myfile_mcmc; 
// myfile_mcmc.open((string(path4)+string("t_proposed_current_index.txt")).c_str(),ios::app);
// myfile_mcmc <<  t_proposed << ","<<index_arg.at(0) << ","<<t_e_arg.at(index_arg.at(0)) << endl;
// myfile_mcmc.close();

//-----------

lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;

for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again

if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {

switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
case 1:{
delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j)) ;
break;
}
case 0:{
delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
break;
}
}//end switch

lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

}
}

//----------
		
switch (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) { // return 1 if proposed subject not one of the original indexes

case 1:{

	switch(t_proposed<t_e_arg.at(index_arg.at(0))){ // original indexes would be replace by the chosen subject
	
	case 1:{

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("test_1_1.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
// myfile_mcmc_out.close();



	index_modified.clear();	
	index_modified.assign(1,subject_proposed);// replace index 
	xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));

	for ( int i =0; i<= (int) (index_arg.size()-1); i++){
	xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
	}	

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
	
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.q_E.at(subject_proposed)=0.0;
	lh_square_modified.g_E.at(subject_proposed)=1.0;
	lh_square_modified.h_E.at(subject_proposed)=1.0;
	lh_square_modified.f_E.at(subject_proposed)=1.0;
	

	
	for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
	
	//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
	
	lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(index_arg.at(i)); // secondary infection plays no role as the infectious times of others are always greater than the infection time of the original index
	lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i))*para_current_arg.eta*stb_arg.at(index_arg.at(i));
	lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
	lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
	}
	
	break;
	}
	
	
	case 0:{
	
	if (t_proposed==t_e_arg.at(index_arg.at(0))){ // addtion of one more index

	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_1_0_a.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
	//myfile_mcmc_out.close();
	
		//if (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()){ // return 1 if proposed subject not one of the original indexes
	
		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
		index_modified.push_back(subject_proposed); // add index 
		xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
		lh_square_modified.q_E.at(subject_proposed)=0.0;
		lh_square_modified.g_E.at(subject_proposed)=1.0;
		lh_square_modified.h_E.at(subject_proposed)=1.0;
		lh_square_modified.f_E.at(subject_proposed)=1.0;
	
		//}
	
	}
	
	if (t_proposed>t_e_arg.at(index_arg.at(0))){ // no shift of cases between xi_E and xi_E_minus



	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_1_0_b.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
	//myfile_mcmc_out.close();

	
		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below


	
		lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
		lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
		lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
	
		log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 


	} // end if t_proposs>t_e_arg.at()
	
	break;
	}
	
	}


break;
}


case 0: { // when chosen subject is one of the indexes



	index_modified.clear();

	int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
	double min_t = t_e_modified.at(first_min); // the minimum time of exposure
	
	int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time

// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("test_0.txt")).c_str(),ios::app);
// 	myfile_mcmc_out << first_min << "," << subject_proposed << "," <<  t_proposed << "," << num_min << endl;
// 	myfile_mcmc_out.close();
	
	switch (num_min>1) {
	case 1: {
	index_modified.reserve(n_CUPDATE);	
	for (int i=0; i<=(n_CUPDATE-1);i++){
	if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
	}
	break;
	}
	case 0:{
	index_modified.assign(1,first_min);
	break;
	}
	
	}

	xi_E_minus_modified = xi_E_arg;


	for (int i=0;i<= (int) (index_modified.size()-1); i++){

	xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i))); // this subject would have to be removed from likelihood function ( new index might be orginally an index, but the the log(lh_square_modified.f_E.at(index_modified.at(i)) will be zero in this case)
	//lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	//lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	//lh_square_modified.q_E.at(subject_proposed)=0.0;
	//lh_square_modified.g_E.at(subject_proposed)=1.0;
	//lh_square_modified.h_E.at(subject_proposed)=1.0;
	//lh_square_modified.f_E.at(subject_proposed)=1.0;

	lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
	lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
	lh_square_modified.q_E.at(index_modified.at(i))=0.0;
	lh_square_modified.g_E.at(index_modified.at(i))=1.0;
	lh_square_modified.h_E.at(index_modified.at(i))=1.0;
	lh_square_modified.f_E.at(index_modified.at(i))=1.0;			
	}

	switch(find(index_modified.begin(),index_modified.end(),subject_proposed) ==index_modified.end() ){ //return 1 when the chosen  subject is NO longer an index
	case 1:{
/*
	ofstream myfile_mcmc_out; 
	myfile_mcmc_out.open((string(path4)+string("test_0_1.txt")).c_str(),ios::app);
	myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(index_modified.at(0)) << endl;
	myfile_mcmc_out.close();*/

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 	

	lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
	lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
	lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
	lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);

	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
	
	break;
	}
	case 0:{

	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_0_0.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(index_modified.at(0)) << endl;
	//myfile_mcmc_out.close();

	break;
	}
	}

break;
}

}
		
		
//--------------------//

switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_I.txt")).c_str(),ios::app);
// myfile_mcmc_out   <<  subject_proposed << "," << t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_I.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------

switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_EnI.txt")).c_str(),ios::app);
// myfile_mcmc_out  <<    subject_proposed << "," <<t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_EnI.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------
		
acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));

//acp_pr =1.0;

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
delta_mat_current_arg = delta_mat_modified;
log_lh_current_arg = log_lh_modified;
t_e_arg= t_e_modified;
index_arg = index_modified;
xi_E_minus_arg = xi_E_minus_modified;
break;
}

case 0: {
break;
}
}

//gsl_rng_free(r_c);

//--------------
		
		//ofstream myfile_mcmc_out; 
		// 
		// myfile_mcmc_out.open((string(path4)+string("t_e_proposed_acp.txt")).c_str(),ios::app);
		// myfile_mcmc_out <<  subject_proposed  << "," <<  acp_pr  <<   "," << t_proposed <<   ","  << t_e_arg.at(subject_proposed) << endl;
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_modified.txt")).c_str(),ios::app);
		// if (index_modified.empty()==0){
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out <<  index_modified.at(i) << "," <<  t_e_modified.at(index_modified.at(i)) <<endl;
		// }
		// }
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg.txt")).c_str(),ios::app);
		// if (index_arg.empty()==0){
		// for (int i=0; i<=(int)(index_arg.size()-1); i++){
		// myfile_mcmc_out << index_arg.at(i) << "," << t_e_arg.at(index_arg.at(i)) <<endl;
		// }
		// }
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("find_xi_E_mibus.txt")).c_str(),ios::app);
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out << 1*(find(xi_E_minus_modified.begin(), xi_E_minus_modified.end(),index_modified.at(i))==xi_E_minus_modified.end()) << endl; // should always equals to 1, i.e., new index has been excluded from xi_E_minus
		// }
		// myfile_mcmc_out.close();
		// 

		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_kt_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.kt_sum_E.at(index_arg.at(0)) << endl; // shouls always be 0
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_k_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.k_sum_E.at(index_arg.at(0)) << endl;  // shouls always be 0
		// myfile_mcmc_out.close();

// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_replct.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << endl;
// 		myfile_mcmc_out.close();

		//---------------

gsl_rng_free(r_c);

}


/*------------------------------------------------*/

void mcmc_UPDATE::t_e_add_del(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg,const vector<int>& infected_source_arg, int iter){

int subject_proposed ;
double t_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

vector <int> xi_U_modified = xi_U_arg;
vector <int> xi_E_modified = xi_E_arg;
vector <int> xi_EnI_modified = xi_EnI_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed

//---------------------


//int action = gsl_rng_uniform_int (r_c, 3);// a int random number in [0,2] will be drawn : 0= replacement of an exsiting infection time; 1=deletion of an infection into gp of xi_EnI (choose from the gp of xi_U); 2= addition of an infection from gp of xi_EnI (add back  to xi_U)

// if (xi_EnI_arg.empty()==1){ // return 1 if xi_EnI is empty
// action = gsl_rng_uniform_int (r_c, 2); // only replacement OR addition is possible
// switch(action){
// case 0:{
// action = 0; //replacement
// break;
// }
// case 1:{
// action = 2; //addition
// break;
// }
// }
// }
// 
// if (xi_U_arg.empty()==1){ // return 1 if xi_U is empty
// action = gsl_rng_uniform_int (r_c, 2);// only replacement OR deletion is possible
// }



double P [2] ={2.0,2.0};
gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (2,P);
int action = gsl_ran_discrete (r_c, g); // a int random number in [0,1] will be drawn according to the weights in P :  0=deletion of an infection into gp of xi_EnI (choose from the gp of xi_U); 1= addition of an infection from gp of xi_EnI (add back  to xi_U)
gsl_ran_discrete_free (g);

if (xi_EnI_arg.empty()==1){ // return 1 if xi_EnI is empty
action = 1; // only addtion possible
}

if (xi_U_arg.empty()==1){ // return 1 if xi_U is empty
action = 0; // only deletion is possible
}


// 	ofstream myfile1_mcmc_out; 
// 	myfile1_mcmc_out.open((string(path4)+string("action.txt")).c_str(),ios::app);
// 	myfile1_mcmc_out << action << endl;
// 	myfile1_mcmc_out.close();


switch (action) {

case 0:{ // deletion
	subject_proposed = xi_EnI_arg.at(gsl_rng_uniform_int (r_c, xi_EnI_arg.size()));
	
// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("subject_proposed_del_before.txt")).c_str(),ios::app);
// 	myfile_mcmc_out <<endl;
// 	myfile_mcmc_out << "sub_del"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 	myfile_mcmc_out << subject_proposed  <<"," <<lh_square_current_arg.kt_sum_U.at(subject_proposed)<<"," <<lh_square_current_arg.f_U.at(subject_proposed)<< "," <<lh_square_current_arg.f_E.at(subject_proposed)<< "," <<lh_square_current_arg.f_EnI.at(subject_proposed)<<endl;
// 	myfile_mcmc_out.close();
	
	xi_EnI_modified.erase(find(xi_EnI_modified.begin(),xi_EnI_modified.end(),subject_proposed));
	xi_E_modified.erase(find(xi_E_modified.begin(),xi_E_modified.end(),subject_proposed));
	xi_U_modified.push_back(subject_proposed);
	
	//----
	
	t_e_modified.at(subject_proposed) = unassigned_time_CUPDATE;
	
	lh_square_modified.kt_sum_U.at(subject_proposed)=0.0;
	
	for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E 
	
		//if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
		
		switch (t_r_arg.at(xi_I_arg.at(j))>=t_max_CUPDATE) {
			case 1:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
			//lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)] ;
			break;
			}
			case 0:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
			break;
			}
		}//end switch
		
		lh_square_modified.kt_sum_U.at(subject_proposed) = lh_square_modified.kt_sum_U.at(subject_proposed) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));
		//}
	}
	
	
	//-----
	
	//log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(subject_proposed));  
	
	// lh_square_modified.q_T.at(subject_proposed) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_U.at(subject_proposed);
    lh_square_modified.q_T.at(subject_proposed) = para_current_arg.alpha*t_r_arg.at(subject_proposed)*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_U.at(subject_proposed);

	lh_square_modified.f_U.at(subject_proposed) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(subject_proposed),1.0);
	
	log_lh_modified = log_lh_modified  + log(lh_square_modified.f_U.at(subject_proposed));  //f_U becomes =/=1
	
	//----
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed));  //f_EnI becomes 1
	lh_square_modified.f_EnI.at(subject_proposed) = 1.0;
	
	//---
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));  //f_E becomes 1
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.q_E.at(subject_proposed)=0.0;
	lh_square_modified.g_E.at(subject_proposed)=1.0;
	lh_square_modified.h_E.at(subject_proposed)=1.0;
	lh_square_modified.f_E.at(subject_proposed)=1.0;
	
		switch(find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) {// returns 1 if it was NOT an index originally
		case 0:{// was an index
	
			index_modified.clear();
		
			int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
			double min_t = t_e_modified.at(first_min); // the minimum time of exposure
			
			int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time
			
				switch (num_min>1) {
				case 1: {
				index_modified.reserve(n_CUPDATE);	
				for (int i=0; i<=(n_CUPDATE-1);i++){
				if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
				}
				break;
				}
				case 0:{
				index_modified.assign(1,first_min);
				break;
				}
				
				}
	
			//xi_E_minus_modified = xi_E_arg;
		
			for (int i=0;i<= (int) (index_modified.size()-1); i++){
			xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));
			log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i)));
			lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
			lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
			lh_square_modified.q_E.at(index_modified.at(i))=0.0;
			lh_square_modified.g_E.at(index_modified.at(i))=1.0;
			lh_square_modified.h_E.at(index_modified.at(i))=1.0;
			lh_square_modified.f_E.at(index_modified.at(i))=1.0;	
			}
					
		break;
		}	
	
		case 1:{// was not an index
			xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
		break;
		}
		}

		//ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed_del_after.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<endl;
// 		myfile_mcmc_out << "sub_del"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 		myfile_mcmc_out << subject_proposed  <<"," <<lh_square_modified.kt_sum_U.at(subject_proposed)<<"," <<lh_square_modified.f_U.at(subject_proposed)<< "," <<lh_square_modified.f_E.at(subject_proposed)<< "," <<lh_square_modified.f_EnI.at(subject_proposed)<<endl;
// 		myfile_mcmc_out.close();

		//-----------------------------------------//
				
		acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(2.0/2.0)*(double)(xi_EnI_arg.size())*(1/t_max_CUPDATE)/((double)(xi_U_arg.size())+1.0) );

// 		myfile_mcmc_out.open((string(path4)+string("acp_del.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "acp_pr" << "," <<"log_lh_modified" << "," << "log_lh_unmodified"  << endl;
// 		myfile_mcmc_out <<  acp_pr << "," <<log_lh_modified << "," << log_lh_current_arg  << endl;
// 		myfile_mcmc_out.close();

		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		lh_square_current_arg = lh_square_modified;
		delta_mat_current_arg = delta_mat_modified;
		log_lh_current_arg = log_lh_modified;
		t_e_arg= t_e_modified;
		index_arg = index_modified;
		xi_E_minus_arg = xi_E_minus_modified;
		xi_U_arg = xi_U_modified;
		xi_E_arg = xi_E_modified;
		xi_EnI_arg = xi_EnI_modified;
		break;
		}
		
		case 0: {
		break;
		}
		}

// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_del.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "log_lh_current_arg"  << "," <<"log_lh_modified" << endl;
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << endl;
// 		myfile_mcmc_out.close();


break;
}// end of deletion

//------------------------------

case 1:{ // addition
	
	subject_proposed = xi_U_arg.at(gsl_rng_uniform_int (r_c, xi_U_arg.size()));
	
// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("subject_proposed_add_before.txt")).c_str(),ios::app);
// 	myfile_mcmc_out <<endl;
// 	myfile_mcmc_out << "sub_add"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 	myfile_mcmc_out << subject_proposed  <<"," <<lh_square_current_arg.kt_sum_U.at(subject_proposed)<<"," <<lh_square_current_arg.f_U.at(subject_proposed)<< "," <<lh_square_current_arg.f_E.at(subject_proposed)<< "," <<lh_square_current_arg.f_EnI.at(subject_proposed)<<endl;
// 	myfile_mcmc_out.close();
	
	xi_U_modified.erase(find(xi_U_modified.begin(),xi_U_modified.end(),subject_proposed)); // alter the gps xi_U, xi_E, xi_EnI
	xi_EnI_modified.push_back(subject_proposed);
	xi_E_modified.push_back(subject_proposed);
	
	//----
	
	t_proposed= gsl_ran_flat(r_c, 0.0, t_max_CUPDATE); // propose an infection time for the added infection
	
	t_e_modified.at(subject_proposed) = t_proposed;
	
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	
	for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E 
	
		if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
		
		switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
			case 1:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
			lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));
			break;
			}
			case 0:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
			break;
			}
		}//end switch
		
		lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j)) ;
		}
	}
	
	//--------
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(subject_proposed));  //f_U becomes 1
	lh_square_modified.q_T.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_U.at(subject_proposed)=0.0;
	lh_square_modified.f_U.at(subject_proposed)=1.0;
	
	//--------
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed));  // f_EnI =/= 1 (since it is from xi_U, it must be not in xi_I before addtion)
	
	lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 
	//-----
	
		switch (t_proposed<t_e_arg.at(index_arg.at(0))) { // 1 when t_proposed<index infection time
		
		case 1:{ // chosen one is the only index
	
	
			index_modified.clear();	
			index_modified.assign(1,subject_proposed);// replace index 
			//xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
		
			for ( int i =0; i<= (int) (index_arg.size()-1); i++){
			xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
			}	
		
	// 		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 		
	// 		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.q_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.g_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.h_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.f_E.at(subject_proposed)=1.0;
			
		
			for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
			
			//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
			
			lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(index_arg.at(i)); // secondary infection plays no role as the infectious times of others are always greater than the infection time of the original index
			lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i))*para_current_arg.eta*stb_arg.at(index_arg.at(i));
			lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
			lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
			}
	
		break;
		}
		
		case 0: {
			if((t_proposed==t_e_arg.at(index_arg.at(0))) ){ // add one more index
	
			index_modified.push_back(subject_proposed); // add index 
	
	// 		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
	// 		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.q_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.g_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.h_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.f_E.at(subject_proposed)=1.0;		
	
			// //xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
		
			}
	
			if((t_proposed>t_e_arg.at(index_arg.at(0))) ){
	
			xi_E_minus_modified.push_back(subject_proposed); 
	
			//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 
			
			lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
			lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(subject_proposed) + para_current_arg.eta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
			lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
			lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); 			
			}
	
		break;
		};
	
	}

		//ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed_add_after.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<endl;
// 		myfile_mcmc_out << "sub_add"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 		myfile_mcmc_out << subject_proposed  <<"," <<lh_square_modified.kt_sum_U.at(subject_proposed)<<"," <<lh_square_modified.f_U.at(subject_proposed)<< "," <<lh_square_modified.f_E.at(subject_proposed)<< "," <<lh_square_modified.f_EnI.at(subject_proposed)<<endl;
// 		myfile_mcmc_out.close();
		
		//-----------------------------------------//
		
		acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(2.0/2.0)*(double)(xi_U_arg.size())*(t_max_CUPDATE)/((double)(xi_EnI_arg.size())+1.0) );

// 		myfile_mcmc_out.open((string(path4)+string("acp_add.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "acp_pr" << "," <<"log_lh_modified" << "," << "log_lh_unmodified_arg"  << endl;
// 		myfile_mcmc_out <<  acp_pr << "," <<log_lh_modified << "," << log_lh_current_arg  << endl;
// 		myfile_mcmc_out.close();

		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		lh_square_current_arg = lh_square_modified;
		delta_mat_current_arg = delta_mat_modified;
		log_lh_current_arg = log_lh_modified;
		t_e_arg= t_e_modified;
		index_arg = index_modified;
		xi_E_minus_arg = xi_E_minus_modified;
		xi_U_arg = xi_U_modified;
		xi_E_arg = xi_E_modified;
		xi_EnI_arg = xi_EnI_modified;
		break;
		}
		
		case 0: {
		break;
		}
		}

// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_add.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "log_lh_current_arg"  << "," <<"log_lh_modified" << endl;
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << endl;
// 		myfile_mcmc_out.close();

break;
}// end of addition


} // end of switch on "action"


gsl_rng_free(r_c);

}


//-----------

// void mcmc_UPDATE::sample_source_distance (const vector< vector<double> >& kernel_mat_current_arg, const vector<int>& xi_E_arg, const vector<int>& xi_I_arg,  const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const para_key& para_current_arg, const para_aux& para_other, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_arg, int iter){

// ofstream myfile_out; 

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c,iter*1000); // set a seed


// // vector <double> u_imputed(xi_E_minus_arg.size()); // the imputed residuals for all infected, exclude index
// vector <int> source_imputed(xi_E_arg.size());
// vector <double> distance_imputed(xi_E_arg.size());


// for (int i=0; i<=((int)xi_E_arg.size()-1);i++){

//     vector<double> ic_link; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/
//     vector<double> pr_link; // probabilities of the infection links between infectious/primary infection and the new infected
    
//     vector<double> ic; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/susceptibles
//     vector<double> pr; // probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
//     vector<double> cdf; // cumulative probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
    
//     double total_ic_link; // the total instantaneous infective challenge to the new infected
//     double total_ic; // the total instantaneous infective challenge to the new infected/susceptibles
    
//     int link_imputed; // link imputed
//     double ic_imputed; // ic of the imputed link
//     int rank_imputed; // the rank of p_imputed among all links; 0 means first element

//     ic_link.reserve(xi_I_arg.size()+2); // 2 more positions reserved: one for primary infection and one for faciliating computation

    
//     ic_link.push_back(0.0); // 1st element in ic
//     ic_link.push_back(para_current_arg.alpha*para_current_arg.eta*stb_arg.at(xi_E_arg.at(i))); // 2nd..

//     int source;
//     double distance=-9999;

//     vector<int> pool_source; // only include those potential infectors (2 elements fewer thh ic_link defined above)

//     total_ic_link = 0.0; 
//     for (int j=0; j<=((int)xi_I_arg.size()-1);j++){

//     // if ((t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_arg.at(i))) & (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_arg.at(i)))){
//     if ((t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_arg.at(i))) & (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_arg.at(i))) &(coordinate_CUPDATE[xi_I_arg.at(j)][0]!=coordinate_CUPDATE[xi_E_arg.at(i)][0])){
      
//         ic_link.push_back(para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*para_current_arg.eta*stb_arg.at(xi_E_arg.at(i))*kernel_mat_current_arg[xi_E_arg.at(i)][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j)));
//         total_ic_link = total_ic_link + para_current_arg.beta_age.at(age_gp_CUPDATE.at(xi_I_arg.at(j)))*para_current_arg.eta*stb_arg.at(xi_E_arg.at(i))*kernel_mat_current_arg[xi_E_arg.at(i)][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));

//         pool_source.push_back(xi_I_arg.at(j));

//     }

//     }

//     total_ic_link = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(xi_E_arg.at(i)) + total_ic_link;

//     pr_link.resize(ic_link.size());

//     pr_link.at(0) = 0.0;
//     pr_link.at(1) = para_current_arg.alpha*para_current_arg.eta*stb_arg.at(xi_E_arg.at(i))/total_ic_link; //pr it was infected by primary infection
    
//     for (int k=2; k<=(int)(pr_link.size()-1);k++){
//     pr_link.at(k) = ic_link.at(k)/total_ic_link;
//     }

//     double *P=&pr_link.at(0);
//     gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)pr_link.size(),P);
//     link_imputed= gsl_ran_discrete (r_c, g); // the link imputed
//     gsl_ran_discrete_free (g);


//     if (link_imputed==1) {
//         source = 9999; // background
//         distance = -9999;
//     }

//     if (link_imputed>1) {
//         source = pool_source.at(link_imputed-2); 
//         distance = -log(kernel_mat_current_arg[xi_E_arg.at(i)][source])/para_current_arg.k_1;
//         // distance = distance_mat_CUPDATE[xi_E_arg.at(i)][source];

//     }

//     source_imputed.at(i) =source;
//     distance_imputed.at(i) = distance;





// //     ic_imputed = ic_link.at(link_imputed);

// //     //-----
    
// //     ic.reserve(ic_link.size()+(xi_I_arg.size()+1)*para_other.n);

// //     total_ic = 0.0; // sum of ic of links from the bunch of infectious to ALL susceptible subject before t_e_arg.at(xi_E_minus_arg.at(i))

// //     for (int iu=0; iu<=((int)para_other.n-1);iu++){

// //         double total_ic_iu =0.0; // sum if ic of links from the bunch of infectious to the particular susceptible subject

// //         switch(t_e_arg.at(iu)>t_e_arg.at(xi_E_minus_arg.at(i))) { // return 1 if the subject is susceptible before t_e_arg.at(xi_E_minus_arg.at(i))note: this excludes the links connected to xi_E_minus.at(i), where the infection actually happened
// //         case 1:{
// //         //ic.push_back(0.0); 
// //         ic.push_back(para_current_arg.alpha); 
    
// //         for (int j=0; j<=((int)xi_I_arg.size()-1);j++){
// //         if ((t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_minus_arg.at(i))) & (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_minus_arg.at(i)))){ // this ensures we are using same bunch of infectious sources
// //         ic.push_back(para_current_arg.beta*stb_arg.at(iu)*kernel_mat_current_arg[iu][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j)));
// //         total_ic_iu = total_ic_iu + stb_arg.at(iu)*kernel_mat_current_arg[iu][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));
// //         }
// //         }
// //         total_ic_iu = para_current_arg.alpha + para_current_arg.beta*total_ic_iu; 
    
// //         break;
// //         }
    
// //         case 0:{
// //         break;
// //         }
// //         }

// //         total_ic = total_ic + total_ic_iu;

// //     }
    
// //     total_ic = total_ic + total_ic_link; // add the ic from the infectious to infected links

// //     ic.insert(ic.begin(),ic_link.begin(),ic_link.end()); // combine ic_link and ic

// //     pr.resize(ic.size());

// //     for (int k=0; k<=((int)pr.size()-1);k++){
// //     pr.at(k) = ic.at(k)/total_ic;
// //     }

// //     sort(pr.begin(),pr.end());

// //     double p_imputed = ic_imputed / total_ic;
    
// //     rank_imputed = distance(  pr.begin(), find(pr.begin(),pr.end(),p_imputed) );


// //     //---------//

// //     int num_dup=1; // number of prs same to p_imputed (including p_imputed itself)
// //     vector<double> dup_test = pr;
// //     dup_test.erase(dup_test.begin()+rank_imputed);

// //     switch(find(dup_test.begin(),dup_test.end(),p_imputed)==dup_test.end()) { // return 1 when no duplicate 
// //     case 0:{ // with dup
// //     num_dup = count(pr.begin(),pr.end(),p_imputed);
// //     break;
// //     }
// //     case 1:{
// //     num_dup=1;
// //     break;
// //     }
// //     }

// //     //---



// //     cdf.resize(pr.size());

// //     cdf.at(0) = 0.0;
// //     for (int k=1; k<=(int)(cdf.size()-1);k++){
// //     cdf.at(k) = pr.at(k) + cdf.at(k-1);
// //     }

// //     double u_i=0.0; // the u_i to be imputed
// //     switch(rank_imputed==0){
// //     case 1:{
// //     u_i = 0.0;
// //     break;
// //     }
// //     case 0:{
// // /*  u_i = gsl_ran_flat(r_c, cdf.at(rank_imputed-1),  cdf.at(rank_imputed-1)+p_imputed); // a uniform rv drawn between cdf.at(rank_imputed-1), and cdf.at(rank_imputed-1)+p_imputed*/
// //     u_i = gsl_ran_flat(r_c, cdf.at(rank_imputed-1),  cdf.at(rank_imputed-1)+num_dup*p_imputed); // a uniform rv drawn between cdf.at(rank_imputed-1), and cdf.at(rank_imputed-1)+num_dup*p_imputed
// //     break;
// //     }
// //     }

// //     u_imputed.at(i) =u_i;

// //     //------------------

// //     pr_link.clear();
// //     ic_link.clear();

// //     pr.clear();
// //     ic.clear();
// //     cdf.clear();

// }




//     myfile_out.open((string(path4)+string("source_imputed.txt")).c_str(),ios::app);
//     myfile_out << endl;
//     for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
//     if (i!= (int) xi_E_arg.size()-1) myfile_out << source_imputed.at(i) << ",";
//     if (i== (int) xi_E_arg.size()-1) myfile_out << source_imputed.at(i);
//     }
//     myfile_out.close();


//     myfile_out.open((string(path4)+string("distance_imputed.txt")).c_str(),ios::app);
//     myfile_out << endl;
//     for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
//     if (i!= (int) xi_E_arg.size()-1) myfile_out << distance_imputed.at(i) << ",";
//     if (i== (int) xi_E_arg.size()-1) myfile_out << distance_imputed.at(i);
//     }
//     myfile_out.close();


//     //---sort output according to t_e ---//

//     // vector<int> t_e_sub(xi_E_arg.size());
//     // for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
//     //     t_e_sub.at(i) = t_e_arg.at(xi_E_arg.at(i));
//     // }

//     // vector<paired_to_t_e> paired(xi_E_arg.size());
//     // for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
//     //     paired[i].source = source_imputed[i];
//     //     paired[i].distance = distance_imputed[i];
//     //     paired[i].t_e = t_e_sub[i];
//     // }


//     // std::sort(paired.begin(), paired.end(), by_t_e());



//     // myfile_out.open((string(path4)+string("source_imputed.txt")).c_str(),ios::app);
//     // myfile_out << endl;
//     // for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
//     // if (i!= (int) xi_E_arg.size()-1) myfile_out << paired[i].source << ",";
//     // if (i== (int) xi_E_arg.size()-1) myfile_out << paired[i].source;
//     // }
//     // myfile_out.close();


//     // myfile_out.open((string(path4)+string("distance_imputed.txt")).c_str(),ios::app);
//     // myfile_out << endl;
//     // for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
//     // if (i!= (int) xi_E_arg.size()-1) myfile_out << paired[i].distance << ",";
//     // if (i== (int) xi_E_arg.size()-1) myfile_out << paired[i].distance;
//     // }
//     // myfile_out.close();


//     //-------------//



// gsl_rng_free(r_c);

// }


// //-----------------------------------------------


