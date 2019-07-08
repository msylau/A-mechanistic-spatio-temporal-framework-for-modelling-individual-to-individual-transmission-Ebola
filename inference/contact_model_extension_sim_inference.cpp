
//include the term 1/r in the likelihood to reflect transformation (incorporated in func_contact)//


#include "contact_model_extension_sim_inference_header.h"





using namespace std;

int main (){

int n_iter = 200000; //number of iterations for MCMC


para_key para_true;

para_aux para_other;

epi_struct epi_final;

// grid_struct grid_data;

grid_lines_struct grid_lines;

vector<int> index;

IO_simpara (para_true, para_other);  //Importing aux/key parameters 

vector < vector<double> > coordinate(para_other.n,vector<double>(2));

vector < vector<double> > pop_grid (para_other.n_row_grid,vector<double>(para_other.n_col_grid));

IO_simdata(para_true, para_other, coordinate, epi_final,index, pop_grid, grid_lines); //Importing simulated/real data 

//--

// ofstream myfile_out_temp; 

// double r_tmp= 5.0;

// vector<set_points_struct> set_points_test = circle_line_intersections (0.3,5, r_tmp, para_other.n_line, grid_lines);

// int n_temp = set_points_test.size();

// myfile_out_temp.open((string(path4)+string("set_points_test.txt")).c_str(),ios::out);
// for (int i=0; i<=(n_temp-1);i++){
// // myfile_out_temp << set_points_test.coor_x[i] << "," << set_points_test.coor_y[i] << "," << set_points_test.theta[i] << endl;
// myfile_out_temp << set_points_test[i].coor_x << "," << set_points_test[i].coor_y << "," << set_points_test[i].theta << endl;

// }
// myfile_out_temp << "n_temp" << endl;
// myfile_out_temp << n_temp;
// myfile_out_temp.close();

// //--
// vector<segments_struct> segments_test = func_segments_attributes ( set_points_test, pop_grid, r_tmp, para_other);

// int m_temp = segments_test.size();

// double total_theta=0;
// double cirumference=0;

// myfile_out_temp.open((string(path4)+string("segments_test.txt")).c_str(),ios::out);
// for (int i=0; i<=(m_temp-1);i++){
// total_theta = total_theta + segments_test[i].theta_abs;
// cirumference = cirumference + segments_test[i].len ;
// myfile_out_temp << segments_test[i].len << "," << segments_test[i].den << "," << segments_test[i].theta_abs << endl;

// }
// myfile_out_temp << "m_temp" << endl;
// myfile_out_temp << m_temp<< endl;
// myfile_out_temp << "total_theta" << endl;
// myfile_out_temp << total_theta << endl;
// myfile_out_temp << "cirumference" << endl;
// myfile_out_temp << cirumference << "," <<  2*M_PI*r_tmp << endl;
// myfile_out_temp.close();


/*----------------------------*/

ifstream myfile_in;
ofstream myfile_out; 

vector<int> xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_InR, xi_null_onset, xi_onset, xi_onset_2; // indices sets indicating the individuals stay in S OR have gone through the other classes (E OR I OR R), and individuals hve gone through E but not I (EnI) and I but not R (InR)


xi_U.reserve(para_other.n); //dynamically updated if necessary
xi_E.reserve(para_other.n);
xi_E_minus.reserve(para_other.n);
xi_I.reserve(para_other.n);
xi_R.reserve(para_other.n);
xi_EnI.reserve(para_other.n);
xi_InR.reserve(para_other.n);
xi_null_onset.reserve(para_other.n); // set of cases without trustworthy onset times (lab tested postive)
xi_onset.reserve(para_other.n); // set of cases with trustworthy onset times (lab tested postive)
xi_onset_2.reserve(para_other.n); // set of cases with trustworthy onset times (lab tested postive); but need to imputed within 1 day range



for (int i=0; i<=(para_other.n-1);i++){
if (epi_final.t_e.at(i)==para_other.unassigned_time) xi_U.push_back(i);
if (epi_final.t_e.at(i)!=para_other.unassigned_time) xi_E.push_back(i);
if (epi_final.t_i.at(i)!=para_other.unassigned_time) xi_I.push_back(i);
if (epi_final.t_r.at(i)!=para_other.unassigned_time & epi_final.t_i.at(i)!=para_other.unassigned_time) xi_R.push_back(i);
if (epi_final.onset_type.at(i)==0) xi_null_onset.push_back(i);
if (epi_final.onset_type.at(i)==1) xi_onset.push_back(i);
if (epi_final.onset_type.at(i)==2) xi_onset_2.push_back(i);

}

int num_null_onset = xi_null_onset.size();
int num_onset = xi_onset.size();
int num_onset_2 = xi_onset_2.size();


xi_E_minus = xi_E;
for (int i=0; i<= (int)(index.size()-1);i++){
xi_E_minus.erase(find(xi_E_minus.begin(),xi_E_minus.end(),index.at(i)));
} // E set excluding index

xi_EnI = xi_E;
for (int i=0; i<= (int)(xi_I.size()-1);i++){
xi_EnI.erase(find(xi_EnI.begin(),xi_EnI.end(),xi_I.at(i)));
} // E set excluding I

xi_InR = xi_I;
for (int i=0; i<= (int)(xi_R.size()-1);i++){
xi_InR.erase(find(xi_InR.begin(),xi_InR.end(),xi_R.at(i)));
} // I set excluding R

myfile_out.open((string(path4)+string("xi_null_onset.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_null_onset.empty()!=1){
for (int i=0; i<=((int)xi_null_onset.size()-1);i++){
myfile_out << xi_null_onset.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_null_onset.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_onset_2.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_onset_2.empty()!=1){
for (int i=0; i<=((int)xi_onset_2.size()-1);i++){
myfile_out << xi_onset_2.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_onset_2.size();
myfile_out.close();


myfile_out.open((string(path4)+string("xi_onset.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_onset.empty()!=1){
for (int i=0; i<=((int)xi_onset.size()-1);i++){
myfile_out << xi_onset.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_onset.size();
myfile_out.close();


// myfile_out.open((string(path4)+string("xi_E_true.txt")).c_str(),ios::out);
// myfile_out << "k" << endl;
// if (xi_E.empty()!=1){
// for (int i=0; i<=((int)xi_E.size()-1);i++){
// myfile_out << xi_E.at(i) << endl;
// }
// }
// myfile_out << "size" << endl;
// myfile_out << xi_E.size();
// myfile_out.close();

myfile_out.open((string(path4)+string("xi_E_minus_true.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_E_minus.empty()!=1){
for (int i=0; i<=((int)xi_E_minus.size()-1);i++){
myfile_out << xi_E_minus.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_E_minus.size();
myfile_out.close();

// myfile_out.open((string(path4)+string("xi_U_true.txt")).c_str(),ios::out);
// myfile_out << "k" << endl;
// if (xi_U.empty()!=1){
// for (int i=0; i<=((int)xi_U.size()-1);i++){
// myfile_out << xi_U.at(i) << endl;
// }
// }
// myfile_out << "size" << endl;
// myfile_out << xi_U.size();
// myfile_out.close();

// myfile_out.open((string(path4)+string("xi_I_true.txt")).c_str(),ios::out);
// myfile_out << "k" << endl;
// if (xi_I.empty()!=1){
// for (int i=0; i<=((int)xi_I.size()-1);i++){
// myfile_out << xi_I.at(i) << endl;
// }
// }
// myfile_out << "size" << endl;
// myfile_out << xi_I.size();
// myfile_out.close();

// myfile_out.open((string(path4)+string("xi_EnI_true.txt")).c_str(),ios::out);
// myfile_out << "k" << endl;
// if (xi_EnI.empty()!=1){
// for (int i=0; i<=((int)xi_EnI.size()-1);i++){
// myfile_out << xi_EnI.at(i) << endl;
// }
// }
// myfile_out << "size" << endl;
// myfile_out << xi_EnI.size();
// myfile_out.close();

// myfile_out.open((string(path4)+string("xi_InR_true.txt")).c_str(),ios::out);
// myfile_out << "k" << endl;
// if (xi_InR.empty()!=1){
// for (int i=0; i<=((int)xi_InR.size()-1);i++){
// myfile_out << xi_InR.at(i) << endl;
// }
// }
// myfile_out << "size" << endl;
// myfile_out << xi_InR.size();
// myfile_out.close();

/*----------------------------*/

lh_SQUARE lh_square; // the object contains the information of the likelihood contribution of each individual, and it is dynamic and changed during MCMC sampling

lh_square.f_U.assign(para_other.n,1.0);
lh_square.q_T.assign(para_other.n,0.0);
lh_square.kt_sum_U.assign(para_other.n,0.0);
lh_square.f_E.assign(para_other.n,1.0);
lh_square.g_E.assign(para_other.n,1.0);
lh_square.h_E.assign(para_other.n,1.0);
lh_square.k_sum_E.assign(para_other.n,0.0);
lh_square.q_E.assign(para_other.n,0.0);
lh_square.kt_sum_E.assign(para_other.n,0.0);
lh_square.f_I.assign(para_other.n,1.0);
lh_square.f_R.assign(para_other.n,1.0);
lh_square.f_EnI.assign(para_other.n,1.0);
lh_square.f_InR.assign(para_other.n,1.0);

lh_square.f_Grid.assign(para_other.n,1.0);
lh_square.f_Arc.assign(para_other.n,1.0);

lh_square.f_Surv=1.0;

/*----------------------------*/


vector < vector<double> > kernel_mat(para_other.n, vector<double>(para_other.n)), delta_mat(para_other.n, vector<double>(para_other.n));; // a dynamic matrix contain the "kernel distance" between each individual (note:named distance_mat in simulation code)
vector <double> norm_const(para_other.n);


//--
// vector < vector<double> > distance_mat(para_other.n, vector<double>(para_other.n));

// for (int i=0;i<=(para_other.n-1);i++) {
//  for (int j=0;j<=(para_other.n-1);j++) {
//  if (i==j) distance_mat[i][j]=0.0;
//  if (i<j) distance_mat[i][j] = func_distance (coordinate[i][0],coordinate[i][1],coordinate[j][0],coordinate[j][1]);
//  if (i>j) distance_mat[i][j]=distance_mat[j][i];
//  }
// }

//---
// myfile_out.open((string(path4)+string("kernel_matrix_before.txt")).c_str(),ios::app);
// for (int i=0;i<=(para_other.n-1);i++) {
// for (int j=0;j<=(para_other.n-1);j++) {
// if (j<(para_other.n-1)) myfile_out << kernel_mat[i][j] << ",";
// if (j==(para_other.n-1)) myfile_out << kernel_mat[i][j] << " " << endl;
// }
// }
// myfile_out.close();

vector<double> t_e = epi_final.t_e, t_i=epi_final.t_i, t_r=epi_final.t_r; 

vector<int> infected_source = epi_final.infected_source;

// vector<int> gp_stb= epi_final.gp_stb;  
vector<double> stb = epi_final.stb;

myfile_out.open((string(path4)+string("stb_true.txt")).c_str(),ios::out);
for (int i=0; i<=((int)stb.size()-1);i++){
myfile_out  << stb.at(i)  << endl;
}
myfile_out.close();

//para_key para_current = para_true; 






// FUNC func;

// func.set_para(para_true, para_other, coordinate, xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_InR,t_e, t_i, t_r, stb, index);
// func.initialize_kernel_mat(kernel_mat, norm_const); // compute the kernel matrix
// func.initialize_delta_mat(delta_mat); // compute the kernel matrix

// func.initialize_lh_square(lh_square, kernel_mat, delta_mat, norm_const); //initialize lh_square

// double log_lh_true = log_lh_func (lh_square, para_other.n); // the log-likelihood value when using true values of the parameters


// myfile_out.open((string(path4)+string("kernel_matrix_after.txt")).c_str(),ios::out);
// for (int i=0;i<=(para_other.n-1);i++) {
// for (int j=0;j<=(para_other.n-1);j++) {
// if (j<(para_other.n-1)) myfile_out << kernel_mat[i][j] << ",";
// if (j==(para_other.n-1)) myfile_out << kernel_mat[i][j] << " " << endl;
// }
// }
// myfile_out.close();
/*
myfile_out.open((string(path4)+string("lh_square_q_E_Minus_epi_final_q.txt")).c_str(),ios::out);
//for (int i=0;i<=(para_other.n-1);i++) {
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.q_E.at(xi_E.at(i)) - epi_final.q.at(xi_E.at(i)) << "," <<endl;
}// it should give 0 for all entries, the small deviation seems to do with rounding error in simulation
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_f_U.txt")).c_str(),ios::out);

for (int i=0;i<=(para_other.n-1);i++) {
myfile_out << lh_square.f_U.at(i) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_f_E.txt")).c_str(),ios::out);

for (int i=0;i<=(para_other.n-1);i++) {
myfile_out << lh_square.f_E.at(i) <<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_g_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.g_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_q_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.q_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_h_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.h_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_k_sum_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.k_sum_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();

myfile_out.open((string(path4)+string("lh_square_kt_sum_E.txt")).c_str(),ios::out);
for (int i=0;i<= (int) (xi_E.size()-1);i++) {
myfile_out << lh_square.kt_sum_E.at(xi_E.at(i)) <<endl;

}
myfile_out.close();*/

// myfile_out.open((string(path4)+string("log_value_true.txt")).c_str(),ios::out);
// myfile_out << log_lh_true <<endl;
// myfile_out.close();

/*----------------------------*/

//lh_plot(para_true,para_other, coordinate, xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_InR, t_e, t_i, t_r, index); // this generates likelihood values at different values of parameters


/*-----------------------------------------------------Start of MCMC sampling------------------------------------------*/

para_key para_current = para_true; 
vector<int> xi_I_current = xi_I;
vector<int> xi_U_current = xi_U;
vector<int> xi_E_current = xi_E;
vector<int> xi_E_minus_current = xi_E_minus;
vector<int> xi_R_current= xi_R;
vector<int> xi_EnI_current = xi_EnI;
vector <int> xi_InR_current =xi_InR;
vector<double> t_e_current =t_e;
vector<double>t_i_current = t_i;
vector<double>t_r_current = t_r;
vector<int>index_current = index;

vector<double> t_obs = t_i; // note that only the entries for valid onset would be used in t_i_update_with_obs

vector<double> stb_current = stb;
// vector<int> gp_stb_current = gp_stb;

//----------------------------------------//
// stb_current.assign(para_other.n,1.0); // when ignoring the stb
//------------------------------------//



vector<int> infected_source_current = infected_source;

vector<double> infected_distance_current(para_other.n);
infected_distance_current.assign(para_other.n,-99);



for (int j=0; j<=((int)xi_E_minus_current.size()-1);j++){

if (infected_source_current.at(xi_E_minus_current.at(j))!=9999){	
	infected_distance_current.at(xi_E_minus_current.at(j)) = func_distance_for_main (coordinate[infected_source_current.at(xi_E_minus_current.at(j))][0],coordinate[infected_source_current.at(xi_E_minus_current.at(j))][1],coordinate[xi_E_minus_current.at(j)][0],coordinate[xi_E_minus_current.at(j)][1]);
}

}


//--


lh_SQUARE lh_square_current ;
vector < vector<double> > kernel_mat_current(para_other.n, vector<double>(para_other.n)), delta_mat_current(para_other.n, vector<double>(para_other.n)); // a dynamic matrix contain the "kernel distance" 
vector <double> norm_const_current(para_other.n);


// para_current.alpha = 0.0; //initialization of parameter to be estimated
// // para_current.beta = 0.2 ; //initialization of parameter to be estimated
// para_current.mu_lat = 5; //initialization of parameter to be estimated
// para_current.var_lat = 20; //initialization of parameter to be estimated
// para_current.c = 2; //initialization of parameter to be estimated
// para_current.d = 1; //initialization of parameter to be estimated
// para_current.k_1 = 0.3; //initialization of parameter to be estimated
// para_current.k_2 = 2; //initialization of parameter to be estimated

// para_current.beta_age.resize(5); //initialization of betas for different age gp
// para_current.beta_age.at(0) = 0.1;
// para_current.beta_age.at(1) = 0.5;
// // para_current.beta_age.at(2) = 0.5;
// // para_current.beta_age.at(3) = 0.5;
// // para_current.beta_age.at(4) = 0.5;

// para_current.omega=0.002;


// para_current.eta=1;



// //--- (exp kernel) set parameters value to be posterior means, when compute 2nd term in DIC4 and DIC8--//


// para_current.alpha = 0.0012; 
// para_current.beta = 3.0573 ; 
// para_current.mu_lat = 4.89;
// para_current.var_lat =  2.553;
// para_current.c = 1.999631;
// para_current.d = 2.05; 
// para_current.k_1 =  0.02985; 
// para_current.k_2 = 1.0;


// //--- (Cauchy kernel) set parameters value to be posterior means, when compute 2nd term in DIC4 and DIC8--//
// para_current.alpha =0.00064; 
// para_current.beta =  1.5839 ; 
// para_current.mu_lat = 4.71;
// para_current.var_lat =1.374;
// para_current.c = 1.99814;
// para_current.d = 2.050879; 
// para_current.k_1 = 1.0; 
// para_current.k_2 = 0.0396; // mode

// //--- (power-law kernel) set parameters value to be posterior means, when compute 2nd term in DIC4 and DIC8--//
// para_current.alpha = 0.0008;
// para_current.beta = 3.822;
// para_current.mu_lat = 4.97;
// para_current.var_lat =  2.93;
// para_current.c = 1.999662;
// para_current.d = 2.052737; 
// para_current.k_1 = 1.0; 
// para_current.k_2 = 2.928;


// //--- (exp kernel & exp latent) set parameters value to be posterior means, when compute 2nd term in DIC4 and DIC8--//
// para_current.alpha = 0.001105 
// para_current.beta = 2.8114;
// para_current.mu_lat = 4.438;
// para_current.var_lat =  2.5;
// para_current.c = 1.999403;
// para_current.d = 2.052137; 
// para_current.k_1 =  0.0313;
// para_current.k_2 = 1.0;


////---- intialization of xi_E and xi_EnI and xi_U ----////

xi_E_current = xi_I_current; // individuals gone through I (assumed known here) would be initialized as infected
xi_EnI_current.clear();
xi_U_current.clear();
for (int i=0; i<= (int)(para_other.n-1);i++){
if(find(xi_E_current.begin(),xi_E_current.end(),i)==xi_E_current.end()){
xi_U_current.push_back(i);
t_e_current.at(i)=para_other.unassigned_time;
}
}


myfile_out.open((string(path4)+string("xi_E_initial.txt")).c_str(),ios::out);
myfile_out << "k"  << endl;
if (xi_E_current.empty()!=1){
for (int i=0; i<=((int)xi_E_current.size()-1);i++){
myfile_out << xi_E_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_E_current.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_I_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_I_current.empty()!=1){
for (int i=0; i<=((int)xi_I_current.size()-1);i++){
myfile_out << xi_I_current.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_I_current.size();
myfile_out.close();



myfile_out.open((string(path4)+string("xi_EnI_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_EnI_current.empty()!=1){
for (int i=0; i<=((int)xi_EnI_current.size()-1);i++){
myfile_out << xi_EnI_current.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_EnI_current.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_U_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_U_current.empty()!=1){
for (int i=0; i<=((int)xi_U_current.size()-1);i++){
myfile_out << xi_U_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_U_current.size();
myfile_out.close();

////---- intialization of index cases and xi_E_minus and and time of exposure----////

// index_current.clear();
// xi_E_minus_current = xi_E_current;

// const gsl_rng_type* T= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r = gsl_rng_alloc (T); // r is pointer points to an object with Type T
// gsl_rng_set (r, 1000); // set a seed

// double min_t_i = *min_element(t_i_current.begin(),t_i_current.end()); // min of t_i (corresponding individual would be used as index, and min_t_i would be the lower bound of all infection times)
// int pst_min_t_i = distance(t_i_current.begin(),min_element(t_i_current.begin(),t_i_current.end())); // the individual with t_i=min_t_i

// myfile_out.open((string(path4)+string("min_t_i.txt")).c_str(),ios::app);
// myfile_out << pst_min_t_i << "," << min_t_i << endl;
// myfile_out.close();

// for (int i=0; i<=(int) (xi_E_current.size()-1);i++){

// 	switch(t_i_current.at(xi_E_current.at(i))==para_other.unassigned_time){
	
// 	case 1:{
	
// 	// t_e_current.at(xi_E_current.at(i)) = gsl_ran_flat(r,0.0, para_other.t_max); // when alpha=/=0.0
	
// 	t_e_current.at(xi_E_current.at(i)) = t_e_current.at(xi_E_current.at(i)); // when alpha=/=0.0

// 	break;
// 	}
	
// 	case 0:{
	
// 		t_e_current.at(xi_E_current.at(i)) = gsl_ran_flat(r,max(0.0,t_i_current.at(xi_E_current.at(i))- 5.0),t_i_current.at(xi_E_current.at(i))); // when alpha=/=0.0
	
		
// 	break;
// 	}


// 	}

// }


// double min_t = *min_element(t_e_current.begin(),t_e_current.end());
// for (int i=0; i<= (int)(xi_E_current.size()-1);i++){
// if (t_e_current.at(xi_E_current.at(i))==min_t) {
// index_current.push_back(xi_E_current.at(i));
// }
// }
// for (int i=0; i<= (int)(index_current.size()-1);i++){
// xi_E_minus_current.erase(find(xi_E_minus_current.begin(),xi_E_minus_current.end(),index_current.at(i)));
// }

// gsl_rng_free(r);




// myfile_out.open((string(path4)+string("initial_time_exposure.txt")).c_str(),ios::app);
// for (int i=0; i<= (int) (xi_E_minus_current.size()-1);i++){
// myfile_out << t_e_current.at(xi_E_minus_current.at(i)) << endl;
// }
// myfile_out.close();

// myfile_out.open((string(path4)+string("index_initial.txt")).c_str(),ios::app);
// for (int i=0; i<= (int) (index_current.size()-1);i++){
// myfile_out << index_current.at(i) << "," << t_e_current.at(index_current.at(i)) << endl;
// }
// myfile_out.close();

// myfile_out.open((string(path4)+string("xi_E_minus_initial.txt")).c_str(),ios::out);
// myfile_out << "k" << endl;
// if (xi_E_minus_current.empty()!=1){
// for (int i=0; i<=((int)xi_E_minus_current.size()-1);i++){
// myfile_out << xi_E_minus_current.at(i) << endl;
// }
// }
// myfile_out << "size" << endl;
// myfile_out << xi_E_minus_current.size();
// myfile_out.close();

//////------------------------------------------------------------------------//////

lh_square_current.f_U.assign(para_other.n,1.0);
lh_square_current.q_T.assign(para_other.n,0.0);
lh_square_current.kt_sum_U.assign(para_other.n,0.0);
lh_square_current.f_E.assign(para_other.n,1.0);
lh_square_current.g_E.assign(para_other.n,1.0);
lh_square_current.h_E.assign(para_other.n,1.0);
lh_square_current.k_sum_E.assign(para_other.n,0.0);
lh_square_current.q_E.assign(para_other.n,0.0);
lh_square_current.kt_sum_E.assign(para_other.n,0.0);
lh_square_current.f_I.assign(para_other.n,1.0);
lh_square_current.f_R.assign(para_other.n,1.0);
lh_square_current.f_EnI.assign(para_other.n,1.0);
lh_square_current.f_InR.assign(para_other.n,1.0);

lh_square_current.f_Grid.assign(para_other.n,1.0);
lh_square_current.f_Arc.assign(para_other.n,1.0);

lh_square_current.f_Surv=1.0;


FUNC func_mcmc;

func_mcmc.set_para(para_current, para_other, coordinate, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_R_current, xi_EnI_current, xi_InR_current, t_e_current, t_i_current, t_r_current,epi_final.pop, stb_current, index_current, epi_final.age_gp,infected_source_current , pop_grid, grid_lines);
func_mcmc.initialize_kernel_mat(kernel_mat_current, norm_const_current); // initialize the kernel matrix
func_mcmc.initialize_delta_mat(delta_mat_current); // initialize the kernel matrix
func_mcmc.initialize_lh_square(lh_square_current, kernel_mat_current, delta_mat_current,norm_const_current); //initialize lh_square

myfile_out.open((string(path4)+string("initial_norm_const.txt")).c_str(),ios::app);
for (int i=0;i<=((int)para_other.n-1);i++){
myfile_out<<norm_const_current.at(i)<<endl;
}
myfile_out.close();

double log_lh_current = log_lh_func (lh_square_current, para_other.n); // initialization of log-likelihood value

myfile_out.open((string(path4)+string("initial_lh.txt")).c_str(),ios::app);
myfile_out<<log_lh_current<<endl; //NOTE: this must be defined, otherwise has to re-initialize some components
myfile_out.close();

//---
myfile_out.open((string(path4)+string("initial_f_I.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_I_current.size()-1);i++){
myfile_out<<lh_square_current.f_I.at(xi_I_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_EnI.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_EnI_current.size()-1);i++){
myfile_out<<lh_square_current.f_EnI.at(xi_EnI_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_kt_sum_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.kt_sum_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_g_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.g_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();myfile_out.open((string(path4)+string("initial_h_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.h_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_q_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.q_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.f_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_U.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_U_current.size()-1);i++){
myfile_out<<lh_square_current.f_U.at(xi_U_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_InR.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_InR_current.size()-1);i++){
myfile_out<<lh_square_current.f_InR.at(xi_InR_current.at(i))<<endl;
}
myfile_out.close();


myfile_out.open((string(path4)+string("initial_f_Grid.txt")).c_str(),ios::app);
for (int i=0;i<=(para_other.n-1);i++){
myfile_out<<lh_square_current.f_Grid.at(i)<<endl;
}
myfile_out.close();



myfile_out.open((string(path4)+string("initial_f_Surv.txt")).c_str(),ios::app);
myfile_out<<lh_square_current.f_Surv<<endl;
myfile_out.close();

//----
// myfile_out.open((string(path4)+string("t_r.txt")).c_str(),ios::app);
// for (int i=0;i<=((int)xi_R_current.size()-1);i++){
// myfile_out<<t_r_current.at(xi_R_current.at(i))<<endl;
// }
// myfile_out.close();

myfile_out.open((string(path4)+string("xi_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_current.size()-1);i++){
myfile_out<<xi_E_current.at(i)<<endl;
}
myfile_out.close();


/*--------------------*/
ifstream myfile1_in, myfile2_in, myfile3_in,myfile4_in, myfile5_in,myfile6_in,myfile7_in;
ofstream myfile1_out, myfile2_out, myfile3_out,myfile4_out ,myfile5_out,myfile6_out,myfile7_out; 



mcmc_UPDATE mcmc_update;
mcmc_update.set_para (para_other, coordinate, epi_final.age_gp, epi_final.pop, pop_grid, grid_lines);



myfile1_out.open((string(path4)+string("parameters_current.txt")).c_str(),ios::app);
// myfile1_out << "alpha" << "," << "beta_1" << ","<< "beta_2" << ","<< "beta_3" << ","<< "beta_4" << "," << "beta_5" << "," << "mu_lat" << "," << "var_lat" << "," << "c" << "," << "d" << "," << "k_1" << "," << "k_2" <<","<< "eta" << endl;
myfile1_out << "alpha" << "," << "beta_1" << ","<< "beta_2" << "," << "mu_lat" << "," << "var_lat" << "," << "c" << "," << "d" << "," << "k_1" << "," << "k_2" << "," << "omega" << endl;

myfile2_out.open((string(path4)+string("log_lh_check.txt")).c_str(),ios::app);
myfile3_out.open((string(path4)+string("log_plh_check.txt")).c_str(),ios::app);

myfile4_out.open((string(path4)+string("index_current.txt")).c_str(),ios::app);

myfile_out.open((string(path4)+string("t_i_mcmc.txt")).c_str(),ios::app);

myfile5_out.open((string(path4)+string("t_e_mcmc.txt")).c_str(),ios::app);

myfile6_out.open((string(path4)+string("infected_source_current.txt")).c_str(),ios::app);

myfile7_out.open((string(path4)+string("infected_distance_current.txt")).c_str(),ios::app);

int freq_t_e_up = 100; // frequency to translate an infection time

for (int i=0;i<=(n_iter-1);i++){


mcmc_update.alpha_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);

mcmc_update.omega_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);


mcmc_update.beta_1_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current,infected_source_current, i);
mcmc_update.beta_2_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current,i);

mcmc_update.mu_lat_update(lh_square_current, log_lh_current,xi_U_current, xi_I_current, xi_EnI_current, t_i_current, t_e_current,t_r_current,  index_current, para_current, infected_source_current,i);
mcmc_update.var_lat_update(lh_square_current, log_lh_current, xi_U_current,xi_I_current, xi_EnI_current,t_i_current, t_e_current, t_r_current, index_current, para_current, infected_source_current,i);

mcmc_update.c_update(lh_square_current,log_lh_current, xi_U_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, index_current, para_current, infected_source_current,i);

mcmc_update.k_1_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current,i);




mcmc_update.source_t_e_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current,infected_distance_current, i);






// myfile1_out << para_current.alpha << "," << para_current.beta_age.at(0) << ","<< para_current.beta_age.at(1) << ","<< para_current.beta_age.at(2) << ","<< para_current.beta_age.at(3) <<  ","<< para_current.beta_age.at(4) << "," << para_current.mu_lat << "," << para_current.var_lat << "," << para_current.c << "," << para_current.d << "," << para_current.k_1 << "," <<  para_current.k_2 << "," <<  para_current.eta << endl;

myfile1_out << para_current.alpha << "," << para_current.beta_age.at(0) << ","<< para_current.beta_age.at(1) << "," << para_current.mu_lat << "," << para_current.var_lat << "," << para_current.c << "," << para_current.d << "," << para_current.k_1 << "," <<  para_current.k_2 << "," <<  para_current.omega << endl;


for (int i_index=0; i_index<=((int)index_current.size()-1);i_index++){
myfile4_out << index_current.at(i_index)<<endl;
}




//-------------------------------//

div_t div_iter_lh;
div_iter_lh = div (i,2000);

if (div_iter_lh.rem==0){

lh_SQUARE lh_square_check;

lh_square_check.f_U.assign(para_other.n,1.0);
lh_square_check.q_T.assign(para_other.n,0.0);
lh_square_check.kt_sum_U.assign(para_other.n,0.0);
lh_square_check.f_E.assign(para_other.n,1.0);
lh_square_check.g_E.assign(para_other.n,1.0);
lh_square_check.h_E.assign(para_other.n,1.0);
lh_square_check.k_sum_E.assign(para_other.n,0.0);
lh_square_check.q_E.assign(para_other.n,0.0);
lh_square_check.kt_sum_E.assign(para_other.n,0.0);
lh_square_check.f_I.assign(para_other.n,1.0);
lh_square_check.f_R.assign(para_other.n,1.0);
lh_square_check.f_EnI.assign(para_other.n,1.0);
lh_square_check.f_InR.assign(para_other.n,1.0);

lh_square_check.f_Grid.assign(para_other.n,1.0);
lh_square_check.f_Arc.assign(para_other.n,1.0);

lh_square_check.f_Surv = 1.0;

vector < vector<double> > kernel_check = kernel_mat_current;
vector < vector<double> > delta_mat_check = delta_mat_current;
vector<double> norm_const_check =  norm_const_current;

vector<int> infected_source_check =  infected_source_current;

FUNC func_check;

func_check.set_para(para_current, para_other, coordinate, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_R_current, xi_EnI_current, xi_InR_current, t_e_current, t_i_current, t_r_current, epi_final.pop, stb_current, index_current, epi_final.age_gp, infected_source_check, pop_grid, grid_lines);
func_check.initialize_kernel_mat(kernel_check, norm_const_check); // initialize the kernel matrix
func_check.initialize_delta_mat(delta_mat_check); // initialize the kernel matrix
func_check.initialize_lh_square(lh_square_check, kernel_check, delta_mat_check, norm_const_check); //initialize lh_square


double log_lh_check = log_lh_func (lh_square_check, para_other.n); // initialization of log-likelihood value

myfile2_out <<  log_lh_current << "," <<   log_lh_check << endl;

double log_plh_current=0.0;
double log_plh_check=0.0;
for (int ik=0; ik<=(int)(para_other.n-1);ik++){
// log_plh_current = log_plh_current + log(lh_square_current.f_U.at(ik)*lh_square_current.f_I.at(ik)*lh_square_current.f_R.at(ik)*lh_square_current.f_InR.at(ik));
// log_plh_check = log_plh_check + log(lh_square_check.f_U.at(ik)*lh_square_check.f_I.at(ik)*lh_square_check.f_R.at(ik)*lh_square_check.f_InR.at(ik));
log_plh_current = log_plh_current + lh_square_current.f_Surv;
log_plh_check = log_plh_check + lh_square_check.f_Surv;
}
myfile3_out << log_plh_current << "," <<   log_plh_check << endl;


for (int js=0;js<=(para_other.n-1);js++){
int rem = (js+1)%para_other.n;
if ((rem!=0) | (js==0)) myfile6_out<< infected_source_current.at(js) << ",";
if ((rem==0) & (js!=0)) myfile6_out <<  infected_source_current.at(js) << " " << endl;
}


for (int js=0;js<=(para_other.n-1);js++){
int rem = (js+1)%para_other.n;
if ((rem!=0) | (js==0)) myfile7_out<< infected_distance_current.at(js) << ",";
if ((rem==0) & (js!=0)) myfile7_out <<  infected_distance_current.at(js) << " " << endl;
}



// myfile_out.open((string(path4)+string("real_time_norm_const.txt")).c_str(),ios::out);
// for (int ii=0;ii<=((int)para_other.n-1);ii++){
// myfile_out<<norm_const_current.at(ii)<<endl;
// }
// myfile_out.close();
// 
// myfile_out.open((string(path4)+string("real_time_time_exposure.txt")).c_str(),ios::out);
// for (int ii=0; ii<= (int) (xi_E_current.size()-1);ii++){
// myfile_out  << t_e_current.at(xi_E_current.at(ii)) << endl;
// }
// myfile_out.close();
// 
// myfile_out.open((string(path4)+string("real_time_stb.txt")).c_str(),ios::out);
// for (int i=0; i<=((int)stb.size()-1);i++){
// myfile_out << gp_stb_current.at(i) << "," << stb_current.at(i)  << endl;
// }
// myfile_out.close();


}

//------------------

} // end of MCMC loop

myfile1_out.close();
myfile2_out.close();
myfile3_out.close();
myfile4_out.close();
myfile_out.close();
myfile5_out.close();
myfile6_out.close();
myfile7_out.close();

myfile1_out.open((string(path4)+string("xi_E_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U"  <<","<< "f_E" << endl;
if (xi_E_current.empty()!=1){
for (int i=0; i<=((int)xi_E_current.size()-1);i++){
myfile1_out << xi_E_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_E_current.at(i)) <<","<< lh_square_current.f_U.at(xi_E_current.at(i)) <<","<< lh_square_current.f_E.at(xi_E_current.at(i)) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_E_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_E_minus_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U" <<","<<  "f_E"<< endl;
if (xi_E_minus_current.empty()!=1){
for (int i=0; i<=((int)xi_E_minus_current.size()-1);i++){
myfile1_out << xi_E_minus_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_E_minus_current.at(i)) <<","<< lh_square_current.f_U.at(xi_E_minus_current.at(i)) <<","<< lh_square_current.f_E.at(xi_E_minus_current.at(i))<< endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_E_minus_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_EnI_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U" <<","<<  "f_E"<< endl;
if (xi_EnI_current.empty()!=1){
for (int i=0; i<=((int)xi_EnI_current.size()-1);i++){
myfile1_out << xi_EnI_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_EnI_current.at(i)) <<","<< lh_square_current.f_U.at(xi_EnI_current.at(i)) <<","<< lh_square_current.f_E.at(xi_EnI_current.at(i)) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_EnI_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_U_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "kt_sum_U" <<","<< "f_U" <<","<<  "f_E"<< endl;
if (xi_U_current.empty()!=1){
for (int i=0; i<=((int)xi_U_current.size()-1);i++){
myfile1_out << xi_U_current.at(i) << "," << lh_square_current.kt_sum_U.at(xi_U_current.at(i)) <<","<< lh_square_current.f_U.at(xi_U_current.at(i)) <<","<< lh_square_current.f_E.at(xi_U_current.at(i)) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_U_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("xi_I_final.txt")).c_str(),ios::out);
myfile1_out << "k" << endl;
if (xi_I_current.empty()!=1){
for (int i=0; i<=((int)xi_I_current.size()-1);i++){
myfile1_out << xi_I_current.at(i) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_I_current.size();
myfile1_out.close();



myfile1_out.open((string(path4)+string("xi_InR_final.txt")).c_str(),ios::out);
myfile1_out << "k" << endl;
if (xi_InR_current.empty()!=1){
for (int i=0; i<=((int)xi_InR_current.size()-1);i++){
myfile1_out << xi_InR_current.at(i) << endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << xi_InR_current.size();
myfile1_out.close();

myfile1_out.open((string(path4)+string("index_final.txt")).c_str(),ios::out);
myfile1_out << "k" << "," << "t_e"<< "," <<  "min_t_e"<<endl;
if (index_current.empty()!=1){
for (int i=0; i<=((int)index_current.size()-1);i++){
myfile1_out << index_current.at(i) << "," << t_e_current.at(index_current.at(0)) << ","<< *min_element(t_e_current.begin(),t_e_current.end())<<endl;
}
}
myfile1_out << "size" << endl;
myfile1_out << index_current.size();
myfile1_out.close();

//------------------

return(0);
}


