
//This source file is to perform some tedious and lengthy operations such as read in simulated data //

#include "contact_model_extension_sim_inference_header.h"


void IO_simpara (para_key& para_true_arg, para_aux& para_other_arg){ // this function is called to input/output parameters from simulation


ifstream myfile_in_simpara;
ofstream myfile_out_simpara; 

string line, field;
int line_count=0, field_count=0;

myfile_in_simpara.open((string(path2)+string("parameters_true.txt")).c_str(),ios::in);
line_count =0;

para_true_arg.beta_age.resize(5); 

//while (myfile_in_simdata.good()) {
while (getline(myfile_in_simpara,line)) {

//getline(myfile_in_simdata,line);
stringstream ss(line);
//string field;
field_count=0;

while (getline(ss, field, ',' )) {
stringstream fs (field);

if ((line_count==1) & (field_count==0)) fs >> para_true_arg.alpha;
if ((line_count==1) & (field_count==1)) fs >> para_true_arg.beta_age.at(0);
if ((line_count==1) & (field_count==2)) fs >> para_true_arg.beta_age.at(1);
if ((line_count==1) & (field_count==3)) fs >> para_true_arg.mu_lat;
if ((line_count==1) & (field_count==4)) fs >> para_true_arg.var_lat;
if ((line_count==1) & (field_count==5)) fs >> para_true_arg.c;
if ((line_count==1) & (field_count==6)) fs >> para_true_arg.d;
if ((line_count==1) & (field_count==7)) fs >> para_true_arg.k_1; //note: k_1 and k_2 were named as par_kernel_1 and par_kernel_2 in simulation
if ((line_count==1) & (field_count==8)) fs >> para_true_arg.k_2;
if ((line_count==1) & (field_count==9)) fs >> para_true_arg.omega; 



field_count = field_count + 1;
}

line_count = line_count +1;
}

myfile_in_simpara.close();



// para_true_arg.mu_lat = para_true_arg.a*para_true_arg.b;
// para_true_arg.var_lat = para_true_arg.a*para_true_arg.b*para_true_arg.b;


myfile_out_simpara.open((string(path4)+string("parameters_true.txt")).c_str(),ios::out);
myfile_out_simpara << "alpha" << "," << "beta_1" << "," << "beta_2" << "," << "mu_lat" << "," << "var_lat" << "," << "c" << "," << "d" << "," << "k_1" << "," <<    "k_2"  << ","<< "omega" <<endl;
myfile_out_simpara<< para_true_arg.alpha << "," << para_true_arg.beta_age.at(0) << "," << para_true_arg.beta_age.at(1) << ","  << para_true_arg.mu_lat << "," << para_true_arg.var_lat << "," << para_true_arg.c << "," << para_true_arg.d << "," << para_true_arg.k_1 << "," <<  para_true_arg.k_2 <<  "," << para_true_arg.omega <<endl;
myfile_out_simpara.close();

//--

myfile_in_simpara.open((string(path2)+string("parameters_other.txt")).c_str(),ios::in);
line_count =0;

while (getline(myfile_in_simpara,line)) {

//getline(myfile_in_simdata,line);
stringstream ss(line);
//string field;
field_count=0;

while (getline(ss, field, ',' )) {
stringstream fs (field);
if ((line_count==1) & (field_count==0)) fs >> para_other_arg.n;
if ((line_count==1) & (field_count==1)) fs >> para_other_arg.dimen_x;
if ((line_count==1) & (field_count==2)) fs >> para_other_arg.dimen_y;
if ((line_count==1) & (field_count==3)) fs >> para_other_arg.t_max;
if ((line_count==1) & (field_count==4)) fs >> para_other_arg.unassigned_time;

if ((line_count==1) & (field_count==5)) fs >> para_other_arg.n_row_grid;
if ((line_count==1) & (field_count==6)) fs >> para_other_arg.n_col_grid;

if ((line_count==1) & (field_count==7)) fs >> para_other_arg.grid_size;
if ((line_count==1) & (field_count==8)) fs >> para_other_arg.n_line;

if ((line_count==1) & (field_count==9)) fs >> para_other_arg.x_min;
if ((line_count==1) & (field_count==10)) fs >> para_other_arg.y_min;
if ((line_count==1) & (field_count==11)) fs >> para_other_arg.x_max;
if ((line_count==1) & (field_count==12)) fs >> para_other_arg.y_max;

if ((line_count==1) & (field_count==13)) fs >> para_other_arg.total_pop_of_grid;

if ((line_count==1) & (field_count==14)) fs >> para_other_arg.T;


field_count = field_count + 1;
}

line_count =line_count + 1;

}

myfile_in_simpara.close();



myfile_out_simpara.open((string(path4)+string("parameters_other.txt")).c_str(),ios::out);
myfile_out_simpara << "n"  << "," << "dimen_x" << "," << "dimen_y" << "," << "t_max" << "," << "unassigned_time" << "," << "n_row_grid" << "," << "n_col_grid" <<"," << "grid_size" << "," << "n_line" <<"," << "x_min" <<"," << "y_min" <<"," << "x_max" <<"," << "y_max" << "," << "total_pop_of_grid"  << "," << "T" << endl;
myfile_out_simpara << para_other_arg.n << "," << para_other_arg.dimen_x << "," << para_other_arg.dimen_y << "," << para_other_arg.t_max << "," << para_other_arg.unassigned_time << "," << para_other_arg.n_row_grid<< ","<< para_other_arg.n_col_grid<< "," << para_other_arg.grid_size << "," << para_other_arg.n_line <<"," << para_other_arg.x_min <<"," << para_other_arg.y_min <<"," << para_other_arg.x_max <<"," << para_other_arg.y_max  << "," <<  para_other_arg.total_pop_of_grid << "," << para_other_arg.T << endl;
myfile_out_simpara.close();


}


/*--------------------------------------------*/
//-------------
void IO_simdata (para_key para_true_arg, para_aux para_other_arg, vector < vector<double> >& coordinate_arg, epi_struct& epi_final_arg,vector<int>& index_arg, vector < vector<double> >& pop_grid_arg, grid_lines_struct& grid_lines){ // this function is called to input/output data from simulation

ifstream myfile_in_simdata;
ofstream myfile_out_simdata; 

string line, field;
int line_count=0, field_count=0;


myfile_in_simdata.open((string(path2)+string("coordinate_cases_trans.txt")).c_str(),ios::in);
//string line;
line_count=0;

//coordinate_arg.resize(para_other_arg.n);

while (getline(myfile_in_simdata,line)) {

 //getline(myfile_in_simdata,line);
 stringstream ss(line);
 //string field;
 field_count=0;

 while (getline(ss, field, ',' )) {
 stringstream fs (field);
 fs >> coordinate_arg[line_count][field_count];
 field_count = field_count + 1;
 }

line_count = line_count + 1;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("coordinate_cases_trans.txt")).c_str(),ios::out);
for (int i=0;i<=(para_other_arg.n-1);i++){
myfile_out_simdata << coordinate_arg[i][0] << "," << coordinate_arg[i][1]<< endl;
}
myfile_out_simdata.close();

/*--------------------------------------------*/


epi_final_arg.k.resize(para_other_arg.n);
// epi_final_arg.q.resize(para_other_arg.n);
epi_final_arg.t_e.resize(para_other_arg.n);
epi_final_arg.t_i.resize(para_other_arg.n);
epi_final_arg.t_r.resize(para_other_arg.n);
// epi_final_arg.status.resize(para_other_arg.n);

// epi_final_arg.coor_x.resize(para_other_arg.n);
// epi_final_arg.coor_y.resize(para_other_arg.n);

// epi_final_arg.gp_stb.resize(para_other_arg.n);
epi_final_arg.stb.resize(para_other_arg.n);

epi_final_arg.onset_type.resize(para_other_arg.n);
epi_final_arg.lab_result.resize(para_other_arg.n);

epi_final_arg.age_gp.resize(para_other_arg.n);
epi_final_arg.sex.resize(para_other_arg.n);

epi_final_arg.infected_source.resize(para_other_arg.n);

epi_final_arg.grid_num_row.resize(para_other_arg.n);
epi_final_arg.grid_num_col.resize(para_other_arg.n);

epi_final_arg.pop.resize(para_other_arg.n);

myfile_in_simdata.open((string(path2)+string("ebola.txt")).c_str(),ios::in);
line_count=0;

while (getline(myfile_in_simdata,line)) {

 stringstream ss(line);
 field_count=0;

 while (getline(ss, field, ',' )) {
 stringstream fs (field);
 if ((line_count>=1) & (field_count==0)) fs >> epi_final_arg.k.at(line_count-1);
 // if ((line_count>=1) & (field_count==1)) fs >> epi_final_arg.coor_x.at(line_count-1); 
 // if ((line_count>=1) & (field_count==2)) fs >> epi_final_arg.coor_y.at(line_count-1);
 if ((line_count>=1) & (field_count==1)) fs >> epi_final_arg.t_e.at(line_count-1);
 if ((line_count>=1) & (field_count==2)) fs >> epi_final_arg.t_i.at(line_count-1);
 if ((line_count>=1) & (field_count==3)) fs >> epi_final_arg.t_r.at(line_count-1);
 if ((line_count>=1) & (field_count==4)) fs >> epi_final_arg.stb.at(line_count-1);
 if ((line_count>=1) & (field_count==5)) fs >> epi_final_arg.onset_type.at(line_count-1);
 if ((line_count>=1) & (field_count==6)) fs >> epi_final_arg.lab_result.at(line_count-1);
 if ((line_count>=1) & (field_count==7)) fs >> epi_final_arg.age_gp.at(line_count-1);
 if ((line_count>=1) & (field_count==8)) fs >> epi_final_arg.sex.at(line_count-1);
 if ((line_count>=1) & (field_count==9)) fs >> epi_final_arg.infected_source.at(line_count-1);
 if ((line_count>=1) & (field_count==10)) fs >> epi_final_arg.pop.at(line_count-1);
 if ((line_count>=1) & (field_count==11)) fs >> epi_final_arg.grid_num_row.at(line_count-1);
 if ((line_count>=1) & (field_count==12)) fs >> epi_final_arg.grid_num_col.at(line_count-1);


 field_count = field_count + 1;
 }

line_count = line_count + 1 ;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("ebola_western.txt")).c_str(),ios::app);
myfile_out_simdata << "k"  << ","  << "t_e" << "," << "t_i"<< "," << "t_r"  << "," << "stb"<< "," << "onset_type"<<  "," << "lab_result" << "," << "age" << ","  << "sex" << ","  << "infected_source" << "," << "pop" << "," << "grid_num_row" <<"," << "grid_num_col" <<  endl;
for (int i=0; i<=(para_other_arg.n-1);i++){
myfile_out_simdata << epi_final_arg.k.at(i) << "," << epi_final_arg.t_e.at(i) << "," << epi_final_arg.t_i.at(i)<< "," << epi_final_arg.t_r.at(i)  << "," << epi_final_arg.stb.at(i) <<"," << epi_final_arg.onset_type.at(i) <<"," << epi_final_arg.lab_result.at(i)<<"," << epi_final_arg.age_gp.at(i)<<"," << epi_final_arg.sex.at(i) << "," << epi_final_arg.infected_source.at(i) << ","<< epi_final_arg.pop.at(i) << ","<< epi_final_arg.grid_num_row.at(i) << ","<< epi_final_arg.grid_num_col.at(i) << endl;
}
myfile_out_simdata.close();

/*--------------------------------------------*/

index_arg.reserve(para_other_arg.n);

myfile_in_simdata.open((string(path2)+string("index.txt")).c_str(),ios::in);
line_count=0;


while (getline(myfile_in_simdata,line)) {

 stringstream ss(line);

 while (getline(ss, field)) {
 stringstream fs (field);
 if (line_count>=1) {
 int ind;
 fs >> ind;
 index_arg.push_back(ind); 

 }
 }

line_count = line_count + 1 ;
}
myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("index_assumed.txt")).c_str(),ios::app);
myfile_out_simdata << "k" << endl;
for (int i=0; i<=((int)index_arg.size()-1);i++){
myfile_out_simdata << index_arg.at(i) << endl;
}
myfile_out_simdata.close();



/*--------------------------------------------*/

// grid_data.k_grid.resize(para_other_arg.n_grid);

// grid_data.coor_x.resize(para_other_arg.n_grid);
// grid_data.coor_y.resize(para_other_arg.n_grid);


// grid_data.coor_x_1.resize(para_other_arg.n_grid);
// grid_data.coor_y_1.resize(para_other_arg.n_grid);

// grid_data.coor_x_2.resize(para_other_arg.n_grid);
// grid_data.coor_y_2.resize(para_other_arg.n_grid);

// grid_data.coor_x_3.resize(para_other_arg.n_grid);
// grid_data.coor_y_3.resize(para_other_arg.n_grid);

// grid_data.coor_x_4.resize(para_other_arg.n_grid);
// grid_data.coor_y_4.resize(para_other_arg.n_grid);

// grid_data.pop.resize(para_other_arg.n_grid);


// myfile_in_simdata.open((string(path2)+string("grid_data.txt")).c_str(),ios::in);
// line_count=0;

// while (getline(myfile_in_simdata,line)) {

//  stringstream ss(line);
//  field_count=0;

//  while (getline(ss, field, ',' )) {
//  stringstream fs (field);

//  if ((line_count>=1) & (field_count==0)) fs >> grid_data.k_grid.at(line_count-1);

//  if ((line_count>=1) & (field_count==1)) fs >> grid_data.coor_x.at(line_count-1);
//  if ((line_count>=1) & (field_count==2)) fs >> grid_data.coor_y.at(line_count-1);


//  if ((line_count>=1) & (field_count==3)) fs >> grid_data.coor_x_1.at(line_count-1);
//  if ((line_count>=1) & (field_count==4)) fs >> grid_data.coor_y_1.at(line_count-1);

//  if ((line_count>=1) & (field_count==5)) fs >> grid_data.coor_x_2.at(line_count-1);
//  if ((line_count>=1) & (field_count==6)) fs >> grid_data.coor_y_2.at(line_count-1);

//  if ((line_count>=1) & (field_count==7)) fs >> grid_data.coor_x_3.at(line_count-1);
//  if ((line_count>=1) & (field_count==8)) fs >> grid_data.coor_y_3.at(line_count-1);

//  if ((line_count>=1) & (field_count==9)) fs >> grid_data.coor_x_4.at(line_count-1);
//  if ((line_count>=1) & (field_count==10)) fs >> grid_data.coor_y_4.at(line_count-1);

//  if ((line_count>=1) & (field_count==11)) fs >> grid_data.pop.at(line_count-1);

//  field_count = field_count + 1;
//  }

// line_count = line_count + 1 ;
// }

// myfile_in_simdata.close();


// myfile_out_simdata.open((string(path4)+string("grid_data.txt")).c_str(),ios::app);
// myfile_out_simdata << "k_grid" << "," << "coor_x"  << ","  << "coor_y"<< "," << "coor_x_1"  << ","  << "coor_y_1" << ","<< "coor_x_2"  << ","  << "coor_y_2" << ","<< "coor_x_3"  << ","  << "coor_y_3" << "," << "coor_x_4"  << ","  << "coor_y_4" << "," <<  "pop" << endl;
// for (int i=0; i<=(para_other_arg.n_grid-1);i++){
// myfile_out_simdata << grid_data.k_grid[i] << "," << grid_data.coor_x[i] << "," << grid_data.coor_y[i]<< "," << grid_data.coor_x_1[i] << "," << grid_data.coor_y_1[i]<< "," << grid_data.coor_x_2[i] << "," << grid_data.coor_y_2[i]<< "," <<grid_data.coor_x_3[i] << "," << grid_data.coor_y_3[i]<< "," << grid_data.coor_x_4[i] << "," << grid_data.coor_y_4[i]<< "," << grid_data.pop[i] << endl;
// }
// myfile_out_simdata.close();

/*--------------------------------------------*/


grid_lines.k_line.resize(para_other_arg.n_line);

grid_lines.orient_line.resize(para_other_arg.n_line);


grid_lines.coor_x_1.resize(para_other_arg.n_line);
grid_lines.coor_y_1.resize(para_other_arg.n_line);

grid_lines.coor_x_2.resize(para_other_arg.n_line);
grid_lines.coor_y_2.resize(para_other_arg.n_line);



myfile_in_simdata.open((string(path2)+string("grid_lines.txt")).c_str(),ios::in);
line_count=0;

while (getline(myfile_in_simdata,line)) {

 stringstream ss(line);
 field_count=0;

 while (getline(ss, field, ',' )) {
 stringstream fs (field);

 if ((line_count>=1) & (field_count==0)) fs >> grid_lines.k_line.at(line_count-1);

 if ((line_count>=1) & (field_count==1)) fs >> grid_lines.orient_line.at(line_count-1);


 if ((line_count>=1) & (field_count==2)) fs >> grid_lines.coor_x_1.at(line_count-1);
 if ((line_count>=1) & (field_count==3)) fs >> grid_lines.coor_y_1.at(line_count-1);

 if ((line_count>=1) & (field_count==4)) fs >> grid_lines.coor_x_2.at(line_count-1);
 if ((line_count>=1) & (field_count==5)) fs >> grid_lines.coor_y_2.at(line_count-1);


 field_count = field_count + 1;
 }

line_count = line_count + 1 ;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("grid_lines.txt")).c_str(),ios::app);
myfile_out_simdata << "k_line" << "," << "orient_line" << ","  << "coor_x_1"  << ","  << "coor_y_1" << ","<< "coor_x_2"  << ","  << "coor_y_2" << endl;
for (int i=0; i<=(para_other_arg.n_line-1);i++){
myfile_out_simdata << grid_lines.k_line[i] << ","  <<  grid_lines.orient_line[i] <<  "," << grid_lines.coor_x_1[i] << "," << grid_lines.coor_y_1[i]<< "," << grid_lines.coor_x_2[i] << "," << grid_lines.coor_y_2[i] << endl;
}
myfile_out_simdata.close();


/*--------------------------------------------*/


myfile_in_simdata.open((string(path2)+string("pop_grid.txt")).c_str(),ios::in);
//string line;
line_count=0;

//coordinate_arg.resize(para_other_arg.n);

while (getline(myfile_in_simdata,line)) {

 //getline(myfile_in_simdata,line);
 stringstream ss(line);
 //string field;
 field_count=0;

 while (getline(ss, field, ',' )) {
 stringstream fs (field);
 fs >> pop_grid_arg[line_count][field_count];
 field_count = field_count + 1;
 }

line_count = line_count + 1;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("pop_grid.txt")).c_str(),ios::out);
for (int i=0;i<=(para_other_arg.n_row_grid-1);i++) {
 for (int j=0;j<=(para_other_arg.n_col_grid-1);j++) {
 	if (j<(para_other_arg.n_col_grid-1)) myfile_out_simdata << pop_grid_arg[i][j] << ",";
	if (j==(para_other_arg.n_col_grid-1)) myfile_out_simdata << pop_grid_arg[i][j] << " " << endl;
}
}
myfile_out_simdata.close();

/*--------------------------------------------*/


}

//---------------------------------------------//

// void lh_plot (para_key para_true_arg, para_aux para_other_arg, vector< vector<double> > coordinate_arg, 
// vector<int> xi_U_arg, vector<int> xi_E_arg, vector<int> xi_E_minus_arg, vector<int> xi_I_arg, vector<int> xi_R_arg, vector<int>xi_EnI_arg,vector<int> xi_InR_arg,  vector<double> t_e_arg, vector<double> t_i_arg, vector<double> t_r_arg, vector<int>index_arg){

// ofstream myfile_out_lh_plot; 

// para_key para_reset = para_true_arg; 
// lh_SQUARE lh_square_internal;
// vector< vector<double> > kernel_mat_internal(para_other_arg.n, vector<double>(para_other_arg.n)), delta_mat_internal(para_other_arg.n, vector<double>(para_other_arg.n));


// FUNC func_internal;

// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg, xi_EnI_arg, xi_InR_arg, t_e_arg, t_i_arg, t_r_arg, index_arg);
// func_internal.initialize_kernel_mat(kernel_mat_internal);
// func_internal.initialize_delta_mat(delta_mat_internal);


// for (int i=1; i<=8;i++){


	// if (i==1){

	// double test_alpha_initial=0.0005;
	// double increment_alpha=0.00001;
	// int num_test_alpha = (int) (0.0015/increment_alpha);
	// vector<double> log_lh_alpha_test(num_test_alpha), test_alpha_vector(num_test_alpha);

        // for (int j=0;j<=(num_test_alpha-1);j++){


	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);	
	
	// test_alpha_vector.at(j) = test_alpha_initial + j*increment_alpha;
	// para_reset.alpha = test_alpha_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg, t_e_arg, t_i_arg, t_r_arg, index_arg);
	// //func_internal.initialize_kernel_mat_internal(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_alpha_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }

	// myfile_out_lh_plot.open((string(path4)+string("log_value_alpha_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_alpha-1);i++) {
	// myfile_out_lh_plot << test_alpha_vector.at(i) << "," << log_lh_alpha_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();

        // }

	// if (i==2){

	// double test_beta_initial=0.8;
	// double increment_beta=0.002;
	// int num_test_beta = (int) (0.2/increment_beta);
	// vector<double> log_lh_beta_test(num_test_beta), test_beta_vector(num_test_beta);
	

        // for (int j=0;j<=(num_test_beta-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_beta_vector.at(j) = test_beta_initial + j*increment_beta;
	// para_reset.beta = test_beta_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// //func_internal.initialize_kernel_mat_internal(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_beta_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);

	// myfile_out_lh_plot.open((string(path4)+string("log_value_beta_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_beta-1);i++) {
	// myfile_out_lh_plot << test_beta_vector.at(i) << "," << log_lh_beta_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();

	// }
        // }

	// if (i==3){

	// double test_a_initial=1.5;
	// double increment_a=0.01;
	// int num_test_a = (int) (1.0/increment_a);
	// vector<double> log_lh_a_test(num_test_a), test_a_vector(num_test_a);
	

        // for (int j=0;j<=(num_test_a-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_a_vector.at(j) = test_a_initial + j*increment_a;
	// para_reset.a = test_a_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// //func_internal.initialize_kernel_mat_internal(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_a_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }
	
	// myfile_out_lh_plot.open((string(path4)+string("log_value_a_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_a-1);i++) {
	// myfile_out_lh_plot << test_a_vector.at(i) << "," << log_lh_a_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();
	
        // }

	// if (i==4){

	// double test_b_initial=0.8;
	// double increment_b=0.005;
	// int num_test_b = (int) (0.4/increment_b);
	// vector<double> log_lh_b_test(num_test_b), test_b_vector(num_test_b);
	
        // for (int j=0;j<=(num_test_b-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_b_vector.at(j) = test_b_initial + j*increment_b;
	// para_reset.b = test_b_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// //func_internal.initialize_kernel_mat_internal(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_b_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }

	
	// myfile_out_lh_plot.open((string(path4)+string("log_value_b_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_b-1);i++) {
	// myfile_out_lh_plot << test_b_vector.at(i) << "," << log_lh_b_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();


        // }

	// if (i==5){

	// double test_c_initial=2.8;
	// double increment_c=0.001;
	// int num_test_c = (int) (0.4/increment_c);
	// vector<double> log_lh_c_test(num_test_c), test_c_vector(num_test_c);
	
        // for (int j=0;j<=(num_test_c-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_c_vector.at(j) = test_c_initial + j*increment_c;
	// para_reset.c = test_c_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// //func_internal.initialize_kernel_mat_internal(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_c_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }

	// myfile_out_lh_plot.open((string(path4)+string("log_value_c_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_c-1);i++) {
	// myfile_out_lh_plot << test_c_vector.at(i) << "," << log_lh_c_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();

        // }

	// if (i==6){

	// double test_d_initial=2.7;
	// double increment_d=0.001;
	// int num_test_d = (int) (0.5/increment_d);
	// vector<double> log_lh_d_test(num_test_d), test_d_vector(num_test_d);
	
        // for (int j=0;j<=(num_test_d-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_d_vector.at(j) = test_d_initial + j*increment_d;
	// para_reset.d = test_d_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// //func_internal.update_initialize_mat_internal(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_d_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }

	// myfile_out_lh_plot.open((string(path4)+string("log_value_d_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_d-1);i++) {
	// myfile_out_lh_plot << test_d_vector.at(i) << "," << log_lh_d_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();

        // }

	// if (i==7){

	// double test_k_1_initial=32.0;
	// double increment_k_1=0.005;
	// int num_test_k_1 = (int) (4.0/increment_k_1);
	// vector<double> log_lh_k_1_test(num_test_k_1), test_k_1_vector(num_test_k_1);
	
        // for (int j=0;j<=(num_test_k_1-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_k_1_vector.at(j) = test_k_1_initial + j*increment_k_1;
	// para_reset.k_1 = test_k_1_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// func_internal.initialize_kernel_mat(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_k_1_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }

	// myfile_out_lh_plot.open((string(path4)+string("log_value_k_1_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_k_1-1);i++) {
	// myfile_out_lh_plot << test_k_1_vector.at(i) << "," << log_lh_k_1_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();

        // }

	// if (i==8){

	// double test_k_2_initial=13.0;
	// double increment_k_2=0.005;
	// int num_test_k_2 = (int) (3.0/increment_k_2);
	// vector<double> log_lh_k_2_test(num_test_k_2), test_k_2_vector(num_test_k_2);
	
        // for (int j=0;j<=(num_test_k_2-1);j++){

	// para_reset = para_true_arg;

	// lh_square_internal.f_U.assign(para_other_arg.n,1.0);
	// lh_square_internal.q_T.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_U.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.g_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.h_E.assign(para_other_arg.n,1.0);
	// lh_square_internal.k_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.q_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.kt_sum_E.assign(para_other_arg.n,0.0);
	// lh_square_internal.f_I.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_R.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_EnI.assign(para_other_arg.n,1.0);
	// lh_square_internal.f_InR.assign(para_other_arg.n,1.0);		
	
	// test_k_2_vector.at(j) = test_k_2_initial + j*increment_k_2;
	// para_reset.k_2 = test_k_2_vector.at(j) ;
	// func_internal.set_para(para_reset, para_other_arg, coordinate_arg, xi_U_arg, xi_E_arg, xi_E_minus_arg, xi_I_arg, xi_R_arg,xi_EnI_arg, xi_InR_arg,  t_e_arg, t_i_arg, t_r_arg, index_arg);
	// func_internal.initialize_kernel_mat(kernel_mat_internal);
	// func_internal.initialize_lh_square(lh_square_internal, kernel_mat_internal,delta_mat_internal);
	// log_lh_k_2_test.at(j) = log_lh_func (lh_square_internal,para_other_arg.n);
	// }

	// myfile_out_lh_plot.open((string(path4)+string("log_value_k_2_test.txt")).c_str(),ios::out);
	// for (int i=0;i<= (int) (num_test_k_2-1);i++) {
	// myfile_out_lh_plot << test_k_2_vector.at(i) << "," << log_lh_k_2_test.at(i) <<endl;
	// }
	// myfile_out_lh_plot.close();

        // }

// }

// }


