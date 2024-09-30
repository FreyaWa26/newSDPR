#include "mcmc.h"
#include <cmath>
#include <random>
#include "gsl/gsl_randist.h"
#include <gsl/gsl_matrix.h>
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <iostream>
#include <chrono>
#include <functional>


using std::log; using std::exp; using std::sqrt;
using std::cout; using std::endl; using std::ref;
using std::cerr; 

#define square(x) ((x)*(x))

void init_state(Dat *dat, MCMC_state *state) {
	//state->r = gsl_rng_alloc(gsl_rng_default);
    state->n_cluster = 20;
	//gsl_vector *V = gsl_vector_alloc(20);
	state->rho = dat->rho_1;
	state->rho_a = dat->rho_3;
	state->rho_b = dat->rho_2;
	state->det_sigma = 1-square(dat->rho_1)-square(dat->rho_2)-square(dat->rho_3)+2*dat->rho_1*dat->rho_2*dat->rho_3;
	
	// Compute the inverse matrix elements
	state->inv[0][0] = (1 - square(state->rho_b)) / state->det_sigma;
	state->inv[0][1] = (state->rho_b * state->rho_a - state->rho) / state->det_sigma;
	state->inv[0][2] = (state->rho * state->rho_b - state->rho_a) / state->det_sigma;

	state->inv[1][0] = state->inv[0][1]; // The matrix is symmetric
	state->inv[1][1] = (1 - square(state->rho_a) )/ state->det_sigma;
	state->inv[1][2] = (state->rho * state->rho_a - state->rho_b) / state->det_sigma;

	state->inv[2][0] = state->inv[0][2]; // The matrix is symmetric
	state->inv[2][1] = state->inv[1][2]; // The matrix is symmetric
	state->inv[2][2] = (1 - square(state->rho)) / state->det_sigma;

    state->beta1 = (double *) calloc(dat->n_snp, sizeof(double));
    state->beta2 = (double *) calloc(dat->n_snp, sizeof(double));
	state->beta3 = (double *) calloc(dat->n_snp, sizeof(double));
	state->b1 = (double *) calloc(dat->n_snp, sizeof(double));
    state->b2 = (double *) calloc(dat->n_snp, sizeof(double));
	state->b3 = (double *) calloc(dat->n_snp, sizeof(double));
	for (size_t i=0; i<dat->n_snp; i++) {
	state->beta1[i] = 0;
	state->beta2[i] = 0;
	state->beta3[i] = 0;
	state->b1[i] = 0;
	state->b2[i] = 0;
	state->b3[i] = 0;
    }

    state->residual = (double *) malloc(sizeof(double)*dat->n_ind);
    state->G = (double *) calloc(dat->n_ind, sizeof(double));//geno*beta sum
	state->det_sigma = 1-square(state->rho)-square(state->rho_a)-square(state->rho_b)+2*state->rho*state->rho_a*state->rho_b;
    state->n_snp = dat->n_snp;

    dat->y = (double *) malloc(sizeof(double)*dat->n_ind);
    for (size_t i=0; i<dat->n_ind; i++) {
	dat->y[i] = dat->pheno[i];
    }

    for (size_t i=0; i<dat->n_ind; i++) {
	state->residual[i] = dat->pheno[i];
    }

    state->assgn = (int *) malloc(sizeof(int)*dat->n_snp);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(1, state->n_cluster);
    
    //std::ifstream infile("../effect_size/sim_1.txt");
    for (size_t i=0; i<dat->n_snp; i++) {
	state->assgn[i] = dist(gen);
	//infile >> state->assgn[i];
    }
	state->p_0 = 0.95;
	state->p = (double *) calloc(state->n_cluster, sizeof(double));
    state->suffstats = (int *) calloc(state->n_cluster, sizeof(int));
    for (size_t i=0; i<state->n_cluster; i++) {
	state->p[i] = 1.0 / state->n_cluster;
	state->suffstats[state->assgn[i]]++;
    }
    state->sumsq = (double *) calloc(state->n_cluster, sizeof(double));
    
    state->cluster_var = (double *) malloc(sizeof(double)*state->n_cluster);
    for (size_t i=1; i<state->n_cluster; i++) {
		std::gamma_distribution<> rgamma(state->suffstats[i]/2.0+state->a0k, \
			1.0/state->b0k);
		state->cluster_var[i] = 1.0/rgamma(gen);
	}
    
    state->tau = 1.0;
    state->log_p = (double *) malloc(sizeof(double)*state->n_cluster);
    state->V = (double *) malloc(sizeof(double)*state->n_cluster);
    
}

void destroy_state(MCMC_state *state) {
    free(state->beta1);
    free(state->beta2);
	free(state->beta3);
    free(state->residual);
    free(state->G);
    free(state->assgn);
    free(state->suffstats);
    free(state->sumsq);
    free(state->cluster_var);
    free(state->p);
    free(state->log_p);
    free(state->V);
    
}

void MCMC_state::calc_b(size_t j, const Dat *dat, const ldmat_data &ldmat_dat) {
	size_t start_i = dat->boundary[j].first;
    size_t end_i = dat->boundary[j].second;
	size_t len_blk = end_i - start_i;
	double *res = residual;
	//cout<<" b"<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
	for (size_t i = 0; i < len_blk; i++) {
        double sum1 = 0.0, sum2 = 0.0;
        for (size_t k = 0; k < dat->n_ind; k++) {
            sum1 += dat->geno1[start_i + i][k] * res[k];
            sum2 += dat->geno2[start_i + i][k] * res[k];
        }
		
        b1[start_i + i] = eta * sum1;
        b2[start_i + i] = eta * sum2;
        b1[start_i + i] += dat->geno1_sq[start_i + i] * beta1[start_i + i] + 
                           dat->geno12_prod[start_i + i] * beta2[start_i + i];
        b2[start_i + i] += dat->geno2_sq[start_i + i] * beta2[start_i + i] + 
                           dat->geno12_prod[start_i + i] * beta1[start_i + i];
		
        b1[start_i + i] *= eta / sigmae2;
        b2[start_i + i] *= eta / sigmae2;
    }
	//cout<<" end"<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
    const double* beta_j = &dat->beta_mrg[start_i];
    gsl_vector_const_view diag = gsl_matrix_const_diagonal(ldmat_dat.B[j]);
	
    gsl_vector* temp = gsl_vector_alloc(len_blk);
	// diag(B) \times beta
	for (int i = 0; i < len_blk; i++) {
		gsl_vector_set(temp, i, beta_j[i] * gsl_vector_get(&diag.vector, i));
	}
	//std::cout << "diagB times beta " << std::endl;

	// eta^2 * (diag(B) \times beta) - eta^2 * B beta
	gsl_vector* beta_view = gsl_vector_alloc(len_blk);
	for (int i = 0; i < len_blk; i++) {
		gsl_vector_set(beta_view, i, beta_j[i]);
	}
	
	gsl_blas_dsymv(CblasUpper, -eta*eta, ldmat_dat.B[j], beta_view, eta*eta, temp);

	// eta^2 * (diag(B) \times beta) - eta^2 * B beta + eta * A^T beta_mrg
	gsl_blas_daxpy(eta, ldmat_dat.calc_b_tmp[j], temp);
	//std::cout << "start copy back " << std::endl;

	// Copy the result back to b3j
	for (int i = 0; i < len_blk; i++) {
		b3[start_i+i] = gsl_vector_get(temp, i);
	}
	
	// Free allocated memory
	gsl_vector_free(temp);
	gsl_vector_free(beta_view);
	//std::cout << "cal b " << j <<" 1 "<< b1[0]<<"2 "<<b2[0]<<"3 "<<b3[0]<< std::endl;
}

void MCMC_state::sample_sigma2() {
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    for (size_t i=1; i<n_cluster; i++) {
	double a = suffstats[i] *1.5 + a0k;
	double b = 1.0 / (sumsq[i] / 2.0 + b0k);
	
	cluster_var[i] = 1.0/gsl_ran_gamma(r, a, b);
	//std::cout << "var"<<cluster_var[i]<< std::endl;
	if (isinf(cluster_var[i])) {
	    cluster_var[i] = 1e5;
	    std::cerr << "Cluster variance is infintie." << std::endl;
	}
	else if (cluster_var[i] == 0) {
	    cluster_var[i] = 1e-10;
	    std::cerr << "Cluster variance is zero." << std::endl;
	}
    }
	gsl_rng_free(r);
}

void MCMC_state::sample_tau(const Dat *dat){
	double shape = dat->n_ind/2.0+a0k;
	double sum_of_squares = 0.0;

	for (size_t i = 0; i < dat->n_ind; i++) {
        sum_of_squares += square(residual[i]);
    }
	
    double rate = sum_of_squares / 2.0 + b0k;
	//cout << shape<<"rate" <<rate<< endl;
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	tau = gsl_ran_gamma(r, shape, 1/rate);
	//cout << "tau done" <<tau<< endl;
}

void MCMC_state::update_suffstats() {
	
    // Reset suffstats and sumsq arrays
    for (int k = 0; k < n_cluster; k++) {
        suffstats[k] = 0;
        sumsq[k] = 0.0;
	}

    for (size_t i = 0; i < n_snp; i++) {
        suffstats[assgn[i]]++;
		sumsq[assgn[i]] += ( (1-square(rho_b))*square(beta1[i]) + (1-square(rho_a))*square(beta2[i]) +(1-square(rho))*square(beta3[i]) \
		- 2 * (rho-rho_a*rho_b) * beta1[i] * beta2[i]\
		- 2 * (rho_a-rho*rho_b) * beta1[i] * beta3[i]\
		- 2 * (rho_b-rho_a*rho) * beta2[i] * beta3[i] ) / det_sigma ;
    }
}

void MCMC_state::sample_V() {
	V[0] = 0;
    vector<double> a(n_cluster);
    a[n_cluster-1] = suffstats[n_cluster-1];
    for (int i=n_cluster-2; i>0; i--) {
	a[i] = suffstats[i] + a[i+1];
    }
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    for (size_t i=1; i<n_cluster; i++) {
	V[i] = gsl_ran_beta(r, 1 + suffstats[i], alpha + a[i]);
    }
    V[n_cluster-1] = 1;
	gsl_rng_free(r);
}

void MCMC_state::update_p() {
    vector<double> cumprod(n_cluster);
    cumprod[0] = 1.0;
    cumprod[1] = 1.0 - V[1];
	
	for (size_t i=2; i<n_cluster; i++) {
	cumprod[i] = cumprod[i-1] * (1 - V[i]);
	if (V[i] == 1) {
	    std::fill(cumprod.begin()+i+1, cumprod.end(), 0.0);
	    break;
	}
    }
    p[0] = 0; 
	
    for (size_t i=1; i<n_cluster; i++) {
		p[i] = cumprod[i-1] * V[i];
    }
	
	double sum = std::accumulate(p, p + n_cluster - 1, 0.0);
    //double sum = std::accumulate(p.begin(), p.end()-1, 0.0);
    if (1 - sum > 0) {
	p[n_cluster-1] = 1 - sum;
    }
    else {
	p[n_cluster-1] = 0;
    }

    for (size_t i=0; i<n_cluster; i++) {
	log_p[i] = logf(p[i] + 1e-40); 
    }
}
void MCMC_state::update_p0(const Dat *dat){
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	p_0 = gsl_ran_beta(r,suffstats[0],dat->n_ind);
	gsl_rng_free(r);
}
void MCMC_state::update_residual(const Dat *dat) {
	for (size_t i = 0; i < dat->n_ind; i++) {
    double G_i = 0.0;
    for (size_t j = 0; j < n_snp; j++) {
        G_i += beta1[j] * dat->geno1[j][i] + beta2[j] * dat->geno2[j][i];
    }
	G[i] =G_i;
    residual[i] = dat->pheno[i]-G_i;//res = Y-aW-X1b1-X2b2
	}
}


// Function to extract a subset of the double** matrix into a GSL matrix
void extract_submatrix(double** large_matrix, gsl_matrix* sub_matrix, size_t start_i, size_t end_i, size_t num_cols) {
    for (size_t i = start_i; i < end_i; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            gsl_matrix_set(sub_matrix, i - start_i, j, large_matrix[i][j]);  // Copy row to submatrix
        }
    }
}

void print_matrix(const gsl_matrix* m) {
    if (!m) {
        std::cout << "Error: Null matrix pointer." << std::endl;
        return;
    }

    size_t rows = m->size1;
    size_t cols = m->size2;

    std::cout << "Matrix (" << rows << "x" << cols << "):" << std::endl;
	if(rows>100){rows = 100;}
	if(cols>100){cols = 100;}
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            cout << gsl_matrix_get(m, i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
/*void adjust_covariance_matrix(gsl_matrix* cov_mat) {
    int n = cov_mat->size1;
    gsl_vector* eval = gsl_vector_alloc(n);
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc(n);

    while (true) {
        // Compute eigenvalues
        gsl_matrix* cov_mat_copy = gsl_matrix_alloc(n, n);
        gsl_matrix_memcpy(cov_mat_copy, cov_mat);
        gsl_eigen_symm(cov_mat_copy, eval, w);
        gsl_matrix_free(cov_mat_copy);

        // Find minimum eigenvalue
        double min_eigenvalue = gsl_vector_min(eval);

        if (min_eigenvalue >= 0) {
            break;  // Exit loop if all eigenvalues are non-negative
        }

        // Adjust the covariance matrix
		cout<<"adjust covariance matrix"<<endl;
        double adjustment = std::max(1.0, -min_eigenvalue);
        for (int i = 0; i < n; ++i) {
            double* element = gsl_matrix_ptr(cov_mat, i, i);
            *element += adjustment;
        }
    }

    gsl_vector_free(eval);
    gsl_eigen_symm_free(w);
}*/
void MCMC_state::sample_beta(size_t j, const Dat *dat, const ldmat_data &ldmat_dat,unsigned sz){
	
	size_t start_i = dat->boundary[j].first;
    size_t end_i = dat->boundary[j].second;

    std::vector <size_t>causal_list;
	
    for (size_t i=start_i; i<end_i; i++) {
		if (assgn[i] != 0) {
			causal_list.push_back(i);
			
		}
    }
	
	if (causal_list.size() == 0) {
		return;
    }
	//std::cout << "causal list"<< std::endl;
	size_t mat_size = 3 * causal_list.size();
	size_t sub_size = causal_list.size();
	gsl_vector *A_vec = gsl_vector_alloc(mat_size);
    //gsl_vector *A_vec2 = gsl_vector_alloc(mat_size);
	//A = res + X1 * beta1 +X2 * beta2
	
	// Initialize A with residual
	gsl_vector_const_view beta1_gsl = gsl_vector_const_view_array(beta1, n_snp);
    gsl_vector_const_view beta1_here = gsl_vector_const_subvector(&beta1_gsl.vector, start_i, end_i-start_i);    
    gsl_vector_const_view beta2_gsl = gsl_vector_const_view_array(beta2, n_snp);
    gsl_vector_const_view beta2_here = gsl_vector_const_subvector(&beta2_gsl.vector, start_i, end_i-start_i);    
    
    // A = res + X1 * beta1 + X2 * beta2
    gsl_vector *A_res = gsl_vector_alloc(dat->n_ind);
    gsl_matrix* X1 = gsl_matrix_alloc(end_i - start_i, dat->n_ind);  // Submatrix for X1
    gsl_matrix* X2 = gsl_matrix_alloc(end_i - start_i, dat->n_ind);  // Submatrix for X2

    // Initialize A with residual
    //double max_a = 0.0;
	for (size_t i = 0; i < dat->n_ind; i++) {
		gsl_vector_set(A_res, i, residual[i]);
		//if(max_a<residual[i]){max_a = residual[i];}
	}
	//cout<<"initial max A "<<max_a<<endl;
    extract_submatrix(dat->geno1, X1, start_i, end_i, dat->n_ind);
    extract_submatrix(dat->geno2, X2, start_i, end_i, dat->n_ind);  
	
    // Perform matrix-vector multiplications
    gsl_blas_dgemv(CblasTrans, 1.0, X1, &beta1_here.vector, 1.0, A_res);  // A += X1^T * beta1
	//std::cout << "A after X1: " << gsl_vector_get(A, 0) << std::endl;
	
    gsl_blas_dgemv(CblasTrans, 1.0, X2, &beta2_here.vector, 1.0, A_res);  // A += X2^T * beta2
	//std::cout << "A after X2: " << gsl_vector_get(A, 0) << std::endl;
	
    double C = square(eta)*sz; 

    gsl_matrix* ptr = gsl_matrix_alloc(mat_size, mat_size);
    
    if (!ptr) {
	std::cerr << "Malloc failed for block " 
	    " may due to not enough memory." << endl;
    }
	gsl_matrix_set_zero(ptr); 

	for (size_t i = 0; i < sub_size; i++) {
		gsl_vector_const_view geno1_i = gsl_vector_const_view_array(dat->geno1[causal_list[i]], dat->n_ind);
		gsl_vector_const_view geno2_i = gsl_vector_const_view_array(dat->geno2[causal_list[i]], dat->n_ind);
		
		for(size_t k = 0; k < sub_size; k++){
			gsl_vector_const_view geno1_k = gsl_vector_const_view_array(dat->geno1[causal_list[k]], dat->n_ind);
			gsl_vector_const_view geno2_k = gsl_vector_const_view_array(dat->geno2[causal_list[k]], dat->n_ind);
			double result1, result2, result12, result21;
			
			gsl_blas_ddot(&geno1_i.vector, &geno1_k.vector, &result1);
			gsl_blas_ddot(&geno2_i.vector, &geno2_k.vector, &result2);
			gsl_blas_ddot(&geno1_i.vector, &geno2_k.vector, &result12);
			gsl_blas_ddot(&geno2_i.vector, &geno1_k.vector, &result21);
			// Upper-left 2x2 block
            gsl_matrix_set(ptr, i, k, tau * result1);
            gsl_matrix_set(ptr, sub_size + i, sub_size + k, tau * result2);
            gsl_matrix_set(ptr, i, sub_size + k, tau * result12);
            gsl_matrix_set(ptr, sub_size + i, k, tau * result21);

            // Lower-right block
            gsl_matrix_set(ptr, 2 * sub_size + i, 2 * sub_size + k,
                           C * ldmat_dat.B[j]->data[ldmat_dat.B[j]->tda * (causal_list[i] - start_i) + (causal_list[k] - start_i)]);
		}

		// A_vec calculations
		double dot_product1 = 0.0, dot_product2 = 0.0;
		for (size_t k = 0; k < dat->n_ind; ++k) {
			dot_product1 += tau * gsl_vector_get(A_res, k) * dat->geno1[causal_list[i]][k];
			dot_product2 += tau * gsl_vector_get(A_res, k) * dat->geno2[causal_list[i]][k];
		}
		gsl_vector_set(A_vec, i, dot_product1);
		gsl_vector_set(A_vec, sub_size+i, dot_product2);
		gsl_vector_set(A_vec, 2*sub_size+i, C*gsl_vector_get(ldmat_dat.calc_b_tmp[j], causal_list[i]-start_i));
	}
	
	//compute the \Sigma_0^-1
	// Step 1: Create diagonal matrix inverse 
    gsl_matrix *diag = gsl_matrix_alloc(sub_size, sub_size);
    gsl_matrix_set_zero(diag);
    for (int i = 0; i < sub_size; i++) {
        gsl_matrix_set(diag, i, i,1.0/ cluster_var[assgn[causal_list[i]]]);
    }
	
    
    // Step 3: Create the large variance matrix
    gsl_matrix *var_mat = gsl_matrix_calloc(mat_size, mat_size);
    
    // Fill var_mat
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            // Define the submatrix view in var_mat where data will be copied and scaled
            gsl_matrix_view submat = gsl_matrix_submatrix(var_mat, i * sub_size, k * sub_size, sub_size, sub_size);

            // Copy the diagonal matrix to the submatrix
            gsl_matrix_memcpy(&submat.matrix, diag);

            // Scale the submatrix by the corresponding value from cov_matrix
           
            gsl_matrix_scale(&submat.matrix, inv[i][k]);
        }
	
    }
	
	// Step 4: Add matrices
	gsl_matrix* result = gsl_matrix_alloc(mat_size, mat_size);
	// Copy var_mat to result
    gsl_matrix_memcpy(result, var_mat);

    // Add ptr_mat to result
    gsl_matrix_add(result, ptr);
	
    gsl_vector *beta_c = gsl_vector_alloc(mat_size);
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    for (size_t i=0; i<mat_size; i++) {
		gsl_vector_set(beta_c, i, gsl_ran_ugaussian(r));
    }
	gsl_rng_free(r);
	//result = result+10-3 I
	for (int i = 0; i < mat_size; i++) {
        gsl_matrix_set(result, i, i, gsl_matrix_get(result, i, i) + 1e-5);
    }
	/*std::cout << j<<"diag"<< std::endl;
	print_matrix(diag);
    //test
	print_matrix(result);*/
	//symmetric and positive-definite for result
    // (B_gamma + \Sigma_0^-1) = L L^T
    gsl_linalg_cholesky_decomp1(result);
	//std::cout << j<<"decomp result"<< std::endl;
	//print_matrix(result);
    // \mu = L^{-1} A_vec
    gsl_blas_dtrsv(CblasLower, CblasNoTrans, \
	    CblasNonUnit, result, A_vec);
	
    // N(\mu, I)
    gsl_blas_daxpy(1.0, A_vec, beta_c);
	
    // X ~ N(\mu, I), L^{-T} X ~ N( L^{-T} \mu, (L L^T)^{-1} )
    gsl_blas_dtrsv(CblasLower, CblasTrans, \
	    CblasNonUnit, result, beta_c);
	
	double max_b = 0.0;
	for (size_t i=0; i<sub_size; i++) {
		beta1[causal_list[i]]= gsl_vector_get(beta_c,i);
		//std::cout << "beta_c" << gsl_vector_get(beta_c,i) << std::endl;
		beta2[causal_list[i]]= gsl_vector_get(beta_c,sub_size+i);
		//std::cout << "beta_c" << gsl_vector_get(beta_c,i)<<" "<<gsl_vector_get(beta_c,sub_size+i) << std::endl;
		if(gsl_vector_get(beta_c,i)>max_b){
			max_b = gsl_vector_get(beta_c,i);
		}
		if(gsl_vector_get(beta_c,sub_size+i)>max_b){
			max_b = gsl_vector_get(beta_c,sub_size+i);
		}
		beta3[causal_list[i]]= gsl_vector_get(beta_c, 2*sub_size+i);
		//std::cout << "beta_c" << gsl_vector_get(beta_c,sub_size+i) << std::endl;
    } 
	
	//update the res here
	/**/
	//std::cout << "A before update: " << gsl_vector_get(A, 2) << std::endl;
	gsl_blas_dgemv(CblasTrans, -1.0, X1, &beta1_here.vector, 1.0, A_res);
	//std::cout << "A after X1 update: " << gsl_vector_get(A, 2) << std::endl;
	gsl_blas_dgemv(CblasTrans, -1.0, X2, &beta2_here.vector, 1.0, A_res);
	//std::cout << "A after X2 update: " << gsl_vector_get(A, 2) << std::endl;
	
	double max_a = 0.0;
	for (size_t i = 0; i < dat->n_ind; i++) {
		residual[i] = gsl_vector_get(A_res, i);
		if(residual[i]>max_a){
			max_a = residual[i];
		}
	}
	if(max_b>1){
		cout<<"max_beta "<<max_b<<" at "<<j<<"and max res "<<max_a<<"tau "<<tau <<endl;
	}
	
	gsl_vector_free(A_vec);
	gsl_vector_free(beta_c);
	gsl_vector_free(A_res);
	
    gsl_matrix_free(diag);
	
    gsl_matrix_free(var_mat);
	
    //gsl_permutation_free(p);
	
	gsl_matrix_free(ptr);
	gsl_matrix_free(X1);
	gsl_matrix_free(X2);
	gsl_matrix_free(result);
	//cout<<"free "<<endl;
}

void MCMC_state::sample_assignment(size_t j, const Dat *dat, const ldmat_data &ldmat_dat,unsigned sz){

    size_t start_i = dat->boundary[j].first;
    size_t end_i = dat->boundary[j].second;
	int num_snp = end_i-start_i;
	
    gsl_matrix*log_prob = gsl_matrix_alloc(n_cluster,num_snp);
    
    //Bj = (tau*b1j,tau*b2j,N3*b3j)
	gsl_matrix* A = gsl_matrix_alloc(3, 3);
    gsl_vector* B_i = gsl_vector_alloc(3);
    gsl_matrix* inv_A = gsl_matrix_alloc(3, 3);
    gsl_permutation* perm = gsl_permutation_alloc(3);
    int signum;
	
	//std::cout << num_snp<< std::endl;
	for (int k = 1; k < n_cluster; ++k) {
        double deno_k = cluster_var[k]*det_sigma;
		double ck1 = (rho*rho_a-rho_b)/deno_k;
        double ck2 = (rho*rho_b-rho_a)/deno_k;
		double ck3_base = (rho_a*rho_b-rho)/deno_k;

        // Compute ak1, ak2, ak3, ck1, ck2, ck3 (using loops and GSL operations)
        for (int i = 0; i < num_snp; ++i) {
            // Fill matrix A_i
			
            gsl_matrix_set(A, 0, 0, tau*eta*eta*dat->geno1_sq[i]+(1-rho_a*rho_a)/deno_k);//ak1 
			gsl_matrix_set(A, 1, 1, tau*eta*eta*dat->geno2_sq[i]+(1-rho_b*rho_b)/deno_k);//ak2 
			double Bii = gsl_matrix_get(ldmat_dat.B[j], i, i);//ak3 = N3*eta^2*B3,ii+...
            gsl_matrix_set(A, 2, 2,sz*eta*eta*Bii+(1-rho*rho)/deno_k );//ak3 
			
			double ck3 = ck3_base+tau*eta*eta*dat->geno12_prod[i];
            gsl_matrix_set(A, 0, 1, ck3);
            gsl_matrix_set(A, 0, 2, ck2);
            gsl_matrix_set(A, 1, 0, ck3);
            gsl_matrix_set(A, 1, 2, ck1);
            gsl_matrix_set(A, 2, 0, ck2);
            gsl_matrix_set(A, 2, 1, ck1);
			
			//cout<<j<<"at i:"<<i<<"var"<<cluster_var[k]<<"tau"<<tau<<"eta"<<eta<<"n " <<dat->n_ind<<endl;
			//print_matrix(A);

            // Compute B_i vector
            gsl_vector_set(B_i, 0, b1[start_i+i]);
            gsl_vector_set(B_i, 1, b2[start_i+i]);
            gsl_vector_set(B_i, 2, sz * b3[start_i+i]);
			
            // Compute inverse of A
            gsl_linalg_LU_decomp(A, perm, &signum);
            gsl_linalg_LU_invert(A, perm, inv_A);
			
            // exp_ele = 0.5 * B_i.T * inv_A * B_i
            double exp_ele;
            gsl_blas_dgemv(CblasNoTrans, 1.0, inv_A, B_i, 0.0, B_i); // inv_A * B_i
            gsl_blas_ddot(B_i, B_i, &exp_ele); // B_i.T * (inv_A * B_i)
            exp_ele *= 0.5;

            // non_exp = -0.5 * log(det(A)) - 1.5 * log(var_k) - 0.5 * log(det_sigma) + log(p1/p0) + log(pi_k)
            double non_exp = -0.5 * std::log(gsl_linalg_LU_det(A, signum)) - 1.5 * std::log(cluster_var[k]) - 0.5 * std::log(det_sigma);
            non_exp += std::log((1-p_0) / p_0) + std::log(p[k] + 1e-40);
            // log_prob[k, i] = exp_ele + non_exp
            gsl_matrix_set(log_prob, k, i, exp_ele + non_exp);
        }
    }

	float *rnd = (float *)malloc(num_snp * sizeof(float));
	for (int i = 0; i < num_snp; ++i) {
		float log_exp_sum = logf(0.0f); // Initialize to -infinity
		for (size_t k = 0; k < n_cluster; k++) {
			log_exp_sum = logf(expf(log_exp_sum) + expf(gsl_matrix_get(log_prob, k, i)));
		}
		assgn[i+start_i] = n_cluster-1;
		
		//generate a uniform variable
        rnd[i] = (float)rand() / RAND_MAX;
		for (size_t k=0; k<n_cluster-1; k++) {
			rnd[i] -= expf(gsl_matrix_get(log_prob, k, i) - log_exp_sum);
			if (rnd[i] < 0) {
			assgn[i+start_i] = k;
			break;
			}
		}
		
	}
	//std::cout << j<<" assign done." << std::endl;
	gsl_matrix_free(A);
    gsl_matrix_free(inv_A);
    gsl_vector_free(B_i);
    gsl_permutation_free(perm);
	gsl_matrix_free(log_prob);
	free(rnd); 
}

void solve_ldmat(const Dat *dat, ldmat_data &ldmat_dat, const double a, unsigned sz) {
    for (size_t i=0; i<dat->ref_ld_mat.size(); i++) {
		size_t size = dat->boundary[i].second - dat->boundary[i].first;
		gsl_matrix *A = gsl_matrix_alloc(size, size);
		gsl_matrix *B = gsl_matrix_alloc(size, size);
		gsl_matrix *L = gsl_matrix_alloc(size, size);
		gsl_matrix_memcpy(A, dat->ref_ld_mat[i]);
		gsl_matrix_memcpy(B, dat->ref_ld_mat[i]);
		gsl_matrix_memcpy(L, dat->ref_ld_mat[i]);
		
		// (R + aNI) / N A = R via cholesky decomp
		// replace aN with a=0.1
		gsl_vector_view diag = gsl_matrix_diagonal(B);
		gsl_vector_add_constant(&diag.vector, a);
		

		gsl_linalg_cholesky_decomp1(B);
		gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, B, A);
		gsl_blas_dtrsm(CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1.0, B, A);

		// B = RA
		// Changed May 21 2021 as A may not be symmetric
		//gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, L, A, 0, B);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, L, A, 0, B);
		
		// L = R %*% R;
		gsl_matrix_mul_elements(L, L);

		// memory allocation for A^T beta_mrg
		// Changed May 21 2021 from A to A^T
		gsl_vector *beta_mrg = gsl_vector_alloc(size);
		for (size_t j=0; j<size; j++) {
			gsl_vector_set(beta_mrg, j, dat->beta_mrg[j+dat->boundary[i].first]);
		}
		gsl_vector *b_tmp = gsl_vector_alloc(size);

		// Changed May 21 2021 from A to A^T
		//b_tmp = A^T * beta_mrg
		gsl_blas_dgemv(CblasTrans, 1.0, A, beta_mrg, 0, b_tmp);

		ldmat_dat.A.push_back(A);
		ldmat_dat.B.push_back(B);
		ldmat_dat.L.push_back(L);
		ldmat_dat.calc_b_tmp.push_back(b_tmp);
		ldmat_dat.beta_mrg.push_back(beta_mrg);
		/*for(size_t j = 0; j < 10; j++){
			for(size_t k=0;k<10;k++){
				cout<< gsl_matrix_get(B,j,k)<<" ";
			}
			cout<<endl;
		}*/
		//print B3
    }
}

void MCMC_state::compute_h2(const Dat *dat) {

    h2_1 = 0.0;
	std::vector<double> result(dat->n_ind, 0.0);
    
    h2_1 = gsl_stats_variance(G, 1, dat->n_ind) / gsl_stats_variance(dat->pheno, 1, dat->n_ind);
    
    double h2_tmp = 0;
    h2_2 = 0;

    if (beta3 == nullptr) {
        cerr << "Error: beta3 is null" << endl;
        return;
    }

    for (size_t j = 0; j < dat->ref_ld_mat.size(); j++) {
        size_t start_i = dat->boundary[j].first;
        size_t end_i = dat->boundary[j].second;
        gsl_vector *tmp = gsl_vector_alloc(end_i - start_i);
    
        gsl_vector_view beta_view = gsl_vector_view_array(beta3, dat->n_snp);
        gsl_vector_view beta_j = gsl_vector_subvector(&beta_view.vector, start_i, end_i - start_i);

        gsl_blas_dsymv(CblasUpper, 1.0, dat->ref_ld_mat[j], &beta_j.vector, 0, tmp);
        gsl_blas_ddot(tmp, &beta_j.vector, &h2_tmp);
        h2_2 += h2_tmp;

        gsl_vector_free(tmp);
    }
}

void mcmc(Dat *dat, std::string out_path, int iter, int burn, double maf,unsigned sz) {
	using namespace std::chrono;
    MCMC_state state;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    std::random_device rd;
    std::mt19937 gen(rd());
    double *res_beta1 = (double *) calloc(dat->n_snp, sizeof(double));
    double *res_beta2 = (double *) calloc(dat->n_snp, sizeof(double));
    double *res_beta3 = (double *) calloc(dat->n_snp, sizeof(double));
	for(int i;i<dat->n_snp;i++){
		res_beta1[i]=0.0;
		res_beta2[i]=0.0;
		res_beta3[i]=0.0;
	}
	

    init_state(dat, &state);
	//std::cout << "Initial state." << std::endl;
	ldmat_data ldmat_dat;

    // set up OLS env
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(dat->n_ind, dat->n_cov);
    gsl_vector *y = gsl_vector_alloc(dat->n_ind);
    gsl_matrix *W = gsl_matrix_alloc(dat->n_ind, dat->n_cov);
    if (dat->n_cov != 0) {
	for (size_t i=0; i<dat->n_ind; i++) {
	    for (size_t j=0; j<dat->n_cov; j++) {
		gsl_matrix_set(W, i, j, dat->covar[j][i]);
	    }
	}
    }
    gsl_matrix *cov = gsl_matrix_alloc(dat->n_cov, dat->n_cov);
    gsl_vector *alpha = gsl_vector_alloc(dat->n_cov);	    
    double chisq = 0;
	int n_total = iter - 1 -  burn;
	// OLS
	if (dat->n_cov != 0) {
		for (size_t i=0; i<dat->n_ind; i++) {
		gsl_vector_set(y, i, dat->y[i] );
		}
		//y = Y  and then regress on W to get a
		gsl_multifit_linear(W, y, alpha, cov, &chisq, work);
		//update y=aW
		gsl_blas_dgemv(CblasNoTrans, 1.0, W, alpha, 0, y);	
		// update pheno and residual
		for (size_t i=0; i<dat->n_ind; i++) {
		dat->pheno[i] = dat->y[i] - gsl_vector_get(y, i);//pheno = Y-aW
		}
	}
	gsl_multifit_linear_free(work);
	gsl_vector_free(y);
	gsl_matrix_free(W);
	gsl_matrix_free(cov);
	gsl_vector_free(alpha);
	// update suffstats
	for (size_t k=0; k<state.n_cluster; k++) {
			state.suffstats[k] = 0;
			state.sumsq[k] = 0;
	}
	
	//solve_ldmat(dat, ldmat_dat, 0.1, sz);//a=0.1
	std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed;
	int ref_size =  dat->ref_ld_mat.size();
	start = std::chrono::high_resolution_clock::now();
	solve_ldmat(dat, ldmat_dat, 0.1, sz);//a=0.1
	state.update_suffstats();
	state.update_residual(dat);
	
    for (size_t n_mcmc=0; n_mcmc<iter; n_mcmc++) {
		
		state.sample_sigma2();
		state.sample_tau(dat);
		//std::cout << "Sigma2. "<<endl;
		for (size_t i=0; i<ref_size; i++) {
			//std::cout << "Starting iteration " << i << " of " << ref_size << std::endl;
			//cout<<i<<"begin"<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
			state.calc_b(i, dat, ldmat_dat);
			//cout<<i<<" b"<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
			//std::cout << "Finished Calculate b for " << i << std::endl;
			state.sample_assignment(i,dat,ldmat_dat,sz);
			
			//cout<<i<<"assig "<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
			//std::cout << "Finished Sample assign for " << i << std::endl;
			state.sample_beta(i, dat, ldmat_dat,sz);
			//std::cout << "Finished Sample beta for " << i << std::endl;
			//cout<<i<<"beta "<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
			
		}
		
		state.update_suffstats();
		state.update_p0(dat);
		state.update_residual(dat);
		//cout<<"suff "<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
		//cout << "Update suffstats" << endl;
		state.sample_V();
		//cout<<"V"<<gsl_matrix_get(ldmat_dat.B[8],1,10)<<endl;
		
		state.update_p();
		
		//cout << "Sample p." << endl;
		
		
		//std::cout << "Sample tau." << std::endl;
		//auto start = steady_clock::now();
		//state.sample_eta(ldmat_dat);
		
		if (n_mcmc % 2 == 0) {
			state.compute_h2(dat);
			
			cout << n_mcmc << " iter. h2_1: " << state.h2_1*square(state.eta) <<" h2_2: " << state.h2_2*square(state.eta);
			end = high_resolution_clock::now();
			elapsed = end - start;
			cout <<" "<< elapsed.count() << " seconds." ;
			start = high_resolution_clock::now();

			double max_beta1 = 0, max_beta2 = 0, max_beta3 =0;
			for (size_t i=0; i<dat->n_snp; i++) {
				if (state.eta * state.beta1[i] > max_beta1) {
				max_beta1 = state.eta * state.beta1[i];
				}
				if (state.eta * state.beta2[i] > max_beta2) {
				max_beta2 = state.eta * state.beta2[i];
				}
				if (state.eta * state.beta3[i] > max_beta3) {
				max_beta3 = state.eta * state.beta3[i];
				}
			}
			cout << " max beta1: " << std::to_string(max_beta1) \
			<<" max beta2: " << std::to_string(max_beta2) \
			<<" max beta3: " << std::to_string(max_beta3) << endl; 
				
		}
		// record the beta
		if (n_mcmc > burn) {
			for (size_t i=0; i<dat->n_snp; i++) {
				res_beta1[i] += state.beta1[i] *state.eta/n_total;
				res_beta2[i] += state.beta2[i] *state.eta/n_total;
				res_beta3[i] += state.beta3[i] *state.eta/n_total;
			}
		}
		
		//std::cout <<n_mcmc<<endl;
	}
	
	std::ofstream out(out_path);
	for (size_t i=0; i<dat->n_snp; i++) {
		out << dat->chr[i] << "\t" << dat->pos[i] << "\t" << \
		dat->rsid[i] << "\t" << dat->ref[i] << "\t" << \
		dat->alt[i] << "\t" << res_beta1[i] << \
		"\t" << res_beta2[i]<< "\t"<<res_beta3[i]<< endl;
	}
	out.close();
		
	destroy_state(&state);
	gsl_rng_free(r);
	free(res_beta1);
	free(res_beta2);
	free(res_beta3);

	for (size_t i=0; i<dat->ref_ld_mat.size(); i++) {
	gsl_matrix_free(ldmat_dat.A[i]);
	gsl_matrix_free(ldmat_dat.B[i]);
	gsl_matrix_free(ldmat_dat.L[i]);
	gsl_vector_free(ldmat_dat.calc_b_tmp[i]);
	gsl_vector_free(ldmat_dat.beta_mrg[i]);
	gsl_matrix_free(dat->ref_ld_mat[i]);
    } 
}

int main(int argc, char *argv[]) {
    
	//start time
	auto start = std::chrono::high_resolution_clock::now();
	Dat dat;

    std::string pheno_path, geno1_path, geno2_path, vcf_path, msp_path, out_path, covar_path, ss_path, ref_dir,valid_path,bim_path;
    //int make_ref = 0, run_mcmc = 0;
    int opt_llk = 1, N=0,n_threads = 1, chr = 0;
	double rho1,rho2,rho3;
    int i = 1, iter = 1000, burn = 500;
    
    while (i < argc) {
	if (strcmp(argv[i], "-pheno") == 0) {
	    pheno_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-geno1") == 0) {
	    geno1_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-geno2") == 0) {
	    geno2_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-vcf") == 0) {
	    vcf_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-msp") == 0) {
	    msp_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-iter") == 0) {
	    iter = std::stoi(argv[i+1]);
	    i +=2;
	}
	else if (strcmp(argv[i], "-burn") == 0) {
	    burn = std::stoi(argv[i+1]);
	    i += 2;
	}
	else if (strcmp(argv[i], "-out") == 0) {
	    out_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-covar") == 0) {
	    covar_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-ss") == 0) {
	    ss_path = argv[i+1];
	    i += 2;
	}
    else if (strcmp(argv[i], "-ref_dir") == 0) {
	    ref_dir = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-valid") == 0) {
	    valid_path = argv[i+1];
	    i += 2;
	}
    else if (strcmp(argv[i], "-N") == 0) {
	    N = std::stoi(argv[i+1]);
	    if ( N <= 0) {
		cout << "Incorrect N: " << argv[i+1] << endl;
		return 0;
	    }
	    if (N <= 1000) {
		cout << "Warning: sample size too small, might" \
		    " not achieve good performance." << endl;
	    }
	    i += 2;
	}
    else if (strcmp(argv[i], "-opt_llk") == 0) {
	    opt_llk = std::stoi(argv[i+1]);
	    if (opt_llk != 1 && opt_llk != 2) {
		cout << "opt_llk must be in 1 or 2." << endl;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-chr") == 0) {   
	    chr = std::stoi(argv[i+1]);
	    if ( chr > 22 || chr < 0) {
		cout << "Incorrect chromosome: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
    else if (strcmp(argv[i], "-bim") == 0) {
	    bim_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-rho1") == 0) {
        rho1 = std::stod(argv[i+1]);
        i += 2;
    }
	else if (strcmp(argv[i], "-rho2") == 0) {
        rho2 = std::stod(argv[i+1]);
        i += 2;
    }
    else if (strcmp(argv[i], "-rho3") == 0) {
        rho3 = std::stod(argv[i+1]);
        i += 2;
    }
	else {
	    cout << "Invalid option: " << argv[i] << endl;
	    return 0;
	}
    }

    get_size_vcf(pheno_path.c_str(), vcf_path.c_str(), &dat);
    
    //read_lanc(vcf_path.c_str(), msp_path.c_str(), &dat);

    read_pheno(pheno_path.c_str(), &dat);

    if (!covar_path.empty()) {
	read_cov(covar_path.c_str(), &dat);
    }
	
    //linear(&dat, out_path.c_str());
	string ref_ldmat = ref_dir + "/chr" + \
	       std::to_string(chr) + ".dat";
	string ref_snpinfo = ref_dir + "/chr" + \
	       std::to_string(chr) + ".snpInfo";

    coord(vcf_path.c_str(), msp_path.c_str(),bim_path.c_str(),ref_snpinfo, ss_path.c_str(), valid_path.c_str(), ref_ldmat, &dat, N, opt_llk);
	
    double maf = 0;

    //check_maf(&dat, maf);
    
    prod_geno(&dat);

	dat.rho_1 = rho1;
	dat.rho_2 = rho2;
	dat.rho_3 = rho3;
	std::cout <<"rho:"<<rho1<<" "<<rho2<<" "<<rho3<< " "<< std::endl; 

	// End time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate elapsed time
    std::chrono::duration<double> elapsed = end - start;

    // Print the time
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;


	//save_geno(&dat);


    maf = 0.01;
	ldmat_data ldmat_dat;

    
	cout<<"Repeat "<< iter<<" with burn in "<<burn<<endl;
    mcmc(&dat, out_path.c_str(), iter, burn, maf,N);

    for (size_t i=0; i<dat.n_snp; i++) {
		free(dat.geno1[i]);
		free(dat.geno2[i]);
    }
    free(dat.geno1);
    free(dat.geno2);
    if (!covar_path.empty()) {
		free(dat.covar);
    }
    free(dat.maf1);
    free(dat.maf2);
    free(dat.n_anc1);
    free(dat.n_anc2);
    free(dat.y);
    free(dat.pheno);
    free(dat.geno1_sq);
    free(dat.geno2_sq);
    free(dat.geno12_prod);
    return 0;
}

