#include "parse_gen.h"
#include "regress.h"
#include "gsl/gsl_rng.h"
typedef struct {
    std::vector<gsl_matrix*> A;
    std::vector<gsl_matrix*> B;
    std::vector<gsl_matrix*> L;
    std::vector<gsl_vector*> beta_mrg;
    std::vector<gsl_vector*> calc_b_tmp;
    size_t n_blocks;  // Number of blocks (size of the outer vectors)
    size_t* block_sizes;
} ldmat_data;

typedef struct {
public:
    double eta = 1.0;
    double h2_1,h2_2;
    double rho ,rho_a, rho_b ;
    double inv[3][3];
    double det_sigma;
    //double sigmae2 = 0.14;
    double sigmae2 = 1.0;
    double a0k = 0.5, b0k = 0.5;
    size_t n_cluster;
    size_t n_snp;
    double tau;
    double *beta1;
    double *beta2;
    double *beta3;
    double *residual;
    double *G;
    int *assgn;
    double *cluster_var;
    double p_0;
    double *p;
    double *log_p;
    double *V;
    int *suffstats;
    double *sumsq;
    //double pi_pop[4] = {0.8, 0.05, 0.05, 0.1};
    //double pi_pop[4] = {0.25, 0.25, 0.25, 0.25};
    double alpha = 1.0;
    /*double *aj1;
    double *aj2;
    double *cj;
    double *mu_j1;
    double *mu_j2;*/
    double *b1;
    double *b2;
    double *b3;
    double bj1 = 0;
    double bj2 = 0;
    //gsl_rng *r;
    void sample_sigma2();
	void calc_b(size_t j, const Dat *dat, const ldmat_data &ldmat_dat);
	void sample_assignment(size_t j, const Dat *dat, \
		        const ldmat_data &ldmat_dat,unsigned sz);
	void update_suffstats();
	void sample_V();
	void update_p();
    void update_p0(const Dat *dat);
	void sample_alpha();
	void sample_beta(size_t j, const Dat *dat, \
                const ldmat_data &ldmat_dat,unsigned sz);
    void update_residual(const Dat *dat);
    void sample_tau(const Dat *dat);
	void compute_h2(const Dat *dat);
    /*
	void sample_eta(const ldmat_data &ldmat_dat);*/
    /*private:
	double a0k;
	double b0k;
	size_t M, n_snp;
	gsl_rng *r;*/
} MCMC_state;



typedef struct {
    double **beta1;
    double **beta2;
    double **beta3;
} MCMC_samples;

void mcmc(Dat *dat, std::string out_path, int iter, int burn, double maf);
