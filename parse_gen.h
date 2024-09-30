#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "gsl/gsl_matrix.h"


#ifndef PARSE_GEN_H
#define PARSE_GEN_H

using std::vector;using std::string;

typedef struct {
    double **geno1;
    double **geno2;
    double *pheno;
    double *y;
    std::vector<size_t> ind_idx;
    size_t n_snp;
    size_t n_ind;
    size_t n_cov = 0;
    double *geno1_sq;
    double *geno2_sq;
    double *geno12_prod;
    double **covar;
    double **proj;
    vector<string> chr;
    vector<string> id;
    vector<string> pos;
    vector<string> rsid;
    vector<string> ref;
    vector<string> alt;
    double *maf1;
    double *maf2;
    double *n_anc1;
    double *n_anc2;
    vector<string> id_st;
    vector<string> A1;
    vector<string> A2;
    vector<double> beta_mrg;
    vector<std::pair<size_t, size_t>> boundary;
    vector<gsl_matrix *> ref_ld_mat;
    vector<string> snplist;
    vector<double> sz;
    vector<int> array;
    double rho_1,rho_2,rho_3;
} Dat;

typedef struct {
    string A1;
    string A2;
    string rsid;
    bool include_ref;
    bool include_ss;
    bool include_geno;
    size_t geno_idx;
    double beta;
    int array;
    double sz;
    bool flip;
} CoordInfo;

#endif

void get_size_vcf(const string &pheno_path, const string &geno_path, Dat *dat);

void read_pheno(const string &pheno_path, Dat *dat);

void read_cov(const string &cov_path,  Dat *dat);

//void read_lanc(const string &vcf_path,  const string &msp_path, Dat *dat,unordered_map<string, string> &bim_dict);

void check_maf(Dat *dat, double maf);

void prod_geno(Dat *dat);

void save_geno(Dat *dat);
 
void coord(const string &vcf_path,  const string &msp_path, const string &bim_path,const string &ref_path, const string &ss_path, const string &valid_path,const string &ldmat_path, Dat *dat, unsigned sz, int opt_llk);
