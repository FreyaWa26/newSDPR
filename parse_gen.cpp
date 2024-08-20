#include "parse_gen.h"
#include <assert.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <unordered_set> 
#include <stdexcept>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_blas.h"
#include "math.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <gsl/gsl_matrix.h>
//#include "H5Cpp.h"
//using namespace H5;

using std::cout; using std::endl; using std::ifstream;
using std::string; using std::getline; using std::istringstream;
using std::vector; using std::unordered_map; 
using std::cout; using std::endl;
using std::pair; using std::find;

// get n_snp and n_ind
void get_size_vcf(const string &pheno_path, const string &geno_path, Dat *dat) {
    size_t n_ind = 0;
    size_t n_invalid = 0;
    size_t n_snp = 0;

    string line;
    ifstream infile1(pheno_path.c_str());
    string id;
    string y;
    size_t i = 0;
    while (infile1 >> id >> id >> y) {
	try {
	    std::stod(y); 
	    dat->ind_idx.push_back(i);
	    n_ind++;
	}
	catch (std::invalid_argument&) {
	    n_invalid++;
	}
	i++;
    }
    dat->n_ind = n_ind;
    cout << "Warning: " + std::to_string(n_invalid) + \
	" individuals with invalid phenotypes." << endl;

    ifstream infile2(geno_path.c_str());
    while(getline(infile2, line)) {
	if (line.find("##") == 0) {
	    continue;
	}
	else if (line.find("#") == 0) {
	    continue;
	}
	else {
	    n_snp++;	
	}
    }
    dat->n_snp = n_snp;
    cout << "In total " + std::to_string(n_snp) + " SNPs and " \
	+ std::to_string(n_ind) + " individuals to be readed." << endl;
}

void read_pheno(const std::string &pheno_path, Dat *dat) {
    ifstream infile(pheno_path.c_str());

    cout << "Reading phenotype file from: " + pheno_path + "." << endl;

    string id;
    string y;
    double *pheno = (double *) malloc(dat->n_ind*sizeof(double));

    size_t i = 0, idx = 0;
    while (infile >> id >> id >> y) {
	if (i == dat->ind_idx[idx]) {
	    pheno[idx] = stod(y);
	    idx++;
	}
	i++;
    }
    dat->pheno = pheno;

    cout << "Readed phenotype from " + std::to_string(idx) + " individuals." << endl;
}

void read_cov(const std::string &cov_path,  Dat *dat) {
    ifstream infile(cov_path.c_str());

    cout << "Reading covariate file from: " + cov_path + "." << endl;

    size_t n_cov = 0, i = 0, idx = 0;

    string line, token;

    while (getline(infile, line)) {
	n_cov = 0;
	std::istringstream iss(line);
	if (i == 0) {
	    while (getline(iss, token, '\t')) {   
		n_cov++;
	    }
	    cout << "Reading " + std::to_string(n_cov) + " covariates." << endl;
	    dat->n_cov = n_cov;
	    dat->covar = (double **) malloc(n_cov*sizeof(double *));
	    for (size_t k=0; k<n_cov; k++) {
		dat->covar[k] = (double *) malloc(dat->n_ind*sizeof(double));
	    }
	    n_cov = 0;
	    iss.clear();
	    iss.str(line);
	}
	if (i == dat->ind_idx[idx]) {
	    while (getline(iss, token, '\t')) {
		dat->covar[n_cov][idx] = stod(token);
		n_cov++;
	    }
	    idx++;
	}
	i++;
    }
    infile.close();
}
void parse_bim(Dat *dat, const string &bim_path,unordered_map<string, string> &bim_dict){
    ifstream infile(bim_path.c_str());
    int n_dup= 0;
    if (!infile) {
	throw std::runtime_error("Error: cannot open .bim file: " + bim_path);
    }
    cout << "Reading bim file from: " + bim_path + "." << endl;

    string line;
    int n =0;
   
    // Read the .bim file line by line
    while (getline(infile, line)) {
        istringstream ss(line);
        string chr, snpid, dummy1, geno_pos,dummy2,dummy3;

        // Read columns
        if (!(ss >> chr >> snpid >> dummy1 >> geno_pos >> dummy2 >> dummy3)) {
            std::cerr << "Error parsing .bim file line: " << line << std::endl;
            continue;
        }
        
        // Store the position and SNPID in the map
        if (!bim_dict.insert(pair<string, string>( geno_pos,snpid)).second) {
    		n_dup++;
            continue; 
	    }
        n = n+1;
    }
    infile.close();
    cout << "Duplicate "<<n_dup<<" snps from bim file."<< endl;
    cout << "Read " <<n <<" snps from bim file." << endl;     
}

void read_lanc(const std::string &vcf_path,  const std::string &msp_path, Dat *dat, unordered_map<string, string> &bim_dict,unordered_map<string, CoordInfo*> &ref_dict) {
    
    string line1, line2;
    string token1, token2;string ref,alt;
    size_t n_match =0, n_bad=0, n_flip=0; 
	unordered_map<string, CoordInfo*>::iterator it;

    ifstream mspfile(msp_path.c_str());
    cout << "Reading RFmix msp file from: " + msp_path + "." << endl;
    
    ifstream infile(vcf_path.c_str());
    cout << "Reading VCF file from: " + vcf_path + "." << endl;

    // skip first two lines of msp file
    getline(mspfile, line1);
    getline(mspfile, line1);

    // skip the header of vcf file
    int n_ind = 0;
    while (getline(infile, line2)) {
	if (line2.find("##") == 0) {
	    continue;
	}
	else if (line2.find("#") == 0) {
	    int idx_2 = 0;
	    std::istringstream iss2(line2);
	    while (getline(iss2, token2, '\t')) {
		idx_2++;
	    }
	    n_ind = idx_2 - 9;
	}
	else {
	    break;
	}
    }

    int *hap_lanc = (int *) malloc(2*n_ind*sizeof(int));
    unsigned spos = 0, epos = 0, pos = 0;
    int chr_vcf = 1, chr_msp = 1;
    
    // read the pos from the first line of vcf file
    std::istringstream iss2(line2);
    int idx2 = 0;
    for (; idx2<9; idx2++) {
	getline(iss2, token2, '\t');
	if (idx2 == 0) {
	    chr_vcf = std::stoi(token2);
	}
	if (idx2 == 1) {
	    pos = std::stoul(token2);
	}
    }

    size_t idx_snp = 0,idx_read=0;
    dat->geno1 = (double **) malloc(dat->n_snp*sizeof(double*));
    dat->geno2 = (double **) malloc(dat->n_snp*sizeof(double*));
    dat->n_anc1 = (double *) calloc(dat->n_snp, sizeof(double));
    dat->n_anc2 = (double *) calloc(dat->n_snp, sizeof(double));
    for (size_t i=0; i<dat->n_snp; i++) {
	dat->geno1[i] = (double *) calloc(dat->n_ind, sizeof(double));
	dat->geno2[i] = (double *) calloc(dat->n_ind, sizeof(double));
	if (!dat->geno1[i] || !dat->geno2[i]) {
	    cout << "Error: memory allocation failed for" + \
		std::to_string(i) + " th SNP.";
	    exit(EXIT_FAILURE);
	}
    }

    while (getline(mspfile, line1)) {
	// read msp file
	std::istringstream iss1(line1);
	for (int idx1=0; idx1<2*n_ind+6; idx1++) {
	    getline(iss1, token1, '\t'); 
	    if (idx1 == 0) {
		chr_msp = std::stoi(token1);
	    }
	    else if (idx1 == 1) {
		spos = std::stoul(token1);
	    }
	    else if (idx1 == 2) {
		epos = std::stoul(token1);
	    }
	    else if (idx1 >= 6) {
		hap_lanc[idx1-6] = std::stoi(token1);
		if (hap_lanc[idx1-6] != 0 && hap_lanc[idx1-6] != 1) {
		    cout << "RFmix field must be either 0 or 1." << endl;
		    return;
		}
	    }
	}

	if ((chr_vcf != chr_msp) || (pos != spos)) {
	    cout << "Inconsistent starting position: chr_vcf: " + std::to_string(chr_vcf) + \
		" chr_msp: " + std::to_string(chr_msp) + " pos: " + \
		std::to_string(pos) + " spos: " + std::to_string(spos) << endl;
	    exit(EXIT_FAILURE);
	} 
	   
	// read vcf file
	while ((chr_vcf == chr_msp && pos >= spos && pos < epos) || idx_read == dat->n_snp-1) {

	    if (idx_read == dat->n_snp-1) {
		assert(chr_vcf == chr_msp && pos == epos);
	    }

	    // reset the stream
	    std::istringstream iss2(line2);

	    size_t k = 0;
		bool flip = false;//reset if the genotype need flip
        for (idx2=0; idx2<n_ind+9; idx2++) {
		getline(iss2, token2, '\t');

		if (idx2 == 0) {
		    dat->chr.push_back(token2);
		}

		if (idx2 == 1) {
		    dat->pos.push_back(token2);
			if (bim_dict.find(token2) != bim_dict.end()) {
                // coordination
                string snp_id = bim_dict[token2];
        	    it = ref_dict.find(snp_id);
        	    if (it != ref_dict.end() && it->second->include_ref) {
                    idx2++;getline(iss2, token2, '\t');//idx2=2
            		dat->id.push_back(token2);
                    idx2++;getline(iss2, ref, '\t');//idx2=3
                    idx2++;getline(iss2, alt, '\t');//idx2=4
            		if (ref == it->second->A1 && alt == it->second->A2) {
						dat->ref.push_back(ref);
						dat->alt.push_back(alt);
						dat->rsid.push_back(snp_id);
                        it->second->include_geno = true;
						it->second->geno_idx = idx_snp;
						idx_snp++;
                    }
					else if (ref == it->second->A2 && alt == it->second->A1){
						dat->ref.push_back(ref);
						dat->alt.push_back(alt);
						dat->rsid.push_back(snp_id);
						it->second->include_geno = true;
						it->second->geno_idx = idx_snp;
						
						idx_snp++;n_flip++;
						it->second->flip=true;//need to flip LD-mat later
						it->second->A1 = ref;
						it->second->A2 = alt;
					}
					else{
						n_bad++;
						dat->chr.pop_back();
                		dat->pos.pop_back();
						break;//skip this line/snp
					}
				}
			}
		}
		if (idx2 >= 9) {

		    // individuals kept in analysis
		    if (idx2-9 != dat->ind_idx[k]) {
			continue;
		    }

		    // check phasing
		    if (token2[1] != '|') {
			cout << "Genotype must be phased." << endl;
			return;
		    }
		    
		    // read the genotype
		    if (token2[0] == '.' || token2[3] == '.') {
			cout << "Missing genotype not supported yet." << endl;
			return;
		    }
			if (hap_lanc[2*(idx2-9)] == 0) {
			dat->geno1[idx_snp-1][k] += std::stod(&token2[0]);
			dat->n_anc1[idx_snp-1]++;
		    }
		    else {
			dat->geno2[idx_snp-1][k] += std::stod(&token2[0]);
			dat->n_anc2[idx_snp-1]++;
		    }

		    if (hap_lanc[2*(idx2-9)+1] == 0) {
			dat->geno1[idx_snp-1][k] += std::stod(&token2[2]);
			dat->n_anc1[idx_snp-1]++;
		    }
		    else  {
			dat->geno2[idx_snp-1][k] += std::stod(&token2[2]);
			dat->n_anc2[idx_snp-1]++;
		    }
			
		    k++;
		}
	    }
	   
	    // read the next line of vcf and update pos
	    assert(idx2 == n_ind+9);
	    
	    if (getline(infile, line2)) {
		iss2.clear();
		iss2.str(line2);
		for (idx2=0; idx2<9; idx2++) {
		    getline(iss2, token2, '\t');
		    if (idx2 == 0) {
			chr_vcf = std::stoi(token2);
		    }
		    if (idx2 == 1) {
			pos = std::stoul(token2);
		    }
		}
		idx_read++;
	    }
	    else {
		break;
	    }
		
    }
	}
    cout << "Read " << std::to_string(idx_read+1) << \
	" SNPs from " << std::to_string(dat->n_ind) << " individuals." << endl;
    
	cout << n_flip<<" SNPs need flip LD matrix between genotypes and reference panel."<< endl;
	//cout << n_bad <<" SNPs removed due to mismatch of allels between genotypes and reference panel."<< endl;
	cout <<"Read in "<<  idx_snp <<" common SNPs among reference and genotypes."<<  endl;
	dat->n_snp = idx_snp;
    free(hap_lanc);
}

void prod_geno(Dat *dat) {
    cout << "Precalculating genotype related products." << endl;

    dat->geno1_sq = (double *) calloc(dat->n_snp, sizeof(double));
    dat->geno2_sq = (double *) calloc(dat->n_snp, sizeof(double));
    dat->geno12_prod = (double *) calloc(dat->n_snp, sizeof(double));

    if (!dat->geno1_sq || !dat->geno2_sq || !dat->geno12_prod) {
	cout << "Error: memory allocation failed for genotype related products." << endl;
	exit(EXIT_FAILURE);
    }

    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    dat->geno1_sq[i] += dat->geno1[i][j]*dat->geno1[i][j];	    
	}
    }

    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    dat->geno2_sq[i] += dat->geno2[i][j]*dat->geno2[i][j];
	}
    }

    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    dat->geno12_prod[i] += dat->geno1[i][j]*dat->geno2[i][j];
	}
    }

    cout << "Finished calculating genotype related products." << endl;
}
void update_geno(unordered_map<size_t, size_t> &geno_idx,Dat *dat){
	size_t n_snp = geno_idx.size();
	std::vector<std::string> chr(n_snp);
    std::vector<std::string> id(n_snp);
    std::vector<std::string> pos(n_snp);
    std::vector<std::string> ref(n_snp);
    std::vector<std::string> alt(n_snp);
	std::vector<std::string> rsid(n_snp);
    double **geno1 = (double **) malloc(n_snp*sizeof(double*));
    double **geno2 = (double **) malloc(n_snp*sizeof(double*));
	
	cout<<"Update genotypes according to LD matrix."<< endl;
	
    size_t k = 0,index =0;
	
    for (size_t i=0; i<dat->n_snp; i++) {
		auto idx = geno_idx.find(i);
    	if (idx != geno_idx.end()) {
		index = idx->second;
		//cout << "i "<<i <<" index "<< index<< endl;
        geno1[index]=dat->geno1[i];
        geno2[index]=dat->geno2[i];
		chr[index]=dat->chr[i];
		id[index] = dat->id[i];       // Copy id[i] to id[index]
		pos[index] = dat->pos[i];     // Copy pos[i] to pos[index]
		ref[index] = dat->ref[i];     // Copy ref[i] to ref[index]
		alt[index] = dat->alt[i];     // Copy alt[i] to alt[index]
		rsid[index] = dat->rsid[i];
	    k++;
        } else {
			if (dat->geno1[i] != nullptr) {
            free(dat->geno1[i]);
        }
        if (dat->geno2[i] != nullptr) {
            free(dat->geno2[i]);
        }
        
    }
    }
	//cout <<"index " <<index << " and k " << k  <<endl; 
    free(dat->geno1);
    free(dat->geno2);
	size_t n_bad = dat->n_snp - n_snp;
    dat->geno1 = geno1;
    dat->geno2 = geno2;
    dat->n_snp = n_snp;
    dat->chr = chr;
    dat->id = id;
    dat->pos = pos;
    dat->ref = ref;
    dat->alt = alt;
	dat->rsid = rsid;
    cout << "Excluded " + std::to_string(n_bad) << " SNPs due to mismatch of summary stats and genoinfo." << endl;
	cout << std::to_string(dat->n_snp )<< " SNPs remain."<< endl;
}


//note the order and delete the bad snps coordinfo later
void check_maf(Dat *dat, double maf) {
    size_t n_bad = 0;
    std::vector<size_t> ok_idx;
    std::vector<std::string> chr;
    std::vector<std::string> id;
    std::vector<std::string> pos;
    std::vector<std::string> ref;
    std::vector<std::string> alt;
	std::vector<std::string> rsid;

    dat->maf1 = (double *) malloc(dat->n_snp*sizeof(double));
    dat->maf2 = (double *) malloc(dat->n_snp*sizeof(double));
    
    for (size_t i=0; i<dat->n_snp; i++) {
	double mean1 = 0, mean2 = 0;
	for (size_t j=0; j<dat->n_ind; j++) {
	    mean1 += dat->geno1[i][j]; 
	    mean2 += dat->geno2[i][j]; 	
	}
	mean1 /= (dat->n_anc1[i]); 	
	mean2 /= (dat->n_anc2[i]);
	if (mean1 > .5) {
	    mean1 = 1 - mean1;
	}
	if (mean2 > .5) {
	    mean2 = 1 - mean2;
	}
	dat->maf1[i] = mean1;
	dat->maf2[i] = mean2;

	if (mean1 < maf || mean1 > 1-maf || mean2 < maf || mean2 > 1-maf) {
	    n_bad++;
	    continue;
	}
	ok_idx.push_back(i);
    }
    
    size_t n_snp = dat->n_snp - n_bad;

    double **geno1 = (double **) malloc(n_snp*sizeof(double*));
    double **geno2 = (double **) malloc(n_snp*sizeof(double*));

    size_t k = 0;
    for (size_t i=0; i<dat->n_snp; i++) {
	if (i == ok_idx[k]) {
	    geno1[k] = dat->geno1[i];
	    geno2[k] = dat->geno2[i];
	    chr.push_back(dat->chr[i]);
	    id.push_back(dat->id[i]);
	    pos.push_back(dat->pos[i]);
	    ref.push_back(dat->ref[i]);
	    alt.push_back(dat->alt[i]);
		rsid.push_back(dat->rsid[i]);
	    k++;
	}
	else {

	    free(dat->geno1[i]);
	    free(dat->geno2[i]);
	}
    }
    free(dat->geno1);
    free(dat->geno2);
    dat->geno1 = geno1;
    dat->geno2 = geno2;
    dat->n_snp = n_snp;
    dat->chr = chr;
    dat->id = id;
    dat->pos = pos;
    dat->ref = ref;
    dat->alt = alt;
	dat->rsid = rsid;
    cout << "Excluded " + std::to_string(n_bad) << " SNPs due to MAF filter." << endl;
}

// Function to save geno1 and geno2 to text files, and pos, chr, id, ref, alt to another file
void save_geno(Dat *dat) {
    std::string chr = dat->chr[0];
    std::string geno1_filename = "geno1_" + chr + ".txt";
    std::string geno2_filename = "geno2_" + chr + ".txt";

    // Create output files
    std::ofstream geno1_file(geno1_filename);
    std::ofstream geno2_file(geno2_filename);
    
    
    if (!geno1_file.is_open() || !geno2_file.is_open()) {
        std::cerr << "Error opening file for writing.\n";
        return;
    }

    // Write geno1 matrix to file
    for (size_t i = 0; i < dat->n_snp; ++i) {
        for (size_t j = 0; j < dat->n_ind; ++j) {
            geno1_file << dat->geno1[i][j];
            if (j < dat->n_ind - 1) geno1_file << "\t";
        }
        geno1_file << "\n";
    }
    geno1_file.close();
    cout << "Write geno1." << endl;
    // Write geno2 matrix to file
    for (size_t i = 0; i < dat->n_snp; ++i) {
        for (size_t j = 0; j < dat->n_ind; ++j) {
            geno2_file << dat->geno2[i][j];
            if (j < dat->n_ind - 1) geno2_file << "\t";
        }
        geno2_file << "\n";
    }
    geno2_file.close();
    cout << "Write geno2." << endl;
    
    std::string boundary_filename = "boundary_" + chr + ".csv";
	std::ofstream boundary_out(boundary_filename);
    for (const auto& pair : dat->boundary) {
        boundary_out << pair.first << "," << pair.second << "\n";
    }
    boundary_out.close();
    
    // Write ref_ld_mat to a binary file
    std::string ldmat_filename = "ldmat_" + chr + ".bin";
    std::ofstream ldmat_out(ldmat_filename, std::ios::binary);
    for (const auto& mat : dat->ref_ld_mat) {
        size_t rows = mat->size1;
        size_t cols = mat->size2;
        ldmat_out.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        ldmat_out.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        
        for (size_t i = 0; i < rows; ++i) {
            ldmat_out.write(reinterpret_cast<const char*>(gsl_matrix_ptr(mat, i, 0)), cols * sizeof(double));
        }
    }
    ldmat_out.close();
	// Open file for pos, chr, id, ref, alt
    std::string sum_filename = "sumst_" + chr + ".txt";
    std::ofstream info_file(sum_filename);
    // Write header
    info_file << "chr\tid\tpos\tref\talt\trsid\tbeta\n";
	if (!info_file.is_open()) {
        std::cerr << "Error opening file for writing.\n";
        return;
    }
    // Write pos, chr, id, ref, alt data to file
    for (size_t i = 0; i < dat->pos.size(); ++i) {
        info_file << dat->chr[i] << "\t"
                  << dat->id[i] << "\t"
                  << dat->pos[i] << "\t" 
                  << dat->ref[i] << "\t"
                  << dat->alt[i] << "\t"
				  << dat->rsid[i] << "\t"
				  << dat->beta_mrg[i]<<"\n";
    }
    info_file.close();
    cout << "Write genotype related information." << endl;
}

double sign(double x) {
    if (x > 0) return 1.0;
    else if (x < 0) return -1.0;
    else return 0;
}

void parse_ref(const string &ref_path, unordered_map<string, CoordInfo*> &ref_dict, vector<pair<size_t, size_t>> &boundary, vector<string> &SNP) {
    ifstream infile(ref_path.c_str());
    string id_st, A1, A2, line;

    if (!infile) {
	throw std::runtime_error("Error: cannot open ref "
		"snpInfo file: " + ref_path);
    }

    int n = 0, section = 1;
    size_t left, right;
   
    // skip the header
    getline(infile, line);
    while (getline(infile, line)) {
	if (line == "") {
	    section++;
	    cout << "Readed " << n << " LD blocks." << endl;
	    n = 0;
	    getline(infile, line); // skip the header
	    continue;
	}

	if (section == 1) {
	    std::istringstream my_stream(line);
	    my_stream >> left >> right;
	    boundary.push_back(std::make_pair(left, right));
	    n++;
	}
	else {
	    std::istringstream my_stream(line);
	    my_stream >> id_st >> A1 >> A2;
	    SNP.push_back(id_st);
	    CoordInfo *ref_info = new CoordInfo;
	    ref_info->A1 = A1; ref_info->A2 = A2;
	    ref_info->include_ref = false; 
	    ref_info->include_ss = false;
        ref_info->include_geno = false;
	    ref_info->beta = 0;
	    if (!ref_dict.insert(pair<string, \
			CoordInfo*>(id_st, ref_info)).second) {
		throw std::runtime_error("Error: duplicate SNP found "
			"in ref snpInfo file: " + id_st);
	    }
	    n++;
	}
    }

    cout << "Readed " << n << " SNPs from reference panel" << endl;

    infile.close();
}

void parse_valid(const string &valid_path, unordered_map<string, CoordInfo*> &ref_dict) {
    ifstream infile(valid_path.c_str());
    if (!infile) {
	throw std::runtime_error("Error: cannot open "
		"ref snpInfo file: " + valid_path);
    }

    string id, A1, A2, header;
    float genPos; unsigned chr, phyPos;
    int n = 0;
    unordered_map<string, CoordInfo*>::iterator idx;

    while (infile >> chr >> id >> genPos >> phyPos >> A1 >> A2) {
	idx = ref_dict.find(id);
	if (idx != ref_dict.end()) {
	    idx->second->include_ref = true;
	    n++;
	}
    }
    cout << n << " common SNPs between reference "
	"and validation datasets." << endl;
    infile.close();
}

void parse_ss(const string &ss_path, unordered_map<string, CoordInfo*> &ref_dict, unsigned sz, int opt_llk) {
    ifstream infile(ss_path.c_str());
    if (!infile) {
	throw std::runtime_error("Error: cannot open "
		"summary statistics: " + ss_path);
    }

    string id, A1, A2, line;
    std::stringstream ss;
    double beta = 0, pval = 0, N = 1, Z = 0;
    int n = 0, array = 0;
    unordered_map<string, CoordInfo*>::iterator idx;

    vector<string> tokens;
    size_t SNP_idx, A1_idx, A2_idx;
    int beta_idx = -1, pval_idx = -1, array_idx = -1; 
    int sz_idx = -1, Z_idx = -1;

    int n_flip = 0, n_bad = 0, nline = 0;
    cout << "Reading GWAS summary statistics from: " + ss_path + "." << endl;
    
    while (getline(infile, line, '\n')) {
	ss.str(line);
	while (ss >> line) {
	    tokens.push_back(line);
	}
	if (nline == 0) {
	    // find corresponding fields in the header
	    vector<string>::iterator token_idx;

	    token_idx = find(tokens.begin(), tokens.end(), "SNP");
	    if (token_idx != tokens.end()) {
		SNP_idx = token_idx - tokens.begin();
	    }
	    else {
		throw std::runtime_error("Error: cannot find SNP column.");
	    }

	    token_idx = find(tokens.begin(), tokens.end(), "A1");
	    if (token_idx != tokens.end()) {
		A1_idx = token_idx - tokens.begin();
	    }
	    else {
		throw std::runtime_error("Error: cannot find A1 column.");
	    }

	    token_idx = find(tokens.begin(), tokens.end(), "A2");
	    if (token_idx != tokens.end()) {
		A2_idx = token_idx - tokens.begin();
	    }
	    else {
		throw std::runtime_error("Error: cannot find A2 column.");
	    }
	    
	    token_idx = find(tokens.begin(), tokens.end(), "BETA");
	    if (token_idx != tokens.end()) {
		beta_idx = token_idx - tokens.begin();
	    }

	    token_idx = find(tokens.begin(), tokens.end(), "P");
	    if (token_idx != tokens.end()) {
		pval_idx = token_idx - tokens.begin();
	    }
	    
	    token_idx = find(tokens.begin(), tokens.end(), "N");
	    if (token_idx != tokens.end()) {
		sz_idx = token_idx - tokens.begin();
	    }

	    token_idx = find(tokens.begin(), tokens.end(), "Z");
	    if (token_idx != tokens.end()) {
		Z_idx = token_idx - tokens.begin();
	    }

	    if (opt_llk == 2) {
		token_idx = find(tokens.begin(), tokens.end(), "ARRAY");
		if (token_idx != tokens.end()) {
		    array_idx = token_idx - tokens.begin();
		}
		else {
		    throw std::runtime_error("Error: cannot find ARRAY column.");
		}
	    }
	}
	else {
	    // parse fields
	    id = tokens[SNP_idx]; 
	    A1 = tokens[A1_idx]; A2  = tokens[A2_idx];
	    if (Z_idx > 0) {
		Z = std::stod(tokens[Z_idx]);
	    }
	    else {
		if (pval_idx < 0 || beta_idx < 0) {
		    throw std::runtime_error("Error: cannot find BETA or P column.");
		}
		pval = std::stod(tokens[pval_idx]); 
		beta = std::stod(tokens[beta_idx]);
	    }
	    if (sz_idx > 0) {
		N = std::stod(tokens[sz_idx]); 
		sz = N;
	    }
	    if (opt_llk == 2) {
		array = std::stoi(tokens[array_idx]);
	    }

	    // coordination
	    idx = ref_dict.find(id);
	    if (pval <= 1e-308) {
		pval = 1e-308;
	    }
	    if (idx != ref_dict.end() && idx->second->include_ref && idx->second->include_geno) {
		if (A1 == idx->second->A1 && A2 == idx->second->A2) {
		    idx->second->include_ss = true;
			idx->second->rsid = id;
		    if (Z_idx > 0) {
			idx->second->beta = Z/sqrt(sz);
		    }
		    else {
			idx->second->beta = 1.0*sign(beta)* \
			fabs(gsl_cdf_ugaussian_Pinv(pval/2.0))/sqrt(sz);
		    }
		    // Added for evaluating likelihood involving Ns
		    if (opt_llk == 2) {
			idx->second->array = array;
			idx->second->sz = sz;
		    }
		    n++;
		}
		else if (A1 == idx->second->A2 && A2 == idx->second->A1) {
		    idx->second->include_ss = true;
			idx->second->rsid = id;
		    if (Z_idx > 0) {
			idx->second->beta = -1.0*Z/sqrt(sz);
		    }
		    else {
			idx->second->beta = -1.0*sign(beta)* \
			fabs(gsl_cdf_ugaussian_Pinv(pval/2.0))/sqrt(sz);
		    }
		    // Added for evaluating likelihood involving Ns
		    if (opt_llk == 2) {
			idx->second->array = array;
			idx->second->sz = sz;
		    }
		    n++;
		    n_flip++;
		}
		else {
		    n_bad++;
		}
	    }
	}
	nline++;
	tokens.clear();
	ss.clear();
    }

    cout << n_flip << " SNPs have flipped alleles between summary statistics and " 
	<< "reference panel." << endl;
    //cout << n_bad << " SNPs removed due to mismatch of allels between summary statistics and reference panel." << endl;
    cout << n << " common SNPs among reference, validation, genotypes"
	" and gwas summary statistics." << endl;
    infile.close();
}


void parse_ld_mat(const string &ldmat_path, unordered_map<string, CoordInfo*> &ref_dict, const vector<pair<size_t, size_t>> &boundary, const vector<string> &SNP, Dat *dat, int opt_llk) {
    
    FILE *fp;
    fp = fopen(ldmat_path.c_str(), "rb");
	size_t n_flip = 0;

    if (!fp) {
        throw std::runtime_error("Error: cannot open LD matrix file: " + ldmat_path);
    }

    unordered_map<string, CoordInfo*>::iterator idx;
    vector<size_t> snp_idx;
	unordered_map<size_t,size_t> geno_idx;
	vector<bool> flip_yes;
    size_t left = 0, right = 0,n=0;
    for (size_t i = 0; i < boundary.size(); i++) {
        snp_idx.clear();
		flip_yes.clear();
        for (size_t j = boundary[i].first; j < boundary[i].second; j++) {
            idx = ref_dict.find(SNP[j - boundary[0].first]);
            if (idx->second->include_ss && idx->second->include_geno) {
                dat->id.push_back(SNP[j - boundary[0].first]);
				dat->rsid.push_back(idx->second->rsid);
                dat->A1.push_back(idx->second->A1);
                dat->A2.push_back(idx->second->A2);
                dat->beta_mrg.push_back(idx->second->beta);
                snp_idx.push_back(j - boundary[i].first);
				
				if (!geno_idx.insert(pair<size_t, size_t>(idx->second->geno_idx,n)).second) {
					throw std::runtime_error("Error: duplicate SNP found "
					"in geno snp: " + idx->second->geno_idx);
	    		}
				flip_yes.push_back(idx->second->flip);
                // Added for evaluating llk involving Ns
                if (opt_llk == 2) {
                    dat->sz.push_back(idx->second->sz);
                    dat->array.push_back(idx->second->array);
                }
                right++;n++;
            }
        }
        gsl_matrix *tmp_mat = gsl_matrix_alloc(boundary[i].second - boundary[i].first, \
                        boundary[i].second - boundary[i].first);
        gsl_matrix_fread(fp, tmp_mat);

        if (left == right) {
            gsl_matrix_free(tmp_mat);
            continue; 
        }
        dat->boundary.push_back(std::make_pair(left, right));
        gsl_matrix *tmp_mat_sub = gsl_matrix_alloc(right - left, right - left);

        // copy rows from original matrix to second matrix with correct SNPs 
        for (size_t j = 0; j < snp_idx.size(); j++) {
            for (size_t k = 0; k < snp_idx.size(); k++) {
                double tmp = gsl_matrix_get(tmp_mat, snp_idx[j], snp_idx[k]);
                gsl_matrix_set(tmp_mat_sub, j, k, tmp);
            }
        }
		
		//apply row flips
		for (size_t j = 0; j < snp_idx.size(); j++) {
			if (flip_yes[j]) {
				for (size_t k = 0; k < snp_idx.size(); k++) {
					double tmp2 = gsl_matrix_get(tmp_mat_sub, j, k);
					gsl_matrix_set(tmp_mat_sub, j, k, -tmp2);
				}
			}
		}
		// Apply column flips
		for (size_t k = 0; k < snp_idx.size(); k++) {
			if (flip_yes[k]) {
				for (size_t j = 0; j < snp_idx.size(); j++) {
					double tmp3 = gsl_matrix_get(tmp_mat_sub, j, k);
					gsl_matrix_set(tmp_mat_sub, j, k, -tmp3);
				}
				n_flip++;
			}
		}
		dat->ref_ld_mat.push_back(tmp_mat_sub);
        left = right;
        gsl_matrix_free(tmp_mat);
        
    }
	cout << n <<" SNPs in LD matrix with " <<n_flip << " flipped."<<endl;
	//copy rows from original matrix to second matrix with correct geno_id
	update_geno(geno_idx,dat);

    fclose(fp);
	
}


void coord(const string &vcf_path, const string &msp_path, const string &bim_path,const string &ref_path, const string &ss_path, const string &valid_path, const string &ldmat_path, Dat *dat, unsigned sz, int opt_llk) {
    unordered_map<string, CoordInfo*> ref_dict;
    unordered_map<string, string> bim_dict;
    vector<pair<size_t, size_t>> boundary;
    vector<string> SNP;
    unordered_map<string, CoordInfo*>::iterator it;

    parse_ref(ref_path, ref_dict, boundary, SNP);
    parse_bim(dat, bim_path, bim_dict);
    for (it=ref_dict.begin(); it != ref_dict.end(); it++) {
	    it->second->flip = false;
	}
    
    if (!valid_path.empty()) {
	parse_valid(valid_path, ref_dict);
    }
    else {
	for (it=ref_dict.begin(); it != ref_dict.end(); it++) {
	    it->second->include_ref = true;
	}
    }
	
    
    read_lanc(vcf_path,  msp_path, dat, bim_dict, ref_dict);
    //cout << "Finished coordinating genotype with bim and reference." << endl;
	bim_dict.clear();
    parse_ss(ss_path, ref_dict, sz, opt_llk);  
	/* size_t n=0;
	for (it=ref_dict.begin(); it != ref_dict.end(); it++) {
		if(it->second->include_geno){
		cout << n<<" "<<it->second->geno_idx<<endl;
		n++;                                             
		}
	}*/
	check_maf(dat,0.0);
    parse_ld_mat(ldmat_path, ref_dict, boundary, SNP, dat, opt_llk);
    for (it=ref_dict.begin(); it != ref_dict.end(); it++) {
	 delete(it->second);
    }
	
    cout << "Coordinate the summary statistics and reference and geno types." << endl;
}


/*
int main() {
    Dat dat;
    coord("test.snpInfo", "/ysm-gpfs/pi/zhao-data//gz222/height/summ_stats/SDPR.txt", \
	    "/ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3.bim", \
	    "test.dat", dat, 252230);
    for (size_t i=0; i<dat->ref_ld_mat.size(); i++) {
	gsl_matrix_free(dat->ref_ld_mat[i]);
    }
    coord_genos(dat, "/gpfs/gibbs/pi/zhao/gz222/SDPR_admix/Real/genotype/UKB/Ukb_imp_v2.bim")
    return 0;
}*/  

