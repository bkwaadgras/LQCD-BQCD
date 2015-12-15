#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#define TEST_MAIN
//#define COMPILE_FOR_TMLQCD
#include "su3.h"
#include "global.h"


/************* Calls to DFE ********************/
void init_dfe(void);
void transfer_spinors_to_dfe (bqcd_spinor *t_in, int address);
void transfer_fix_spinors_to_dfe (bqcd_fix_spinor *t_in, int address);
void transfer_gauges_to_dfe (su3 *t_in, int address);
void transfer_fix_gauges_to_dfe (fix_su3 *t_in, int address);
void transfer_clovers_to_dfe (bqcd_clover *t_in, int address);
void transfer_fix_clovers_to_dfe (bqcd_fix_clover *t_in, int address);
/*void apply_dirac (int in_address, int out_address, int gauge_address, int p_address,
		int clover_address, int ieo, int doSub, int isign, int apply_dslash, int apply_clover,
		uint64_t *out_exp64, uint64_t *out_sq64);*/
void apply_dirac (dirac_run_parameters *run_params);
void apply_misc_ops (dirac_run_parameters *run_params);
void transfer_spinors_to_cpu (bqcd_spinor *out, int address);
void transfer_fix_spinors_to_cpu (bqcd_fix_spinor *out, int address);
void unload_dfe(void);

/************* Verifying Functions ***************/
int verify_results (bqcd_spinor* dfe_out, bqcd_spinor *expected_out, int V);
int verify_results_extended (bqcd_spinor* dfe_out, bqcd_spinor *expected_out, int V);
void output_spinor_analysis (char *filename, bqcd_spinor* dfe_out,
		bqcd_spinor *expected_out, int V);
int AreNotSameComplex(complex float a, float complex b);
int compare_spinor (bqcd_spinor *a, bqcd_spinor *b);
float calc_normsq (bqcd_spinor* t_in, int V);

/************ Utility Functions *********************/
void create_random_input(spinor* s, su3* u);
void create_random_spinor(spinor * s);
void create_random_su3vector(su3_vector *v);
void create_random_su3(su3 *s);
void create_random_complex(float complex *a);
void print_spinors (spinor* s);
void print_fix_spinors (bqcd_fix_spinor* s);
void print_gauges (su3* s);

/****************** Functions for reading data from and writing data to files ***************/
void conv_bqcd_spinor_to_bqcd_fix_spinor(bqcd_fix_spinor *out, bqcd_spinor *t_in,
														int exponent, int frac_bits);
void conv_bqcd_fix_spinor_to_bqcd_spinor(bqcd_spinor *out, bqcd_fix_spinor *t_in,
														int exp_shift, int frac_bits);
void conv_su3_to_fix_su3(fix_su3 *out, su3 *t_in, int exponent, int frac_bits);
void conv_bqcd_clover_to_bqcd_fix_clover(bqcd_fix_clover *out, bqcd_clover *t_in,
				int exp_shift_diag, int exp_shift_offdiag, int frac_bits);
void conv_enc_fix_spinor_to_bqcd_spinor(bqcd_spinor *out, bqcd_fix_spinor *t_in,
																		int frac_bits);
void conv_bqcd_spinor_to_enc_fix_spinor(bqcd_fix_spinor *out, bqcd_spinor *in,
															int frac_bits);
void conv_fix_su3_to_su3(su3 *out, fix_su3 *t_in, int exponent, int frac_bits);
void read_spinor(char * filename, bqcd_spinor *out, int *maximum_exponent);
void read_gauge(char * filename, su3 *s);
void read_clover(char * filename, bqcd_clover *result,
							int *max_exp_diagonal, int *max_exp_offdiagonal);
void write_spinor(char *filename, bqcd_spinor *to_write, char *format);
void write_gauge(char *filename, su3 *s);
void write_clover(char * filename, bqcd_clover *to_write);
void conv_float_to_fix24(fix24 *out, float t_in,  int exponent, int frac_bits);
float conv_fix24_to_float(fix24 t_in, int exponent, int frac_bits);

/****************** Functions for reordering data in memory ***************/
//inline size_t index_from_coords(int x, int y, int z, int t);
void reorganize_gauge(su3 const * const t_in, su3 * const out, int ieo);
void new_reorganize_gauge(su3 const * const t_in, su3 * const out, int ieo);
void collate_gauges(su3* u0, su3* u1, su3 *ures, int ieo, int *maximum_exponent);

static max_file_t *maxfile;
static max_engine_t *engine;
static int burstsPerGaugeTimeSlice, burstsPerSpinorTimeSlice, burstsPerCloverTimeSlice;
static double beta_s, beta_t_b, beta_t_f, mass;
//static float bqcd_b;
//static float cg_beta;
//static float rtrold;
static double elapsedTime;

int min(int a, int b) { return a > b ? b : a; }

int main(void)
{
	LT = LQCD_T;
	LX = LQCD_LX;
	LY = LQCD_LY;
	LZ  = LQCD_LZ;

	VOLUME = LZ * LY * LX * LT;

	NUM_PIPES = LQCD_numPipes;
	LOOP_OFFSET = LQCD_loopOffset;

	int n_frac_bits=-LQCD_realStoreOffset;
	int n_frac_bits_outsq=32;

	/*beta_s   = -0.5;
	beta_t_f = 0.3;
	beta_t_b = 0.3;
	mass     = 0.1;*/
	beta_s = 1.;
	beta_t_f = 1.;
	beta_t_b = 1.;
	mass = 1.;
	float cg_beta = 0.123455;

	float bqcd_b = -1.44720900000000017E-002;

	float rtrold = 337.361268;

	setbuf(stdout, NULL);

	printf("Allocating memory for data ...");

	bqcd_spinor *in_x = malloc(14 * VOLUME/2 * sizeof(bqcd_spinor));
	bqcd_spinor *in_p = &in_x[VOLUME / 2];
	bqcd_spinor *in_r = &in_p[VOLUME / 2];
	bqcd_spinor *in_u = &in_r[VOLUME / 2];
	bqcd_spinor *out_x_dfe = &in_u[VOLUME / 2];
	bqcd_spinor *out_r_dfe = &out_x_dfe[VOLUME / 2];
	bqcd_spinor *out_s_dfe = &out_r_dfe[VOLUME / 2];
	bqcd_spinor *out_p_dfe = &out_s_dfe[VOLUME / 2];
	bqcd_spinor *out_u_dfe = &out_p_dfe[VOLUME / 2];
	bqcd_spinor *out_x_expected = &out_u_dfe[VOLUME / 2];
	bqcd_spinor *out_r_expected = &out_x_expected[VOLUME / 2];
	bqcd_spinor *out_s_expected = &out_r_expected[VOLUME / 2];
	bqcd_spinor *out_p_expected = &out_s_expected[VOLUME / 2];
	bqcd_spinor *out_u_expected = &out_p_expected[VOLUME / 2];

	su3 *gauge0 = malloc (4 * VOLUME/2 * 8 * sizeof(su3));
	su3 *gauge1 = &gauge0[VOLUME/2 * 8];
	su3 *gauge_merged0 = &gauge1[VOLUME/2 * 8];
	su3 *gauge_merged1 = &gauge_merged0[VOLUME/2 * 8];

	bqcd_clover *clover0 = malloc(2 * VOLUME/2 * sizeof(bqcd_clover));
	bqcd_clover *clover1 = &clover0[VOLUME/2];

	bqcd_fix_spinor *in_x_fix = malloc(9 * VOLUME/2 * sizeof(bqcd_fix_spinor));
	bqcd_fix_spinor *in_p_fix = &in_x_fix[VOLUME/2];
	bqcd_fix_spinor *in_r_fix = &in_p_fix[VOLUME/2];
	bqcd_fix_spinor *in_u_fix = &in_r_fix[VOLUME/2];
	bqcd_fix_spinor *out_x_fix_dfe = &in_u_fix[VOLUME/2];
	bqcd_fix_spinor *out_r_fix_dfe = &out_x_fix_dfe[VOLUME/2];
	bqcd_fix_spinor *out_s_fix_dfe = &out_r_fix_dfe[VOLUME/2];
	bqcd_fix_spinor *out_p_fix_dfe = &out_s_fix_dfe[VOLUME/2];
	bqcd_fix_spinor *out_u_fix_dfe = &out_p_fix_dfe[VOLUME/2];

	fix_su3 *gauge_merged0_fix = malloc (2 * VOLUME/2 * 8 * sizeof(fix_su3));
	fix_su3 *gauge_merged1_fix = &gauge_merged0_fix[VOLUME/2 * 8];
	bqcd_fix_clover *clover0_fix = malloc(2 * VOLUME/2 * sizeof(bqcd_fix_clover));
	bqcd_fix_clover *clover1_fix = &clover0_fix[VOLUME/2];

	bqcd_spinor *out_tmp_expected = malloc(VOLUME/2 * sizeof(bqcd_spinor));
	bqcd_spinor *out_tmp_dfe = malloc(VOLUME/2 * sizeof(bqcd_spinor));
	bqcd_fix_spinor *out_tmp_fix_dfe = malloc(VOLUME/2 * sizeof(bqcd_spinor));

	printf("Done!\n");

	printf("Reading spinor and gauge inputs ...");

	int maximum_exponent_x;
	int maximum_exponent_p;
	int maximum_exponent_r;
	int maximum_exponent_u;
	int maximum_exponent_out_x;
	int maximum_exponent_out_r;
	int maximum_exponent_out_s;
	int maximum_exponent_out_p;
	int maximum_exponent_out_u;

	read_spinor("x_prev_spinor.txt", in_x, &maximum_exponent_x);
	read_spinor("p_prev_spinor.txt", in_p, &maximum_exponent_p);
	read_spinor("r_prev_spinor.txt", in_r, &maximum_exponent_r);
	read_spinor("u_prev_spinor.txt", in_u, &maximum_exponent_u);
	read_spinor("x_spinor.txt", out_x_expected, &maximum_exponent_out_x);
	read_spinor("r_spinor.txt", out_r_expected, &maximum_exponent_out_r);
	read_spinor("s_spinor.txt", out_s_expected, &maximum_exponent_out_s);
	read_spinor("p_spinor.txt", out_p_expected, &maximum_exponent_out_p);
	read_spinor("aap_spinor.txt", out_u_expected, &maximum_exponent_out_u);

	read_gauge("cg_mtil_in_gauge0.txt", gauge0);
	read_gauge("cg_mtil_in_gauge1.txt", gauge1);

	int max_exp_clover0_offdiag;
	int max_exp_clover0_diag;
	int max_exp_clover1_offdiag;
	int max_exp_clover1_diag;

	read_clover("cg_mtil_clover0_f.txt", clover0,
			&max_exp_clover0_diag, &max_exp_clover0_offdiag);
	read_clover("cg_mtil_clover1_f.txt", clover1,
			&max_exp_clover1_diag, &max_exp_clover1_offdiag);

	printf("Max exp clover 0 diagonal = %d\n", max_exp_clover0_diag);
	printf("Max exp clover 0 offdiag = %d\n", max_exp_clover0_offdiag);
	printf("Max exp clover 1 diagonal = %d\n", max_exp_clover1_diag);
	printf("Max exp clover 1 offdiag = %d\n", max_exp_clover1_offdiag);
	if(max_exp_clover0_offdiag > 0) { max_exp_clover0_offdiag=0; }
	if(max_exp_clover0_diag > 0) { max_exp_clover0_diag=0; }
	if(max_exp_clover1_offdiag > 0) { max_exp_clover0_offdiag=0; }
	if(max_exp_clover1_diag > 0) { max_exp_clover0_diag=0; }

	//conv_bqcd_spinor_to_bqcd_fix_spinor(in_x_fix, in_x, -maximum_exponent_x, n_frac_bits);
	conv_bqcd_spinor_to_enc_fix_spinor(in_x_fix, in_x, n_frac_bits);
	//conv_bqcd_spinor_to_bqcd_fix_spinor(in_p_fix, in_p, -maximum_exponent_p, n_frac_bits);
	conv_bqcd_spinor_to_enc_fix_spinor(in_p_fix, in_p, n_frac_bits);
	conv_bqcd_spinor_to_enc_fix_spinor(in_r_fix, in_r, n_frac_bits);
	conv_bqcd_spinor_to_enc_fix_spinor(in_u_fix, in_u, n_frac_bits);
//	conv_bqcd_spinor_to_bqcd_fix_spinor(in_r_fix, in_r, -maximum_exponent_r, n_frac_bits);
//	conv_bqcd_spinor_to_bqcd_fix_spinor(in_u_fix, in_u, -maximum_exponent_u, n_frac_bits);

	printf("Done!\n");
	printf("Data reordering and adding necessary halos ...\n");

	int max_exp_gauge_merged0;
	int max_exp_gauge_merged1;

	collate_gauges(gauge0, gauge1, gauge_merged0, 0, &max_exp_gauge_merged0);
	collate_gauges(gauge0, gauge1, gauge_merged1, 1, &max_exp_gauge_merged1);
	if(max_exp_gauge_merged0 > 0) { max_exp_gauge_merged0=0; }
	if(max_exp_gauge_merged1 > 0) { max_exp_gauge_merged1=0; }

	write_gauge("debug_gauge_merged0.txt", gauge_merged0);
	write_gauge("debug_gauge_merged1.txt", gauge_merged1);

	printf("merged gauge 0 max exp = %d\n", max_exp_gauge_merged0);
	printf("merged gauge 1 max exp = %d\n", max_exp_gauge_merged1);
	conv_su3_to_fix_su3(gauge_merged0_fix, gauge_merged0, -max_exp_gauge_merged0, n_frac_bits);
	conv_su3_to_fix_su3(gauge_merged1_fix, gauge_merged1, -max_exp_gauge_merged1, n_frac_bits);

	conv_bqcd_clover_to_bqcd_fix_clover(clover0_fix, clover0,
			-max_exp_clover0_diag, -max_exp_clover0_offdiag, n_frac_bits);
	conv_bqcd_clover_to_bqcd_fix_clover(clover1_fix, clover1,
			-max_exp_clover1_diag, -max_exp_clover1_offdiag, n_frac_bits);

	printf("Done!\n");

	init_dfe();

	int address_g0		= 0;
	int address_g1		= address_g0 		+ LT * burstsPerGaugeTimeSlice;
	int address_x		= address_g1  		+ LT * burstsPerGaugeTimeSlice;
	int address_r		= address_x 		+ LT * burstsPerSpinorTimeSlice;
	int address_p1		= address_r			+ LT * burstsPerSpinorTimeSlice;
	int address_s		= address_p1		+ LT * burstsPerSpinorTimeSlice;
	int address_p0		= address_s			+ LT * burstsPerSpinorTimeSlice;
	int address_u		= address_p0		+ LT * burstsPerSpinorTimeSlice;
	int address_clover0	= address_u			+ LT * burstsPerSpinorTimeSlice;
	int address_clover1	= address_clover0 	+ LT * burstsPerCloverTimeSlice;

	printf("Transferring Gauge to Memory  ...");
	transfer_fix_gauges_to_dfe(gauge_merged1_fix, address_g1);
	transfer_fix_gauges_to_dfe(gauge_merged0_fix, address_g0);
	printf("Done!\n");
	fflush(stdout);

	printf("Transferring Input Spinors to Memory...");
	transfer_fix_spinors_to_dfe (in_x_fix, address_x);
	transfer_fix_spinors_to_dfe (in_p_fix, address_p0);
	transfer_fix_spinors_to_dfe (in_r_fix, address_r);
	transfer_fix_spinors_to_dfe (in_u_fix, address_u);
	printf("Done!\n");
	fflush(stdout);

	printf("Transferring Input Clovers to Memory...");
	transfer_fix_clovers_to_dfe (clover0_fix, address_clover0) ;
	transfer_fix_clovers_to_dfe (clover1_fix, address_clover1) ;
	printf("Done!\n");
	fflush(stdout);

	elapsedTime = 0.0;

	printf("Calculating on DFE ...");
	fflush(stdout);
	float cg_alpha = 0.686326;
	float rtr=rtrold;

	dirac_run_parameters d_params_misc_run;
	d_params_misc_run.in_spinor0_address  = address_p0;
	d_params_misc_run.in_spinor1_address  = address_x;
	d_params_misc_run.out_spinor0_address = address_r;
	d_params_misc_run.out_spinor1_address = address_x;
	d_params_misc_run.gauge_address       = address_r;
	d_params_misc_run.clover_address      = address_u;

	dirac_run_parameters d_params_1st_run;
	d_params_1st_run.in_spinor0_address   = address_p0;
	d_params_1st_run.in_spinor1_address   = address_r;
	d_params_1st_run.out_spinor0_address  = address_s;
	d_params_1st_run.out_spinor1_address  = address_p1;
	d_params_1st_run.gauge_address        = address_g1;
	d_params_1st_run.clover_address       = address_clover1;
	d_params_1st_run.apply_axpy_pre       = 1;
	d_params_1st_run.apply_axpy_post      = 0;
	d_params_1st_run.enable_second_input  = 1;
	d_params_1st_run.apply_dslash         = 1;
	d_params_1st_run.apply_clover         = 1;
	d_params_1st_run.enable_clover_input  = 1;
	d_params_1st_run.enable_second_output = 1;
	d_params_1st_run.even_odd             = 1;
	d_params_1st_run.dagger               = 0;
	d_params_1st_run.misc_mode            = 0;
	d_params_1st_run.clover_diag_shift = -max_exp_clover1_diag;
	d_params_1st_run.clover_offdiag_shift = -max_exp_clover1_offdiag;
	d_params_1st_run.gauge_shift = -max_exp_gauge_merged1;
	dirac_run_parameters d_params_2nd_run;
	d_params_2nd_run.in_spinor0_address   = address_s;
	d_params_2nd_run.in_spinor1_address   = address_p1;
	d_params_2nd_run.out_spinor0_address  = address_s;
	d_params_2nd_run.gauge_address        = address_g0;
	d_params_2nd_run.clover_address       = address_clover0;
	d_params_2nd_run.apply_axpy_pre       = 0;
	d_params_2nd_run.apply_axpy_post      = 1;
	d_params_2nd_run.enable_second_input  = 1;
	d_params_2nd_run.apply_dslash         = 1;
	d_params_2nd_run.apply_clover         = 1;
	d_params_2nd_run.enable_clover_input  = 1;
	d_params_2nd_run.enable_second_output = 0;
	d_params_2nd_run.even_odd             = 0;
	d_params_2nd_run.dagger               = 0;
	d_params_2nd_run.misc_mode            = 0;
	d_params_2nd_run.clover_diag_shift = -max_exp_clover0_diag;
	d_params_2nd_run.clover_offdiag_shift = -max_exp_clover0_offdiag;
	d_params_2nd_run.gauge_shift = -max_exp_gauge_merged0;
	dirac_run_parameters d_params_3rd_run;
	d_params_3rd_run.in_spinor0_address   = address_s;
	d_params_3rd_run.out_spinor0_address  = address_u;
	d_params_3rd_run.gauge_address        = address_g1;
	d_params_3rd_run.clover_address       = address_clover0;
	d_params_3rd_run.apply_axpy_pre       = 0;
	d_params_3rd_run.apply_axpy_post      = 0;
	d_params_3rd_run.enable_second_input  = 0;
	d_params_3rd_run.apply_dslash         = 1;
	d_params_3rd_run.apply_clover         = 1;
	d_params_3rd_run.enable_clover_input  = 1;
	d_params_3rd_run.enable_second_output = 0;
	d_params_3rd_run.even_odd             = 1;
	d_params_3rd_run.dagger               = 1;
	d_params_3rd_run.misc_mode            = 0;
	d_params_3rd_run.clover_diag_shift = -max_exp_clover0_diag;
	d_params_3rd_run.clover_offdiag_shift = -max_exp_clover0_offdiag;
	d_params_3rd_run.gauge_shift = -max_exp_gauge_merged1;
	dirac_run_parameters d_params_4th_run;
	d_params_4th_run.in_spinor0_address   = address_u;
	d_params_4th_run.in_spinor1_address   = address_s;
	d_params_4th_run.out_spinor0_address  = address_u;
	d_params_4th_run.gauge_address        = address_g0;
	d_params_4th_run.clover_address       = address_clover1;
	d_params_4th_run.apply_axpy_pre       = 0;
	d_params_4th_run.apply_axpy_post      = 1;
	d_params_4th_run.enable_second_input  = 1;
	d_params_4th_run.apply_dslash         = 1;
	d_params_4th_run.apply_clover         = 1;
	d_params_4th_run.enable_clover_input  = 1;
	d_params_4th_run.enable_second_output = 0;
	d_params_4th_run.even_odd             = 0;
	d_params_4th_run.dagger               = 1;
	d_params_4th_run.misc_mode            = 0;
	d_params_4th_run.clover_diag_shift = -max_exp_clover1_diag;
	d_params_4th_run.clover_offdiag_shift = -max_exp_clover1_offdiag;
	d_params_4th_run.gauge_shift = -max_exp_gauge_merged0;

	uint64_t out_sq_raw;

	for(int niter=0; niter<1/*1200000*/; ++niter)
	{
		if(niter%2==0)
		{
			d_params_misc_run.in_spinor0_address  = address_p0;
			d_params_1st_run.in_spinor0_address   = address_p0;
			d_params_1st_run.out_spinor1_address  = address_p1;
			d_params_2nd_run.in_spinor1_address   = address_p1;
		}
		else
		{
			d_params_misc_run.in_spinor0_address  = address_p1;
			d_params_1st_run.in_spinor0_address   = address_p1;
			d_params_1st_run.out_spinor1_address  = address_p0;
			d_params_2nd_run.in_spinor1_address   = address_p0;
		}
		printf("iteration %d: cg_alpha = %f, cg_beta = %f, rtr = %E\n",
				niter, cg_alpha, cg_beta, rtr);

		d_params_misc_run.axpy_pre_mult        = cg_alpha;
		d_params_misc_run.axpy_post_mult       = -cg_alpha;
		apply_misc_ops(&d_params_misc_run);
		out_sq_raw = d_params_misc_run.out_sq;
		int outSqShift = d_params_misc_run.out_sq_shift;
		rtr=(float)out_sq_raw/pow(2., n_frac_bits_outsq - 2*outSqShift);


		if(sqrtf(rtr) <= 0.00000000057719421/* || 1*/)
		{
			break;
		}

		cg_beta = rtr/rtrold;
		rtrold=rtr;

		d_params_1st_run.axpy_pre_mult        = cg_beta;

		apply_dirac(&d_params_1st_run);

		if(niter<45)
		{
			transfer_fix_spinors_to_cpu (out_tmp_fix_dfe, address_s);
			conv_enc_fix_spinor_to_bqcd_spinor(out_tmp_dfe, out_tmp_fix_dfe, n_frac_bits);
			char *buffer2 = malloc(100*sizeof(char));
			sprintf(buffer2, "out_tmp_spinor_dfe_(k=%d).txt", niter);
			write_spinor(buffer2, out_tmp_dfe, "( %21.17E , %21.17E )\n");
			sprintf(buffer2, "./spinors-by-iteration/dslash-c/out_tmp_spinor_dslash-c_(k=%d).txt", niter);
			int temp2;
			read_spinor(buffer2, out_tmp_expected, &temp2);
			printf("verify tmp: "); verify_results(out_tmp_dfe, out_tmp_expected, VOLUME/2);
			sprintf(buffer2, "out_tmp_spinor_dfe_analysis_(k=%02d).txt", niter);
			output_spinor_analysis(buffer2, out_tmp_dfe, out_tmp_expected, VOLUME/2);
			free(buffer2);
		}

		d_params_2nd_run.axpy_post_mult       = bqcd_b;
		apply_dirac(&d_params_2nd_run);
		out_sq_raw=d_params_2nd_run.out_sq;

		outSqShift = d_params_2nd_run.out_sq_shift;
		float sts=(float)out_sq_raw/pow(2., n_frac_bits_outsq - 2*outSqShift);

		cg_alpha = rtrold/sts;

		apply_dirac(&d_params_3rd_run);

		if(niter<45)
		{
			transfer_fix_spinors_to_cpu (out_tmp_fix_dfe, address_s);
			conv_enc_fix_spinor_to_bqcd_spinor(out_tmp_dfe, out_tmp_fix_dfe, n_frac_bits);
			char *buffer2 = malloc(100*sizeof(char));
			sprintf(buffer2, "out_dag_tmp_spinor_dfe_(k=%d).txt", niter);
			write_spinor(buffer2, out_tmp_dfe, "( %21.17E , %21.17E )\n");
			sprintf(buffer2, "./spinors-by-iteration/dslash-c/out_dag_tmp_spinor(k=%d).txt", niter);
			int temp2;
			read_spinor(buffer2, out_tmp_expected, &temp2);
			printf("verify tmp: "); verify_results(out_tmp_dfe, out_tmp_expected, VOLUME/2);
			sprintf(buffer2, "out_dag_tmp_spinor_dfe_analysis_(k=%02d).txt", niter);
			output_spinor_analysis(buffer2, out_tmp_dfe, out_tmp_expected, VOLUME/2);
			free(buffer2);
		}


		d_params_4th_run.axpy_post_mult       = bqcd_b;
		apply_dirac(&d_params_4th_run);

		transfer_fix_spinors_to_cpu (out_x_fix_dfe, address_x);
		transfer_fix_spinors_to_cpu (out_r_fix_dfe, address_r);
		transfer_fix_spinors_to_cpu (out_s_fix_dfe, address_s);
		if(niter%2==0)
		{
			transfer_fix_spinors_to_cpu (out_p_fix_dfe, address_p1);
		}
		else
		{
			transfer_fix_spinors_to_cpu (out_p_fix_dfe, address_p0);
		}
		transfer_fix_spinors_to_cpu (out_u_fix_dfe, address_u);
		conv_enc_fix_spinor_to_bqcd_spinor(out_x_dfe, out_x_fix_dfe, n_frac_bits);
		conv_enc_fix_spinor_to_bqcd_spinor(out_r_dfe, out_r_fix_dfe, n_frac_bits);
		conv_enc_fix_spinor_to_bqcd_spinor(out_s_dfe, out_s_fix_dfe, n_frac_bits);
		conv_enc_fix_spinor_to_bqcd_spinor(out_p_dfe, out_p_fix_dfe, n_frac_bits);
		conv_enc_fix_spinor_to_bqcd_spinor(out_u_dfe, out_u_fix_dfe, n_frac_bits);

		char *buffer = malloc(100*sizeof(char));
		int temp;
		if(niter<45)
		{

			sprintf(buffer, "./spinors-by-iteration/dslash-c/out_r_spinor_dslash-c_(k=%d).txt", niter);
			read_spinor(buffer, out_r_expected, &temp);
			buffer[36]='s';
			read_spinor(buffer, out_s_expected, &temp);
			buffer[36]='x';
			read_spinor(buffer, out_x_expected, &temp);
			buffer[36]='p';
			read_spinor(buffer, out_p_expected, &temp);
			buffer[36]='u';
			read_spinor(buffer, out_u_expected, &temp);
			printf("verify r: "); verify_results(out_r_dfe, out_r_expected, VOLUME/2);
			printf("verify s: "); verify_results(out_s_dfe, out_s_expected, VOLUME/2);
			printf("verify x: "); verify_results(out_x_dfe, out_x_expected, VOLUME/2);
			printf("verify p: "); verify_results(out_p_dfe, out_p_expected, VOLUME/2);
			printf("verify u: "); verify_results(out_u_dfe, out_u_expected, VOLUME/2);
		}

		printf("x analysis: "); verify_results_extended(out_x_dfe, out_x_expected, VOLUME/2);
		printf("r analysis: "); verify_results_extended(out_r_dfe, out_r_expected, VOLUME/2);
		printf("p analysis: "); verify_results_extended(out_p_dfe, out_p_expected, VOLUME/2);
		printf("s analysis: "); verify_results_extended(out_s_dfe, out_s_expected, VOLUME/2);
		printf("u analysis: "); verify_results_extended(out_u_dfe, out_u_expected, VOLUME/2);

		sprintf(buffer, "out_r_spinor_dfe_(k=%d).txt", niter);
		write_spinor(buffer, out_r_dfe, "( %21.17E , %21.17E )\n");
		buffer[4]='s';
		write_spinor(buffer, out_s_dfe, "( %21.17E , %21.17E )\n");
		buffer[4]='x';
		write_spinor(buffer, out_x_dfe, "( %21.17E , %21.17E )\n");
		buffer[4]='p';
		write_spinor(buffer, out_p_dfe, "( %21.17E , %21.17E )\n");
		buffer[4]='u';
		write_spinor(buffer, out_u_dfe, "( %21.17E , %21.17E )\n");

		sprintf(buffer, "out_x_spinor_dfe_analysis_(k=%02d).txt", niter);
		output_spinor_analysis(buffer, out_x_dfe, out_x_expected, VOLUME/2);
		buffer[4]='r';
		output_spinor_analysis(buffer, out_r_dfe, out_r_expected, VOLUME/2);
		buffer[4]='p';
		output_spinor_analysis(buffer, out_p_dfe, out_p_expected, VOLUME/2);
		buffer[4]='s';
		output_spinor_analysis(buffer, out_s_dfe, out_s_expected, VOLUME/2);
		buffer[4]='u';
		output_spinor_analysis(buffer, out_u_dfe, out_u_expected, VOLUME/2);

		free(buffer);


	}

	dirac_run_parameters d_params_test;
	d_params_test.in_spinor0_address   = address_p0;
	d_params_test.in_spinor1_address   = address_r;
	d_params_test.out_spinor0_address  = address_s;
	d_params_test.out_spinor1_address  = address_p1;
	d_params_test.gauge_address        = address_g1;
	d_params_test.clover_address       = address_clover1;
	d_params_test.apply_axpy_pre       = 1;
	d_params_test.apply_axpy_post      = 0;
	d_params_test.enable_second_input  = 1;
	d_params_test.apply_dslash         = 1;
	d_params_test.apply_clover         = 1;
	d_params_test.enable_clover_input  = 1;
	d_params_test.enable_second_output = 1;
	d_params_test.even_odd             = 1;
	d_params_test.dagger               = 0;
	d_params_test.misc_mode            = 0;
	d_params_test.clover_diag_shift = -max_exp_clover1_diag;
	d_params_test.clover_offdiag_shift = -max_exp_clover1_offdiag;
	d_params_test.gauge_shift = -max_exp_gauge_merged1;



	printf("DFE time elapsed: %f ms\n", elapsedTime);
//	printf("DFE Throughput: %g FLOPS\n", (1320.0*(double)VOLUME/2*4*Iterations)*1000.0/elapsedTime);

	return  0;
}
/************* Calls to DFE ********************/
void init_dfe(void) {
	maxfile = LQCD_init();
	engine = max_load(maxfile, "*");
	int burstSize = max_get_burst_size(maxfile, 0);
	burstsPerGaugeTimeSlice  = (LX*LY*LZ/2*8*sizeof(fix_su3)) / burstSize;
	burstsPerSpinorTimeSlice = (LX*LY*LZ/2*sizeof(bqcd_fix_spinor)) / burstSize;
	burstsPerCloverTimeSlice = (LX*LY*LZ/2*sizeof(bqcd_fix_clover)) / burstSize;
}

void transfer_spinors_to_dfe (bqcd_spinor *t_in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "data_in", t_in, VOLUME/2 * sizeof(spinor));
	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);
    max_route(act, "toLmemMux_fromCPU", "toLmemMux");
    max_lmem_set_interrupt_on(act, "toLmem0");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel1");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void transfer_fix_spinors_to_dfe (bqcd_fix_spinor *t_in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "data_in", t_in, VOLUME/2 * sizeof(bqcd_fix_spinor));

	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);
    max_route(act, "toLmemMux_fromCPU", "toLmemMux");
    max_lmem_set_interrupt_on(act, "toLmem0");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel1");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}


void transfer_gauges_to_dfe (su3 *t_in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "data_in", t_in, VOLUME/2 * 8 * sizeof(su3));
	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerGaugeTimeSlice/LQCD_gCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_gCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerGaugeTimeSlice/LQCD_gCmdSize);
    max_route(act, "toLmemMux_fromCPU", "toLmemMux");
	max_lmem_set_interrupt_on(act, "toLmem0");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel1");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void transfer_fix_gauges_to_dfe (fix_su3 *t_in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "data_in", t_in, VOLUME/2 * 8 * sizeof(fix_su3));
	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerGaugeTimeSlice/LQCD_gCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_gCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerGaugeTimeSlice/LQCD_gCmdSize);
    max_route(act, "toLmemMux_fromCPU", "toLmemMux");
	max_lmem_set_interrupt_on(act, "toLmem0");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel1");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void transfer_clovers_to_dfe (bqcd_clover *t_in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "data_in", t_in, VOLUME/2 * sizeof(bqcd_clover));
	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerCloverTimeSlice/LQCD_cCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_cCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerCloverTimeSlice/LQCD_cCmdSize);
    max_route(act, "toLmemMux_fromCPU", "toLmemMux");
    max_lmem_set_interrupt_on(act, "toLmem0");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel1");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void transfer_fix_clovers_to_dfe (bqcd_fix_clover *t_in, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_input(act, "data_in", t_in, VOLUME/2 * sizeof(bqcd_fix_clover));
	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerCloverTimeSlice/LQCD_cCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_cCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerCloverTimeSlice/LQCD_cCmdSize);
    max_route(act, "toLmemMux_fromCPU", "toLmemMux");
    max_lmem_set_interrupt_on(act, "toLmem0");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel1");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel1");
    max_ignore_route(act, "spfromLmem0Demux");

	max_run(engine, act);
	max_actions_free(act);
}

void apply_dirac (dirac_run_parameters *run_params) {

	max_actions_t *act = max_actions_init(maxfile, 0);

	max_set_uint64t(act, "diracKernel", "ieo", run_params->even_odd);
	max_set_uint64t(act, "diracKernel", "enableAxpyPost", run_params->apply_axpy_post);
	max_set_uint64t(act, "diracKernel", "dagger", run_params->dagger);
	max_set_uint64t(act, "diracKernel", "applyDslash", run_params->apply_dslash);
	max_set_uint64t(act, "diracKernel", "applyClover", run_params->apply_clover);
	max_set_uint64t(act, "diracKernel", "miscMode", run_params->misc_mode);
	max_set_uint64t(act, "diracKernel", "enableAxpyPre", run_params->apply_axpy_pre);
	max_set_uint64t(act, "diracKernel", "enableOutput2", run_params->enable_second_output);
	max_set_uint64t(act, "diracKernel", "cloverDiagShift", run_params->clover_diag_shift);
	max_set_uint64t(act, "diracKernel", "cloverOffdiagShift", run_params->clover_offdiag_shift);
	max_set_uint64t(act, "diracKernel", "gaugeShift", run_params->gauge_shift);
	//max_set_uint64t(act, "diracKernel", "s0ExpDiff", run_params->in0_exp_difference);
	if (run_params->apply_axpy_post)
	{
		max_set_double(act, "diracKernel", "axpyPostMultiplier", run_params->axpy_post_mult);
		/*max_set_double(act, "diracKernel", "alpha", 4+mass);*/
		/*max_set_double(act, "diracKernel", "beta_s",   beta_s / (16+4*mass));
		max_set_double(act, "diracKernel", "beta_t_b", beta_t_b / (16+4*mass));
		max_set_double(act, "diracKernel", "beta_t_f", beta_t_f / (16+4*mass));*/
	}
	else
	{
		max_ignore_scalar(act, "diracKernel", "axpyPostMultiplier");
		/*max_ignore_scalar(act, "diracKernel", "alpha");*/
		/*max_set_double(act, "diracKernel", "beta_s",   beta_s);
		max_set_double(act, "diracKernel", "beta_t_b", beta_t_b);
		max_set_double(act, "diracKernel", "beta_t_f", beta_t_f);*/
	}
	if (run_params->misc_mode)
	{
	//	max_set_uint64t(act, "diracKernel", "s2ExpDiff", run_params->in2_exp_difference);
	//	max_set_uint64t(act, "diracKernel", "s3ExpDiff", run_params->in3_exp_difference);
	}
	else
	{
	//	max_ignore_scalar(act, "diracKernel", "s2ExpDiff");
	//	max_ignore_scalar(act, "diracKernel", "s3ExpDiff");
	}
	if(run_params->apply_axpy_pre)
	{
		max_set_double(act, "diracKernel", "axpyPreMultiplier", run_params->axpy_pre_mult);
	}
	else
	{
		max_ignore_scalar(act, "diracKernel", "axpyPreMultiplier");
	}


	max_set_ticks(act, "diracKernel", 16/NUM_PIPES * (LT+2)*LX*LY*LZ/2 );

	max_set_ticks(act, "readCmdKernel0", (LT)*burstsPerGaugeTimeSlice/LQCD_gCmdSize);
	max_set_uint64t(act, "readCmdKernel0", "startAddress", run_params->gauge_address);
	max_set_uint64t(act, "readCmdKernel0", "halos", /*1*/0);
	max_set_uint64t(act, "readCmdKernel0", "cmdSize", LQCD_gCmdSize);
	max_set_uint64t(act, "readCmdKernel0", "burstsPerTSlice", burstsPerGaugeTimeSlice/LQCD_gCmdSize);

	max_set_ticks(act, "readCmdKernel1", (LT+2)*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "startAddress", run_params->in_spinor0_address);
	max_set_uint64t(act, "readCmdKernel1", "halos", 1);
	max_set_uint64t(act, "readCmdKernel1", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	if (run_params->enable_second_input == 0) {
		//max_ignore_scalar(act, "diracKernel", "s1ExpDiff");
		max_ignore_kernel(act, "readCmdKernel2");
	} else  {
	//	max_set_uint64t(act, "diracKernel", "s1ExpDiff", -run_params->in1_exp_difference);
		max_set_ticks(act, "readCmdKernel2", (LT + (run_params->apply_axpy_pre ? 2 : 0))
														*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
		max_set_uint64t(act, "readCmdKernel2", "startAddress", run_params->in_spinor1_address);
		max_set_uint64t(act, "readCmdKernel2", "halos", run_params->apply_axpy_pre ? 1 : 0);
		max_set_uint64t(act, "readCmdKernel2", "cmdSize", LQCD_spCmdSize);
		max_set_uint64t(act, "readCmdKernel2", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	}
	if (run_params->enable_clover_input == 0)
	{
		max_ignore_kernel(act, "readCmdKernel3");
	}
	else
	{
		max_set_ticks(act, "readCmdKernel3", (LT+2)*burstsPerCloverTimeSlice/LQCD_cCmdSize);
		max_set_uint64t(act, "readCmdKernel3", "startAddress", run_params->clover_address);
		max_set_uint64t(act, "readCmdKernel3", "halos", run_params->dagger ? 1 : 0);
		max_set_uint64t(act, "readCmdKernel3", "cmdSize", LQCD_cCmdSize);
		max_set_uint64t(act, "readCmdKernel3", "burstsPerTSlice", burstsPerCloverTimeSlice/LQCD_cCmdSize);
	}

	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", run_params->out_spinor0_address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	if(run_params->enable_second_output)
	{
		max_set_ticks(act, "writeCmdKernel1", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
		max_set_uint64t(act, "writeCmdKernel1", "startAddress", run_params->out_spinor1_address);
		max_set_uint64t(act, "writeCmdKernel1", "halos", 0);
		max_set_uint64t(act, "writeCmdKernel1", "cmdSize", LQCD_spCmdSize);
		max_set_uint64t(act, "writeCmdKernel1", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	}
	else
	{
		max_ignore_kernel(act, "writeCmdKernel1");
	}

	max_route(act, "toLmemMux_fromKernel", "toLmemMux");
	max_route(act, "spfromLmem0Demux", "spfromLmem0Demux_toKernel");
	max_lmem_set_interrupt_on(act, "toLmem0");

	struct timeval start_time, end_time;
	uint64_t t_out_sq;
	uint64_t t_out_sq_shift;
	//uint64_t t_out_exp0;
	//uint64_t t_out_exp1;
	max_get_uint64t(act, "diracKernel", "outSq", &t_out_sq);
	max_get_uint64t(act, "diracKernel", "outSqShift", &t_out_sq_shift);
	//max_get_uint64t(act, "diracKernel", "expShift0", &t_out_exp0);
	//max_get_uint64t(act, "diracKernel", "expShift1", &t_out_exp1);

	gettimeofday(&start_time, NULL);
	max_run(engine, act);
	gettimeofday(&end_time, NULL);

	//(run_params->out0_exp)=(int)t_out_exp0;
	//(run_params->out1_exp)=(int)t_out_exp1;
	(run_params->out_sq)=t_out_sq;
	(run_params->out_sq_shift)=(int)t_out_sq_shift;

	elapsedTime += (end_time.tv_sec - start_time.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (end_time.tv_usec - start_time.tv_usec) / 1000.0;   // us to ms
	max_actions_free(act);

}

void apply_misc_ops (dirac_run_parameters *run_params)
{
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_set_uint64t(act, "diracKernel", "applyDslash", 0);
	max_set_uint64t(act, "diracKernel", "applyClover", 0);

	max_set_uint64t(act, "diracKernel", "enableAxpyPost", 1);
	max_set_uint64t(act, "diracKernel", "enableAxpyPre", 1);
	max_set_uint64t(act, "diracKernel", "miscMode", 1);
	max_set_uint64t(act, "diracKernel", "enableOutput2", 1);

	max_set_double(act, "diracKernel", "axpyPostMultiplier", run_params->axpy_post_mult);
	max_set_double(act, "diracKernel", "axpyPreMultiplier", run_params->axpy_pre_mult);

	//max_set_uint64t(act, "diracKernel", "s0ExpDiff", -run_params->in0_exp_difference);
	//max_set_uint64t(act, "diracKernel", "s1ExpDiff", -run_params->in1_exp_difference);
	//max_set_uint64t(act, "diracKernel", "s2ExpDiff", -run_params->in2_exp_difference);
	//max_set_uint64t(act, "diracKernel", "s3ExpDiff", -run_params->in3_exp_difference);

	max_ignore_scalar(act, "diracKernel", "ieo");
	max_ignore_scalar(act, "diracKernel", "dagger");
	max_ignore_scalar(act, "diracKernel", "cloverDiagShift");
	max_ignore_scalar(act, "diracKernel", "cloverOffdiagShift");
	max_ignore_scalar(act, "diracKernel", "gaugeShift");

	max_set_ticks(act, "diracKernel", LT*LX*LY*LZ/2 );

	max_set_ticks(act, "readCmdKernel0", LT*burstsPerGaugeTimeSlice/LQCD_gCmdSize);
	max_set_uint64t(act, "readCmdKernel0", "startAddress", run_params->gauge_address);
	max_set_uint64t(act, "readCmdKernel0", "halos", 0);
	max_set_uint64t(act, "readCmdKernel0", "cmdSize", LQCD_gCmdSize);
	max_set_uint64t(act, "readCmdKernel0", "burstsPerTSlice", burstsPerGaugeTimeSlice/LQCD_gCmdSize);

	max_set_ticks(act, "readCmdKernel1", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "startAddress", run_params->in_spinor0_address);
	max_set_uint64t(act, "readCmdKernel1", "halos", 0);
	max_set_uint64t(act, "readCmdKernel1", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	max_set_ticks(act, "readCmdKernel2", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel2", "startAddress", run_params->in_spinor1_address);
	max_set_uint64t(act, "readCmdKernel2", "halos", 0);
	max_set_uint64t(act, "readCmdKernel2", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel2", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	max_set_ticks(act, "readCmdKernel3", LT*burstsPerCloverTimeSlice/LQCD_cCmdSize);
	max_set_uint64t(act, "readCmdKernel3", "startAddress", run_params->clover_address);
	max_set_uint64t(act, "readCmdKernel3", "halos", 0);
	max_set_uint64t(act, "readCmdKernel3", "cmdSize", LQCD_cCmdSize);
	max_set_uint64t(act, "readCmdKernel3", "burstsPerTSlice", burstsPerCloverTimeSlice/LQCD_cCmdSize);

	max_set_ticks(act, "writeCmdKernel0", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "startAddress", run_params->out_spinor0_address);
	max_set_uint64t(act, "writeCmdKernel0", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel0", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel0", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	max_set_ticks(act, "writeCmdKernel1", LT*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel1", "startAddress", run_params->out_spinor1_address);
	max_set_uint64t(act, "writeCmdKernel1", "halos", 0);
	max_set_uint64t(act, "writeCmdKernel1", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "writeCmdKernel1", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	max_route(act, "toLmemMux_fromKernel", "toLmemMux");
	max_route(act, "spfromLmem0Demux", "spfromLmem0Demux_toKernel");
	max_lmem_set_interrupt_on(act, "toLmem0");

	uint64_t t_out_sq;
	uint64_t t_out_sq_shift;
	//uint64_t t_out_exp0;
	//uint64_t t_out_exp1;
	max_get_uint64t(act, "diracKernel", "outSq", &t_out_sq);
	max_get_uint64t(act, "diracKernel", "outSqShift", &t_out_sq_shift);
	//max_get_uint64t(act, "diracKernel", "expShift0", &t_out_exp0);
	//max_get_uint64t(act, "diracKernel", "expShift1", &t_out_exp1);

	max_run(engine, act);

	//(run_params->out0_exp)=(int)t_out_exp0;
	//(run_params->out1_exp)=(int)t_out_exp1;
	(run_params->out_sq)=t_out_sq;
	(run_params->out_sq_shift)=(int)t_out_sq_shift;

	max_actions_free(act);
}

void transfer_spinors_to_cpu (bqcd_spinor *out, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_output(act, "data_out", out,  VOLUME/2 * sizeof(bqcd_spinor));

	max_set_ticks(act, "readCmdKernel1", (LT)*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "startAddress", address);
	max_set_uint64t(act, "readCmdKernel1", "halos", 0);
	max_set_uint64t(act, "readCmdKernel1", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	max_route(act, "spfromLmem0Demux", "spfromLmem0Demux_toCPU");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel0");
	max_ignore_kernel(act, "writeCmdKernel1");
	max_ignore_route(act, "toLmemMux");

	max_run(engine, act);
	max_actions_free(act);

}

void transfer_fix_spinors_to_cpu (bqcd_fix_spinor *out, int address) {
	max_actions_t *act = max_actions_init(maxfile, 0);
	max_queue_output(act, "data_out", out,  VOLUME/2 * sizeof(bqcd_fix_spinor));

	max_set_ticks(act, "readCmdKernel1", (LT)*burstsPerSpinorTimeSlice/LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "startAddress", address);
	max_set_uint64t(act, "readCmdKernel1", "halos", 0);
	max_set_uint64t(act, "readCmdKernel1", "cmdSize", LQCD_spCmdSize);
	max_set_uint64t(act, "readCmdKernel1", "burstsPerTSlice", burstsPerSpinorTimeSlice/LQCD_spCmdSize);

	max_route(act, "spfromLmem0Demux", "spfromLmem0Demux_toCPU");

	max_ignore_kernel(act, "diracKernel");
	max_ignore_kernel(act, "readCmdKernel0");
	max_ignore_kernel(act, "readCmdKernel2");
	max_ignore_kernel(act, "readCmdKernel3");
	max_ignore_kernel(act, "writeCmdKernel0");
	max_ignore_kernel(act, "writeCmdKernel1");
	max_ignore_route(act, "toLmemMux");

	max_run(engine, act);
	max_actions_free(act);

}


void unload_dfe(void) {
	max_unload(engine);
}

/************* Verifying Functions ***************/
float calc_normsq (bqcd_spinor* t_in, int V)
{
	FILE * outFile=fopen("calc_normsq_out.txt", "w");
	float sum=0.;
	for(int i=0; i<V; i++)
	{
		float to_sum = creal(t_in[i].s0c0)*creal(t_in[i].s0c0) + cimag(t_in[i].s0c0)*cimag(t_in[i].s0c0)
			+  creal(t_in[i].s0c1)*creal(t_in[i].s0c1) + cimag(t_in[i].s0c1)*cimag(t_in[i].s0c1)
			+  creal(t_in[i].s0c2)*creal(t_in[i].s0c2) + cimag(t_in[i].s0c2)*cimag(t_in[i].s0c2)
			+  creal(t_in[i].s1c0)*creal(t_in[i].s1c0) + cimag(t_in[i].s1c0)*cimag(t_in[i].s1c0)
			+  creal(t_in[i].s1c1)*creal(t_in[i].s1c1) + cimag(t_in[i].s1c1)*cimag(t_in[i].s1c1)
			+  creal(t_in[i].s1c2)*creal(t_in[i].s1c2) + cimag(t_in[i].s1c2)*cimag(t_in[i].s1c2)
			+  creal(t_in[i].s2c0)*creal(t_in[i].s2c0) + cimag(t_in[i].s2c0)*cimag(t_in[i].s2c0)
			+  creal(t_in[i].s2c1)*creal(t_in[i].s2c1) + cimag(t_in[i].s2c1)*cimag(t_in[i].s2c1)
			+  creal(t_in[i].s2c2)*creal(t_in[i].s2c2) + cimag(t_in[i].s2c2)*cimag(t_in[i].s2c2)
			+  creal(t_in[i].s3c0)*creal(t_in[i].s3c0) + cimag(t_in[i].s3c0)*cimag(t_in[i].s3c0)
			+  creal(t_in[i].s3c1)*creal(t_in[i].s3c1) + cimag(t_in[i].s3c1)*cimag(t_in[i].s3c1)
			+  creal(t_in[i].s3c2)*creal(t_in[i].s3c2) + cimag(t_in[i].s3c2)*cimag(t_in[i].s3c2);
		fprintf(outFile, "i=%d, to_sum=%f\n", i, to_sum);
		sum += to_sum;
	}
	return (float)sum;
}


int verify_results (bqcd_spinor* dfe_out, bqcd_spinor *expected_out, int V) {
	for (int i=0 ; i<V ; i++ ) {
		int error = compare_spinor(expected_out+i , dfe_out+i);
		if (error) {
			int idx=i;
			int t   = idx / ((LZ)*(LY)*(LX)/2);
			idx = idx % ((LZ)*(LY)*(LX)/2);
			int z   = idx / ((LY)*(LX)/2);
			idx = idx % ((LY)*(LX)/2);
			int y   = idx / ((LX)/2);
			int x   = idx % ((LX)/2);

			printf("Mismatch at spinor %d (x=%d, y=%d, z=%d, t=%d), value %d !\n",
					i, x, y, z, t, error);
			bqcd_spinor *a = expected_out + i;
			bqcd_spinor *b = dfe_out + i;
			printf("Expected: ( %E , %E ); Got: ( %E , %E )\n", creal(a->s0c0), cimag(a->s0c0),
					                  creal(b->s0c0), cimag(b->s0c0));
			return 1;
		}
	}
	printf("Good.\n");
	return 0;

}

float *calc_analysis_numbers(bqcd_spinor* dfe_out, bqcd_spinor *expected_out, int V)
{
	float *to_ret=malloc(8*sizeof(float));
	float max_abs_error=0;
	float max_rel_error=0;
	float avg_abs_error=0;
	float avg_rel_error=0;

	float max_abs_error_v1=0;
	float max_abs_error_v2=0;
	float max_rel_error_v1=0;
	float max_rel_error_v2=0;

	for (int i=0 ; i<V ; i++ ) {
		bqcd_spinor *a = expected_out + i;
		bqcd_spinor *b = dfe_out + i;
		for(int j=0; j<24; j++)
		{
			complex float ae;
			complex float be;
			switch(j/2)
			{
				case 0:  ae=a->s0c0; be=b->s0c0; break;
				case 1:  ae=a->s0c1; be=b->s0c1; break;
				case 2:  ae=a->s0c2; be=b->s0c2; break;
				case 3:  ae=a->s1c0; be=b->s1c0; break;
				case 4:  ae=a->s1c1; be=b->s1c1; break;
				case 5:  ae=a->s1c2; be=b->s1c2; break;
				case 6:  ae=a->s2c0; be=b->s2c0; break;
				case 7:  ae=a->s2c1; be=b->s2c1; break;
				case 8:  ae=a->s2c2; be=b->s2c2; break;
				case 9:  ae=a->s3c0; be=b->s3c0; break;
				case 10: ae=a->s3c1; be=b->s3c1; break;
				case 11: ae=a->s3c2; be=b->s3c2; break;
			}
			float av, bv;
			switch(j%2)
			{
				case 0:  av=creal(ae); bv=creal(be); break;
				case 1:  av=cimag(ae); bv=cimag(be); break;
			}

			float t_abs_error = fabs(av - bv);
			float t_rel_error = fabs((av - bv)/av);
/*			printf("av=%f; bv=%f\n", av, bv);
			printf("t_abs_err=%f; t_rel_err=%f\n", t_abs_error, t_rel_error);
*/
			if( (t_abs_error>max_abs_error) | ((i==0) & (j==0)) )
			{
				max_abs_error=t_abs_error;
				max_abs_error_v1 = av;
				max_abs_error_v2 = bv;
			}
			if( (t_rel_error>max_rel_error) | ((i==0) & (j==0)) )
			{
				max_rel_error=t_rel_error;
				max_rel_error_v1 = av;
				max_rel_error_v2 = bv;
			}

			avg_abs_error+=t_abs_error/(24*V);
			avg_rel_error+=t_rel_error/(24*V);
		}
	}

	to_ret[0]=max_abs_error;
	to_ret[1]=max_abs_error_v1;
	to_ret[2]=max_abs_error_v2;
	to_ret[3]=max_rel_error;
	to_ret[4]=max_rel_error_v1;
	to_ret[5]=max_rel_error_v2;
	to_ret[6]=avg_abs_error;
	to_ret[7]=avg_rel_error;

	return to_ret;
}

int verify_results_extended (bqcd_spinor* dfe_out, bqcd_spinor *expected_out, int V) {

	float *numbers=calc_analysis_numbers(dfe_out, expected_out, V);
	float max_abs_error=numbers[0];
	float max_abs_error_v1=numbers[1];
	float max_abs_error_v2=numbers[2];
	float max_rel_error=numbers[3];
	float max_rel_error_v1=numbers[4];
	float max_rel_error_v2=numbers[5];
	float avg_abs_error=numbers[6];
	float avg_rel_error=numbers[7];
	free(numbers);


/*	for (int i=0 ; i<V ; i++ ) {
		bqcd_spinor *a = expected_out + i;
		bqcd_spinor *b = dfe_out + i;
		for(int j=0; j<24; j++)
		{
			complex float ae;
			complex float be;
			switch(j/2)
			{
				case 0:  ae=a->s0c0; be=b->s0c0; break;
				case 1:  ae=a->s0c1; be=b->s0c1; break;
				case 2:  ae=a->s0c2; be=b->s0c2; break;
				case 3:  ae=a->s1c0; be=b->s1c0; break;
				case 4:  ae=a->s1c1; be=b->s1c1; break;
				case 5:  ae=a->s1c2; be=b->s1c2; break;
				case 6:  ae=a->s2c0; be=b->s2c0; break;
				case 7:  ae=a->s2c1; be=b->s2c1; break;
				case 8:  ae=a->s2c2; be=b->s2c2; break;
				case 9:  ae=a->s3c0; be=b->s3c0; break;
				case 10: ae=a->s3c1; be=b->s3c1; break;
				case 11: ae=a->s3c2; be=b->s3c2; break;
			}
			float av, bv;
			switch(j%2)
			{
				case 0:  av=creal(ae); bv=creal(be); break;
				case 1:  av=cimag(ae); bv=cimag(be); break;
			}

			float t_abs_error = fabs(av - bv);
			float t_rel_error = fabs((av - bv)/av);

			if( (t_abs_error>max_abs_error) | ((i==0) & (j==0)) )
			{
				max_abs_error=t_abs_error;
				max_abs_error_v1 = av;
				max_abs_error_v2 = bv;
			}
			if( (t_rel_error>max_rel_error) | ((i==0) & (j==0)) )
			{
				max_rel_error=t_rel_error;
				max_rel_error_v1 = av;
				max_rel_error_v2 = bv;
			}

			avg_abs_error+=t_abs_error/(24*V);
			avg_rel_error+=t_rel_error/(24*V);
		}
	}
*/
	printf ("max abs error = %E; occurred for av=%E and bv=%E.\n",
			max_abs_error, max_abs_error_v1, max_abs_error_v2);
	printf ("max rel error = %E; occurred for av=%E and bv=%E.\n",
			max_rel_error, max_rel_error_v1, max_rel_error_v2);
	printf ("avg abs error = %E; avg rel error=%E.\n",
			avg_abs_error, avg_rel_error);

	if(max_abs_error>.001)
	{
		printf("Spinor is bad!\n");
		verify_results(dfe_out, expected_out, V);
		return 1;
	}
	printf("Spinor is good.\n");
	return 0;
}

void output_spinor_analysis (char *filename, bqcd_spinor* dfe_out,
		bqcd_spinor *expected_out, int V)
{
	float *numbers=calc_analysis_numbers(dfe_out, expected_out, V);
	float max_abs_error=numbers[0];
	float max_abs_error_v1=numbers[1];
	float max_abs_error_v2=numbers[2];
	float max_rel_error=numbers[3];
	float max_rel_error_v1=numbers[4];
	float max_rel_error_v2=numbers[5];
	float avg_abs_error=numbers[6];
	float avg_rel_error=numbers[7];
	free(numbers);

	FILE *outfile=fopen(filename, "w");
	fprintf (outfile, "max abs error = %E; occurred for av=%E and bv=%E.\n",
			max_abs_error, max_abs_error_v1, max_abs_error_v2);
	fprintf (outfile, "max rel error = %E; occurred for av=%E and bv=%E.\n",
			max_rel_error, max_rel_error_v1, max_rel_error_v2);
	fprintf (outfile, "avg abs error = %E; avg rel error=%E.\n",
			avg_abs_error, avg_rel_error);
	if(max_abs_error>.001)
	{
		fprintf(outfile, "Spinor is bad!\n");
		verify_results(dfe_out, expected_out, V);
	}
	else
	{
		fprintf(outfile, "Spinor is good.\n");
	}

	fclose(outfile);
}


int AreNotSameComplex(complex float a, float complex b)
{
    return ( fabs(creal(a) - creal(b)) > 0.001 ) ||
    	   ( cimag(creal(a) - cimag(b)) > 0.001 );
}

int compare_spinor (bqcd_spinor *a, bqcd_spinor *b) {
	if (AreNotSameComplex(a->s0c0 ,b->s0c0)) return 1;
	if (AreNotSameComplex(a->s0c1 ,b->s0c1)) return 2;
	if (AreNotSameComplex(a->s0c2 ,b->s0c2)) return 3;
	if (AreNotSameComplex(a->s1c0 ,b->s1c0)) return 4;
	if (AreNotSameComplex(a->s1c1 ,b->s1c1)) return 5;
	if (AreNotSameComplex(a->s1c2 ,b->s1c2)) return 6;
	if (AreNotSameComplex(a->s2c0 ,b->s2c0)) return 7;
	if (AreNotSameComplex(a->s2c1 ,b->s2c1)) return 8;
	if (AreNotSameComplex(a->s2c2 ,b->s2c2)) return 9;
	if (AreNotSameComplex(a->s3c0 ,b->s3c0)) return 10;
	if (AreNotSameComplex(a->s3c1 ,b->s3c1)) return 11;
	if (AreNotSameComplex(a->s3c2 ,b->s3c2)) return 12;
	return 0;
}

/************ Utility Functions *********************/

void create_random_input(spinor* s, su3* u) {
	for (int i = 0; i < VOLUME / 2; i++) {
		create_random_spinor(s + i);
	}
	for (int i = 0; i < VOLUME * 4; i++) {
		create_random_su3(u + i);
	}
}

void create_random_spinor(spinor * s) {
	create_random_su3vector(&s->s0);
	create_random_su3vector(&s->s1);
	create_random_su3vector(&s->s2);
	create_random_su3vector(&s->s3);
}

void create_random_su3vector(su3_vector *v) {
	create_random_complex(&v->c0);
	create_random_complex(&v->c1);
	create_random_complex(&v->c2);
}

void create_random_su3(su3 *s) {
	create_random_complex(&s->c00);
	create_random_complex(&s->c01);
	create_random_complex(&s->c02);
	create_random_complex(&s->c10);
	create_random_complex(&s->c11);
	create_random_complex(&s->c12);
	create_random_complex(&s->c20);
	create_random_complex(&s->c21);
	create_random_complex(&s->c22);
}

void create_random_complex(float complex *a) {
	float r[2];
	for (int i = 0; i < 2; i++) {
		r[i] = (float) (rand()) / RAND_MAX * 2;
	}
	*a = r[0] + I * r[1];
}

void print_spinors (spinor* s) {
	for (int i=0 ; i < VOLUME/2 ; i++ ) {
		printf("%f %f\n", creal(s->s0.c0), cimag(s->s0.c0));
		s++;
	}
}

void print_fix_spinors (bqcd_fix_spinor* s)
{
	for (int i=0 ; i < VOLUME/2 ; i++ ) {
		printf("(%02x%02x%02x, %02x%02x%02x)\n", s[i].s0c0.real.bytes[0], s[i].s0c0.real.bytes[1],
		  s[i].s0c0.real.bytes[2], s[i].s0c0.imag.bytes[0], s[i].s0c0.imag.bytes[1],
		  s[i].s0c0.imag.bytes[2]);
		s++;
	}
}


void print_gauges (su3* s) {
	for (int i=0 ; i < VOLUME/2 * 4 ; i++ ) {
		printf("%f %f\n", creal(s->c00), cimag(s->c00));
		s++;
	}
}

/****************** Functions for reading data from files ***************/
void conv_float_to_fix24(fix24 *out, float t_in,  int exponent, int frac_bits)
{
	if(exponent<0)
	{
		//in/=(1<<abs((int)exponent));
		t_in/=pow(2.,abs(exponent));
	}
	else if (exponent>0)
	{
		//in*=(1<<(int)exponent);
		t_in*=pow(2.,exponent);
	}
	//int data = (int)round(in*((unsigned long long)1<<frac_bits));
	int data = (int)round(t_in*pow(2.,frac_bits));


//	printf("input= %f; data as int: %d\n", in, data);
	out->bytes[2]=(unsigned char)((data>>16) & 0xFF);
	out->bytes[1]=(unsigned char)((data>>8) & 0xFF);
	out->bytes[0]=(unsigned char)( data & 0xFF);
}

float conv_fix24_to_float(fix24 t_in, int exponent, int frac_bits)
{
	int in0=(int)t_in.bytes[0];
	int in1=(int)t_in.bytes[1];
	int in2=(int)t_in.bytes[2];
	unsigned int data = (unsigned int)in0 + ( (unsigned int)in1 << 8 )
								+ ( (unsigned int)in2 << 16 );

	int in32=(int)data;
	if(in32>8388607)
	{
		in32=in32-16777216;
	}

	float toret=(float)in32/pow(2.,exponent + frac_bits);

	return toret;
}

void conv_bqcd_spinor_to_bqcd_fix_spinor(bqcd_fix_spinor *out, bqcd_spinor *t_in, int exponent, int frac_bits)
{
	for(int i=0; i<VOLUME/2; ++i)
	{
		conv_float_to_fix24(&out[i].s0c0.real, creal(t_in[i].s0c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s0c0.imag, cimag(t_in[i].s0c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s0c1.real, creal(t_in[i].s0c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s0c1.imag, cimag(t_in[i].s0c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s0c2.real, creal(t_in[i].s0c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s0c2.imag, cimag(t_in[i].s0c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s1c0.real, creal(t_in[i].s1c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s1c0.imag, cimag(t_in[i].s1c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s1c1.real, creal(t_in[i].s1c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s1c1.imag, cimag(t_in[i].s1c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s1c2.real, creal(t_in[i].s1c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s1c2.imag, cimag(t_in[i].s1c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s2c0.real, creal(t_in[i].s2c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s2c0.imag, cimag(t_in[i].s2c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s2c1.real, creal(t_in[i].s2c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s2c1.imag, cimag(t_in[i].s2c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s2c2.real, creal(t_in[i].s2c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s2c2.imag, cimag(t_in[i].s2c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s3c0.real, creal(t_in[i].s3c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s3c0.imag, cimag(t_in[i].s3c0), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s3c1.real, creal(t_in[i].s3c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s3c1.imag, cimag(t_in[i].s3c1), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s3c2.real, creal(t_in[i].s3c2), exponent, frac_bits);
		conv_float_to_fix24(&out[i].s3c2.imag, cimag(t_in[i].s3c2), exponent, frac_bits);
	}
}

void conv_bqcd_fix_spinor_to_bqcd_spinor(bqcd_spinor *out, bqcd_fix_spinor *t_in, int exponent, int frac_bits)
{
	for(int i=0; i<VOLUME/2; ++i)
	{
		out[i].s0c0=	conv_fix24_to_float(t_in[i].s0c0.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s0c0.imag, exponent, frac_bits);
		out[i].s0c1=	conv_fix24_to_float(t_in[i].s0c1.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s0c1.imag, exponent, frac_bits);
		out[i].s0c2=	conv_fix24_to_float(t_in[i].s0c2.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s0c2.imag, exponent, frac_bits);
		out[i].s1c0=	conv_fix24_to_float(t_in[i].s1c0.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s1c0.imag, exponent, frac_bits);
		out[i].s1c1=	conv_fix24_to_float(t_in[i].s1c1.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s1c1.imag, exponent, frac_bits);
		out[i].s1c2=	conv_fix24_to_float(t_in[i].s1c2.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s1c2.imag, exponent, frac_bits);
		out[i].s2c0=	conv_fix24_to_float(t_in[i].s2c0.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s2c0.imag, exponent, frac_bits);
		out[i].s2c1=	conv_fix24_to_float(t_in[i].s2c1.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s2c1.imag, exponent, frac_bits);
		out[i].s2c2=	conv_fix24_to_float(t_in[i].s2c2.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s2c2.imag, exponent, frac_bits);
		out[i].s3c0=	conv_fix24_to_float(t_in[i].s3c0.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s3c0.imag, exponent, frac_bits);
		out[i].s3c1=	conv_fix24_to_float(t_in[i].s3c1.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s3c1.imag, exponent, frac_bits);
		out[i].s3c2=	conv_fix24_to_float(t_in[i].s3c2.real, exponent, frac_bits)
					 +I*conv_fix24_to_float(t_in[i].s3c2.imag, exponent, frac_bits);
	}
}

void shave_fix24_to_fix23(fix24 *out, fix24 *in)
{
	out->bytes[0]=in->bytes[0] & 0xFE;
	out->bytes[1]=in->bytes[1];
    out->bytes[2]=in->bytes[2];
}

void conv_fix24_to_fix23(fix24 *out, fix24 *in)
{
	unsigned char lsb=in->bytes[0] & 1;
	out->bytes[0]=in->bytes[0] & 0xFE;
	out->bytes[1]=in->bytes[1];
    out->bytes[2]=in->bytes[2];

    if(lsb)
    {
		if(out->bytes[0] < 254)
		{
			out->bytes[0]+=2;
		}
		else
		{
			if(out->bytes[1] < 255)
			{
				out->bytes[0]+=2;
				out->bytes[1]+=1;
			}
			else
			{
				out->bytes[0]+=2;
				out->bytes[1]+=1;
				out->bytes[2]+=1;
			}
		}
    }
}

unsigned char comp_fix24(fix24 a, fix24 b)
{
	int a_in0=(int)a.bytes[0];
	int a_in1=(int)a.bytes[1];
	int a_in2=(int)a.bytes[2];
	unsigned int a_data = (unsigned int)a_in0 + ( (unsigned int)a_in1 << 8 )
								+ ( (unsigned int)a_in2 << 16 );
	int a_in32=(int)a_data;
	if(a_in32>8388607)
	{
		a_in32=a_in32-16777216;
	}
	int b_in0=(int)b.bytes[0];
	int b_in1=(int)b.bytes[1];
	int b_in2=(int)b.bytes[2];
	unsigned int b_data = (unsigned int)b_in0 + ( (unsigned int)b_in1 << 8 )
								+ ( (unsigned int)b_in2 << 16 );
	int b_in32=(int)b_data;
	if(b_in32>8388607)
	{
		b_in32=b_in32-16777216;
	}

	return (abs(a_in32) > abs(b_in32));
}

void conv_enc_fix_spinor_to_bqcd_spinor(bqcd_spinor *out, bqcd_fix_spinor *in, int frac_bits)
{
	for(int k=0; k<VOLUME/2; ++k)
	{
		bqcd_fix_spinor *t_in = &in[k];

		fix24 data[24] = { 	t_in->s0c0.real, t_in->s0c0.imag, t_in->s1c0.real, t_in->s1c0.imag,
							t_in->s2c0.real, t_in->s2c0.imag, t_in->s3c0.real, t_in->s3c0.imag,
							t_in->s0c1.real, t_in->s0c1.imag, t_in->s1c1.real, t_in->s1c1.imag,
							t_in->s2c1.real, t_in->s2c1.imag, t_in->s3c1.real, t_in->s3c1.imag,
							t_in->s0c2.real, t_in->s0c2.imag, t_in->s1c2.real, t_in->s1c2.imag,
							t_in->s2c2.real, t_in->s2c2.imag, t_in->s3c2.real, t_in->s3c2.imag };
		float to_comp_float[24];

		for(int i=0; i<24; ++i)
		{
			fix24 temp;
			shave_fix24_to_fix23(&temp, &data[i]);
			to_comp_float[i]=conv_fix24_to_float(temp, 0, frac_bits);
			to_comp_float[i]=fabs(to_comp_float[i]);
		}

		int max_ids_step0 [12];
		float max_vals_step0 [12];
		for(int i=0; i<12; ++i)
		{
			if(to_comp_float[2*i]>to_comp_float[2*i+1])
			{
				max_ids_step0[i]=2*i;
				max_vals_step0[i]=to_comp_float[2*i];
			}
			else
			{
				max_ids_step0[i]=2*i+1;
				max_vals_step0[i]=to_comp_float[2*i+1];
			}
		}

		unsigned char bits[6];
		for(int i=0; i<6; ++i)
		{
			int tribute_id;

			if(max_vals_step0[2*i] > max_vals_step0[2*i+1])
			{
				tribute_id = max_ids_step0[2*i];
			}
			else
			{
				tribute_id = max_ids_step0[2*i+1];
			}
			bits[i]=data[tribute_id].bytes[0]&1;
		}

		int biased_exp_shift=(bits[5]<<5)+(bits[4]<<4)+(bits[3]<<3)+(bits[2]<<2)
				+(bits[1]<<1)+(bits[0]<<0);
		if(biased_exp_shift>31)
		{
			biased_exp_shift = biased_exp_shift - 64;
		}
		int exp_shift=biased_exp_shift-16;

		//printf("exp shift: %d\n", exp_shift);

		float fdata[24];
		for(int i=0; i<24; ++i)
		{
			fdata[i]=conv_fix24_to_float(data[i], -exp_shift, frac_bits);
		}

		bqcd_spinor *t_out = &out[k];
		t_out->s0c0 = fdata[0]	+ I*fdata[1];
		t_out->s1c0 = fdata[2]	+ I*fdata[3];
		t_out->s2c0 = fdata[4]	+ I*fdata[5];
		t_out->s3c0 = fdata[6]	+ I*fdata[7];
		t_out->s0c1 = fdata[8]	+ I*fdata[9];
		t_out->s1c1 = fdata[10]	+ I*fdata[11];
		t_out->s2c1 = fdata[12]	+ I*fdata[13];
		t_out->s3c1 = fdata[14]	+ I*fdata[15];
		t_out->s0c2 = fdata[16]	+ I*fdata[17];
		t_out->s1c2 = fdata[18]	+ I*fdata[19];
		t_out->s2c2 = fdata[20]	+ I*fdata[21];
		t_out->s3c2 = fdata[22]	+ I*fdata[23];
	}
}

void conv_bqcd_spinor_to_enc_fix_spinor(bqcd_fix_spinor *out, bqcd_spinor *in, int frac_bits)
{
	for(int k=0; k<VOLUME/2; ++k)
	{
		bqcd_spinor *t_in = &in[k];

		float f_data[24] = {
				creal(t_in->s0c0), cimag(t_in->s0c0), creal(t_in->s1c0), cimag(t_in->s1c0),
				creal(t_in->s2c0), cimag(t_in->s2c0), creal(t_in->s3c0), cimag(t_in->s3c0),
				creal(t_in->s0c1), cimag(t_in->s0c1), creal(t_in->s1c1), cimag(t_in->s1c1),
				creal(t_in->s2c1), cimag(t_in->s2c1), creal(t_in->s3c1), cimag(t_in->s3c1),
				creal(t_in->s0c2), cimag(t_in->s0c2), creal(t_in->s1c2), cimag(t_in->s1c2),
				creal(t_in->s2c2), cimag(t_in->s2c2), creal(t_in->s3c2), cimag(t_in->s3c2)
		};
		int exponent = -128;
		for(int i=0; i<24; ++i)
		{
			int tlog = floor(log2(fabs(f_data[i])));
			if(tlog>exponent) { exponent = tlog; }
		}
		exponent=-exponent;
		fix24 data[24];
		fix24 to_comp_data[24];
		for(int i=0; i<24; ++i)
		{
			conv_float_to_fix24(&data[i], f_data[i], exponent, frac_bits);
			shave_fix24_to_fix23(&to_comp_data[i], &data[i]);
		}

		int max_ids_step0 [12];
		fix24 max_vals_step0 [12];
		for(int i=0; i<12; ++i)
		{
			if(comp_fix24(to_comp_data[2*i], to_comp_data[2*i+1]))
			{
				max_ids_step0[i]=2*i;
				max_vals_step0[i]=to_comp_data[2*i];
			}
			else
			{
				max_ids_step0[i]=2*i+1;
				max_vals_step0[i]=to_comp_data[2*i+1];
			}
		}

		unsigned char bits[6];
		int biased_exponent = -exponent + 16;
		bits[0]=biased_exponent & 1;
		bits[1]=(biased_exponent & 2)>>1;
		bits[2]=(biased_exponent & 4)>>2;
		bits[3]=(biased_exponent & 8)>>3;
		bits[4]=(biased_exponent & 16)>>4;
		bits[5]=(biased_exponent & 32)>>5;

		for(int i=0; i<6; ++i)
		{
			int tribute_id;
			if(comp_fix24(max_vals_step0[2*i], max_vals_step0[2*i+1]))
			{
				tribute_id = max_ids_step0[2*i];
			}
			else
			{
				tribute_id = max_ids_step0[2*i+1];
			}
			fix24 temp;
			shave_fix24_to_fix23(&temp, &data[tribute_id]);
			temp.bytes[0]+=bits[i];
			data[tribute_id]=temp;
		}

		bqcd_fix_spinor *t_out = &out[k];
		t_out->s0c0.real = data[0];	t_out->s0c0.imag = data[1];
		t_out->s1c0.real = data[2];	t_out->s1c0.imag = data[3];
		t_out->s2c0.real = data[4];	t_out->s2c0.imag = data[5];
		t_out->s3c0.real = data[6];	t_out->s3c0.imag = data[7];
		t_out->s0c1.real = data[8];	t_out->s0c1.imag = data[9];
		t_out->s1c1.real = data[10];t_out->s1c1.imag = data[11];
		t_out->s2c1.real = data[12];t_out->s2c1.imag = data[13];
		t_out->s3c1.real = data[14];t_out->s3c1.imag = data[15];
		t_out->s0c2.real = data[16];t_out->s0c2.imag = data[17];
		t_out->s1c2.real = data[18];t_out->s1c2.imag = data[19];
		t_out->s2c2.real = data[20];t_out->s2c2.imag = data[21];
		t_out->s3c2.real = data[22];t_out->s3c2.imag = data[23];
	}

}

void conv_fix_su3_to_su3(su3 *out, fix_su3 *t_in, int exponent, int frac_bits)
{
	for(int i=0; i<8*VOLUME/2; ++i)
	{
		out[i].c00=	conv_fix24_to_float(t_in[i].c00.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c00.imag, exponent, frac_bits);
		out[i].c01=	conv_fix24_to_float(t_in[i].c01.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c01.imag, exponent, frac_bits);
		out[i].c02=	conv_fix24_to_float(t_in[i].c02.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c02.imag, exponent, frac_bits);
		out[i].c10=	conv_fix24_to_float(t_in[i].c10.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c10.imag, exponent, frac_bits);
		out[i].c11=	conv_fix24_to_float(t_in[i].c11.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c11.imag, exponent, frac_bits);
		out[i].c12=	conv_fix24_to_float(t_in[i].c12.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c12.imag, exponent, frac_bits);
		out[i].c20=	conv_fix24_to_float(t_in[i].c20.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c20.imag, exponent, frac_bits);
		out[i].c21=	conv_fix24_to_float(t_in[i].c21.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c21.imag, exponent, frac_bits);
		out[i].c22=	conv_fix24_to_float(t_in[i].c22.real, exponent, frac_bits)
							 +I*conv_fix24_to_float(t_in[i].c22.imag, exponent, frac_bits);
	}
}

void conv_su3_to_fix_su3(fix_su3 *out, su3 *t_in, int exponent, int frac_bits)
{
	for(int i=0; i<8*VOLUME/2; ++i)
	{
		conv_float_to_fix24(&out[i].c00.real, creal(t_in[i].c00), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c00.imag, cimag(t_in[i].c00), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c01.real, creal(t_in[i].c01), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c01.imag, cimag(t_in[i].c01), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c02.real, creal(t_in[i].c02), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c02.imag, cimag(t_in[i].c02), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c10.real, creal(t_in[i].c10), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c10.imag, cimag(t_in[i].c10), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c11.real, creal(t_in[i].c11), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c11.imag, cimag(t_in[i].c11), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c12.real, creal(t_in[i].c12), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c12.imag, cimag(t_in[i].c12), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c20.real, creal(t_in[i].c20), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c20.imag, cimag(t_in[i].c20), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c21.real, creal(t_in[i].c21), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c21.imag, cimag(t_in[i].c21), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c22.real, creal(t_in[i].c22), exponent, frac_bits);
		conv_float_to_fix24(&out[i].c22.imag, cimag(t_in[i].c22), exponent, frac_bits);
	}
}


void conv_bqcd_clover_to_bqcd_fix_clover(bqcd_fix_clover *out, bqcd_clover *t_in,
				int max_exp_diag, int max_exp_offdiag, int frac_bits)
{
	for(int i=0; i<VOLUME/2; ++i)
	{
		for(int k=0; k<18; ++k)
		{
			if(k==0 || k==10 || k==16)
			{
				conv_float_to_fix24(&out[i].c0[k].real, creal(t_in[i].c0[k]), max_exp_diag, frac_bits);
				conv_float_to_fix24(&out[i].c0[k].imag, cimag(t_in[i].c0[k]), max_exp_diag, frac_bits);
				conv_float_to_fix24(&out[i].c1[k].real, creal(t_in[i].c1[k]), max_exp_diag, frac_bits);
				conv_float_to_fix24(&out[i].c1[k].imag, cimag(t_in[i].c1[k]), max_exp_diag, frac_bits);
			}
			else
			{
				conv_float_to_fix24(&out[i].c0[k].real, creal(t_in[i].c0[k]), max_exp_offdiag, frac_bits);
				conv_float_to_fix24(&out[i].c0[k].imag, cimag(t_in[i].c0[k]), max_exp_offdiag, frac_bits);
				conv_float_to_fix24(&out[i].c1[k].real, creal(t_in[i].c1[k]), max_exp_offdiag, frac_bits);
				conv_float_to_fix24(&out[i].c1[k].imag, cimag(t_in[i].c1[k]), max_exp_offdiag, frac_bits);
			}
		}
	}
}

void read_spinor(char * filename, bqcd_spinor *out, int *maximum_exponent) {
	FILE * fptr = fopen(filename, "r" );
	char temp[64];
	*maximum_exponent=-128;

	for (int t=0 ; t<LT ; t++ ) {
		for (int z=0 ; z<LZ ; z++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int x=0 ; x<LX/2 ; x++ ) {

					bqcd_spinor *s = &out [ (((((t*LZ)+z)*LY)+y)*LX/2)+x ];

					fgets (temp, 64, fptr);

					float a[12], b[12];
					for(int i=0; i<12; i++)
					{
						fscanf(fptr, " ( %f , %f )\n", &a[i], &b[i]);
						int t_expa = floor(log2(a[i]));
						int t_expb = floor(log2(b[i]));
						if(t_expa>*maximum_exponent)
						{
							*maximum_exponent=t_expa;
						}
						if(t_expb>*maximum_exponent)
						{
							*maximum_exponent=t_expb;
						}
					}
					s->s0c0 = a[0] + I * b[0];
					s->s1c0 = a[1] + I * b[1];
					s->s2c0 = a[2] + I * b[2];
					s->s3c0 = a[3] + I * b[3];

					s->s0c1 = a[4] + I * b[4];
					s->s1c1 = a[5] + I * b[5];
					s->s2c1 = a[6] + I * b[6];
					s->s3c1 = a[7] + I * b[7];

					s->s0c2 = a[8] + I * b[8];
					s->s1c2 = a[9] + I * b[9];
					s->s2c2 = a[10] + I * b[10];
					s->s3c2 = a[11] + I * b[11];
				}
			}
		}
	}
	fclose(fptr);
}

void read_gauge(char * filename, su3 *s) {
	FILE * fptr = fopen(filename, "r" );
	char temp[64];
	for (int i = 0 ; i < VOLUME/2 * 8 ; i++ ) {
		fgets (temp, 64, fptr);

		for(int j=0; j<9; ++j)
		{
			float a, b;
			fscanf(fptr, "%f\n%f\n", &a, &b);
			switch(j)
			{
				case 0: s->c00 = a + I * b; break;
				case 1: s->c01 = a + I * b; break;
				case 2: s->c02 = a + I * b; break;
				case 3: s->c10 = a + I * b; break;
				case 4: s->c11 = a + I * b; break;
				case 5: s->c12 = a + I * b; break;
				case 6: s->c20 = a + I * b; break;
				case 7: s->c21 = a + I * b; break;
				case 8: s->c22 = a + I * b; break;
			}
		}
		s++;
	}
	fclose(fptr);
}

void read_clover(char * filename, bqcd_clover *result,
		int *max_exp_diagonal, int *max_exp_offdiagonal)
{
	FILE * fptr = fopen(filename, "r" );
	char temp[64];
	*max_exp_diagonal = -128;
	*max_exp_offdiagonal = -128;
	for(int i=0; i < VOLUME/2; i++)
	{
		for(int k=0; k<2; k++)
		{
			fgets (temp, 64, fptr);
			_Complex float *tC = (k==0 ? result[i].c0 : result[i].c1);
			for(int j=0; j<18; j++)
			{
				float a, b;
				fscanf(fptr, " ( %f , %f )\n", &a, &b);
				//printf("( %d , %d, %d )\n", i, k, j);
				tC[j]=a+I*b;
				int t_expa = floor(log2(fabs(a)));
				int t_expb = floor(log2(fabs(b)));
				if(j==0 || j==10 || j==16)
				{
					if(t_expa>*max_exp_diagonal) { *max_exp_diagonal=t_expa; }
					if(t_expb>*max_exp_diagonal) { *max_exp_diagonal=t_expb; }
					/*if(fabs(a)>max_diagonal) { max_diagonal = fabs(a); }
					if(fabs(b)>max_diagonal) { max_diagonal = fabs(b); }*/
				}
				else
				{
					if(t_expa>*max_exp_offdiagonal) { *max_exp_offdiagonal=t_expa; }
					if(t_expb>*max_exp_offdiagonal) { *max_exp_offdiagonal=t_expb; }
					/*if(fabs(a)>max_offdiagonal) { max_offdiagonal = fabs(a); }
					if(fabs(b)>max_offdiagonal) { max_offdiagonal = fabs(b); }*/
				}
			}
		}
	}
	fclose(fptr);
}

void write_spinor(char *filename, bqcd_spinor *to_write, char *format)
{
	FILE * fptr = fopen(filename, "w" );

	for (int t=0 ; t<LT ; t++ ) {
		for (int z=0 ; z<LZ ; z++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int x=0 ; x<LX/2 ; x++ ) {

					bqcd_spinor *s = &to_write [ (((((t*LZ)+z)*LY)+y)*LX/2)+x ];

					fprintf(fptr,"x=%d,y=%d,z=%d,t=%d\n", x, y, z, t);

					fprintf(fptr, format, creal(s->s0c0), cimag(s->s0c0));
					fprintf(fptr, format, creal(s->s1c0), cimag(s->s1c0));
					fprintf(fptr, format, creal(s->s2c0), cimag(s->s2c0));
					fprintf(fptr, format, creal(s->s3c0), cimag(s->s3c0));

					fprintf(fptr, format, creal(s->s0c1), cimag(s->s0c1));
					fprintf(fptr, format, creal(s->s1c1), cimag(s->s1c1));
					fprintf(fptr, format, creal(s->s2c1), cimag(s->s2c1));
					fprintf(fptr, format, creal(s->s3c1), cimag(s->s3c1));

					fprintf(fptr, format, creal(s->s0c2), cimag(s->s0c2));
					fprintf(fptr, format, creal(s->s1c2), cimag(s->s1c2));
					fprintf(fptr, format, creal(s->s2c2), cimag(s->s2c2));
					fprintf(fptr, format, creal(s->s3c2), cimag(s->s3c2));

				}
			}
		}
	}
	fclose(fptr);

}

void write_gauge(char *filename, su3 *s)
{
	FILE * fptr = fopen(filename, "w" );
	int i=0;

	for (int t=0 ; t<LT ; t++ )
	{
		for (int z=0 ; z<LZ ; z++ )
		{
			for (int y=0 ; y<LY ; y++ )
			{
				for (int x=0 ; x<LX/2 ; x++ )
				{
					for(int mu=0; mu<8; mu++)
					{
						su3 *ts=&s[i++];
						fprintf(fptr,"=== x=%d,y=%d,z=%d,t=%d, mu=%d ===\n", x, y, z, t, mu);
						fprintf(fptr,"   %10.7f+I*%10.7f    %10.7f+I*%10.7f    %10.7f+I*%10.7f\n",
								creal(ts->c00),cimag(ts->c00),creal(ts->c01),cimag(ts->c01),creal(ts->c02),cimag(ts->c02));
						fprintf(fptr,"(  %10.7f+I*%10.7f    %10.7f+I*%10.7f    %10.7f+I*%10.7f  )\n",
								creal(ts->c10),cimag(ts->c10),creal(ts->c11),cimag(ts->c11),creal(ts->c12),cimag(ts->c12));
						fprintf(fptr,"   %10.7f+I*%10.7f    %10.7f+I*%10.7f    %10.7f+I*%10.7f\n",
								creal(ts->c20),cimag(ts->c20),creal(ts->c21),cimag(ts->c21),creal(ts->c22),cimag(ts->c22));
					}
				}
			}
		}
	}
	fclose(fptr);
}

void write_clover(char * filename, bqcd_clover *to_write)
{
	FILE * fptr = fopen(filename, "w" );
	for(int i=0; i<VOLUME/2; i++)
	{
		for(int k=0; k<2; k++)
		{
			fprintf(fptr, "i=%d, k=%d\n", i, k);

			_Complex float *tC = k==0 ? to_write[i].c0 : to_write[i].c1;
			for(int j=0; j<18; j++)
			{
				fprintf(fptr, "( %10.7f , %10.7f )\n", creal(tC[j]), cimag(tC[j]));
			}
		}
	}
	fclose(fptr);
}

/****************** Functions for reordering data in memory ***************/
/*inline size_t index_from_coords(int x, int y, int z, int t)
{
	return ((t*LZ + z)*LY + y)*(LX/2) + x;

	// LZ*LY*(LX/2)*t + LY*(LX/2)*z + (LX/2)*y + x;
}*/


//void new_reorganize_gauge(su3 const * const in, su3 * const out, int ieo)
//{
//	for(int t=0; t<LT; t++)
//	{
//		for(int z=0; z<LZ; z++)
//		{
//			for(int y=0; y<LY; y++)
//			{
//				for(int x=0; x<LX/2; x++)
//				{
//					for(int mu=0; mu<4; mu++)
//					{
//						for(int dir=0; dir<2; dir++)
//						{
//							int offset = (dir==0) ? -1 : 1;
//							if(mu==0)
//							{
//								if(offset<0) { offset = 0; }
//								if(1 - ((t & 1) ^ (z & 1) ^ (y & 1) ^ ieo))
//								{
//									offset = 1-offset;
//								}
//							}
//							if(mu>1)
//							{
//								offset *= -1;
//							}
//
//							int nt = (mu == 3 ) ? (t+LT+offset)%LT	   : t;
//							int nz = (mu == 2 ) ? (z+LZ+offset)%LZ	   : z;
//							int ny = (mu == 1 ) ? (y+LY+offset)%LY	   : y;
//							int nx = (mu == 0 ) ? (x+LX+offset)%(LX/2) : x;
//
//							size_t out_index=8*index_from_coords(x,y,z,t) + 2*mu + dir;
//							size_t in_index=8*index_from_coords(nx,ny,nz,nt) + 2*mu + 1;
//
//
//							/*printf("(t=%d, z=%d, y=%d, x=%d, mu=%d) <- gauge(x=%d, y=%d, z=%d, t=%d, mu=%d)\n",
//								t, z, y, x, mu*2+dir,     nx, ny, nz, nt, nmu
//							);*/
//
//
//							out[out_index] = in[in_index];
//						}
//					}
//				}
//			}
//		}
//	}
//
//}

/*void new_reorganize_gauge(su3 const * const in, su3 * const out, int ieo)
{
	for(int t=0; t<LT; t++)
	{
		for(int z=0; z<LZ; z++)
		{
			for(int y=0; y<LY; y++)
			{
				for(int x=0; x<LX/2; x++)
				{
					for(int mu=0; mu<4; mu++)
					{
						for(int dir=0; dir<2; dir++)
						{
							int ndir= (mu == 1) ? dir : 1-dir;
							int offset = ndir == 0 ? -1 : 0;
							if(ieo==0) { ++offset; }

//							if(mu==0)
//							{
//								if(ieo==1)
//								{
//									if((t & 1) ^ (z & 1) ^ (y & 1) ^ ieo)
//									{
//										offset = ndir == 0 ? 1 : 0;
//									}
//									else
//									{
//										offset = 0;
//									}
//								}
//								else
//								{
//									if((t & 1) ^ (z & 1) ^ (y & 1) ^ ieo)
//									{
//										offset = 0;
//									}
//									else
//									{
//										offset = ndir == 0 ? 0 : 1;
//									}
//								}
//							}
							if(mu==0)
							{
								if(1 - ((t & 1) ^ (z & 1) ^ (y & 1)))
								{
									offset = ndir == 0 ? 1 : 0;
									if(ieo==0) { offset = 1-offset; }
								}
								else
								{
									offset = 0;
								}
							}
							int nt = (mu == 3 ) ? (t+LT+offset)%LT	   : t;
							int nz = (mu == 2 ) ? (z+LZ+offset)%LZ	   : z;
							int ny = (mu == 1 ) ? (y+LY+offset)%LY	   : y;
							int nx = (mu == 0 ) ? (x+LX+offset)%(LX/2) : x;

							size_t out_index=8*index_from_coords(x,y,z,t) + 2*mu + dir;
							size_t in_index=8*index_from_coords(nx,ny,nz,nt) + 2*mu + 1-ndir;

							if(ieo==0)
							{
								printf("(x=%d, y=%d, z=%d, t=%d, mu=%d) <- gauge(x=%d, y=%d, z=%d, t=%d, mu=%d)\n",
									x, y, z, t, 2*mu + dir,     nx, ny, nz, nt, 2*mu + 1-ndir
								);
							}


							out[out_index] = in[in_index];
						}
					}
					if(ieo==0)
					{
						printf("\n");
					}
				}
			}
		}
	}

}*/

//void reorganize_gauge(su3 const * const in, su3 * const out, int ieo) {
//	int i = 0;
//	for (int t=0 ; t<LT ; t++ ) {
//		for (int z=0 ; z<LZ ; z++ ) {
//			for (int y=0 ; y<LY ; y++ ) {
//				for (int x=0 ; x<LX/2 ; x++ ) {
//					for (int mu=0; mu<4 ; mu++ ) {
//						for (int f=-1; f<=1 ; f+= 2) {
//							su3 tmp = in[i];
//
//							int isOddRow = (t & 1) ^ (z & 1) ^ (y & 1) ^ ieo;
//
//							/*int mu_ = (mu+1)%4;
//							int t_ = t;               // converting from checkerboarded
//							int z_ = z/2;             // coordinates of qphix along x-axis
//							int y_ = y;               // to tmLQCD checkerboarding along
//							int x_ = (2*x)+isOddRow;  // along y-axis*/
//
//							if (f == 1) {
//
//								int xx = (mu==0)?( isOddRow ? x+1 :x ):x;
//								int yy = (mu==1)?y+1:y;
//								int zz = (mu==2)?z+1:z;
//								int tt = (mu==3)?t+1:t;
//
//								tt = (tt+LT) % LT;
//								zz = (zz+LZ) % LZ;
//								yy = (yy+LY) % LY;
//								xx = (xx+LX) % (LX/2);
//
//								out[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*8+mu*2+1 ] = tmp;
//
//							} else {
//
//								int xx = (mu==0)?( isOddRow ? x : x-1 ):x;
//								int yy = (mu==1)?y-1:y;
//								int zz = (mu==2)?z-1:z;
//								int tt = (mu==3)?t-1:t;
//
//								tt = (tt+LT) % LT;
//								zz = (zz+LZ) % LZ;
//								yy = (yy+LY) % LY;
//								xx = (xx+LX) % (LX/2);
//
//								out[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*8+mu*2+0 ] = tmp;
//
//							}
//
//							i++;
//						}
//					}
//				}
//			}
//		}
//	}
//}

void collate_gauges(su3* u0, su3* u1, su3 *ures, int ieo, int *maximum_exponent)
{
	*maximum_exponent=-128;
	for(size_t i=0; i<(size_t)VOLUME/2; i++)
	{
		for (size_t mu=0 ; mu<4 ; ++mu ) {
			for (size_t dir=0 ; dir<2 ; ++dir ) {
				if(dir==0)
				{
					ures[8*i+2*mu+dir] = u1[8*i+2*mu+ieo];
				}
				else
				{
					ures[8*i+2*mu+dir] = u0[8*i+2*mu+ieo];
				}
				for(int j=0; j<9; ++j)
				{
					complex float c;
					switch(j)
					{
						case 0: c=ures[8*i+2*mu+dir].c00; break;
						case 1: c=ures[8*i+2*mu+dir].c01; break;
						case 2: c=ures[8*i+2*mu+dir].c02; break;
						case 3: c=ures[8*i+2*mu+dir].c10; break;
						case 4: c=ures[8*i+2*mu+dir].c11; break;
						case 5: c=ures[8*i+2*mu+dir].c12; break;
						case 6: c=ures[8*i+2*mu+dir].c20; break;
						case 7: c=ures[8*i+2*mu+dir].c21; break;
						case 8: c=ures[8*i+2*mu+dir].c22; break;
					}
					int t_expa = floor(log2(fabs(creal(c))));
					int t_expb = floor(log2(fabs(cimag(c))));
					if(t_expa > *maximum_exponent)
					{
						*maximum_exponent = t_expa;
					}
					if(t_expb > *maximum_exponent)
					{
						*maximum_exponent = t_expb;
					}
				}

			}
		}
	}
}


/************* These functions are not needed anymore ************************/

void reorganize_ueven (su3 *out, su3 *in) {
	int i = 0;
	for (int t=0 ; t<LT ; t++ ) {
		for (int x=0 ; x<LX ; x++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int z=0 ; z<LZ/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1);
					for (int mu=0; mu<4 ; mu++ ) {
						int tt = (mu==0)?t+1:t;
						int xx = (mu==1)?x+1:x;
						int yy = (mu==2)?y+1:y;
						int zz = (mu==3)?( (isOddRow)?(z):z+1 ):z;

						tt = (tt+LT) % LT;
						xx = (xx+LX) % LX;
						yy = (yy+LY) % LY;
						zz = (zz+LZ) % (LZ/2);

						out[i] = in [ ((((((tt*LX)+xx)*LY)+yy)*LZ/2)+zz)*4+mu ];
						i++;
					}
				}
			}
		}
	}
}

void reorganize_back_ueven (su3 *out, su3 *in) {
	int i = 0;
	for (int t=0 ; t<LT ; t++ ) {
		for (int x=0 ; x<LX ; x++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int z=0 ; z<LZ/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1);
					for (int mu=0; mu<4 ; mu++ ) {
						int tt = (mu==0)?t+1:t;
						int xx = (mu==1)?x+1:x;
						int yy = (mu==2)?y+1:y;
						int zz = (mu==3)?( (isOddRow)?(z):z+1 ):z;

						tt = (tt+LT) % LT;
						xx = (xx+LX) % LX;
						yy = (yy+LY) % LY;
						zz = (zz+LZ) % (LZ/2);

						out[((((((tt*LX)+xx)*LY)+yy)*LZ/2)+zz)*4+mu ] = in [ i ];
						i++;
					}
				}
			}
		}
	}
}

/* Converting from qphix style gauge for the whole lattice, to even/odd separated
 * tmLQCD style gauge fields
 */
void devide_gauge_to_oddeven(su3 const * const in, su3 * const even, su3 * const odd, int ieo) {
	int i = 0;
	for (int t=0 ; t<LT ; t++ ) {
		for (int z=0 ; z<LZ ; z++ ) {
			for (int y=0 ; y<LY ; y++ ) {
				for (int x=0 ; x<LX/2 ; x++ ) {
					for (int mu=0; mu<4 ; mu++ ) {
						for (int f=-1; f<=1 ; f+= 2) {
							su3 tmp = in[i];

							int isOddRow = (t & 1) ^ (z & 1) ^ (y & 1) ^ ieo;

							/*int mu_ = (mu+1)%4;
							int t_ = t;               // converting from checkerboarded
							int z_ = z/2;             // coordinates of qphix along x-axis
							int y_ = y;               // to tmLQCD checkerboarding along
							int x_ = (2*x)+isOddRow;  // along y-axis*/

							if (f == 1) {

								int xx = (mu==0)?( isOddRow ? x+1 :x ):x;
								int yy = (mu==1)?y+1:y;
								int zz = (mu==2)?z+1:z;
								int tt = (mu==3)?t+1:t;

								tt = (tt+LT) % LT;
								zz = (zz+LZ) % LZ;
								yy = (yy+LY) % LY;
								xx = (xx+LX) % (LX/2);

								even[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*4+mu ] = tmp;

							} else {

								int xx = (mu==0)?( isOddRow ? x : x-1 ):x;
								int yy = (mu==1)?y-1:y;
								int zz = (mu==2)?z-1:z;
								int tt = (mu==3)?t-1:t;

								tt = (tt+LT) % LT;
								zz = (zz+LZ) % LZ;
								yy = (yy+LY) % LY;
								xx = (xx+LX) % (LX/2);

								odd[ ((((((tt*LZ)+zz)*LY)+yy)*LX/2)+xx)*4+mu ] = tmp;

							}

							i++;
						}
					}
				}
			}
		}
	}
}


void add_1d_halos_spinor(spinor* with_halos, spinor* orig, int halos) {

	for (int t = -halos ; t < LT+halos ; t++ ) {
		for (int z = 0 ; z < LZ ; z++ ) {
			for (int y = 0 ; y < LY ; y++ ) {
				for (int x = 0 ; x < LX/2 ; x++ ) {
					int tt = (t+LT)%LT;

					with_halos[ ((( (t+halos) * LZ + z) * LY + y) * LX/2 ) + x] =
							orig[ (((tt*LZ + z) * LY + y) * LX/2 ) + x];
				}
			}
		}
	}

}

void add_1d_halos_gauge(su3* with_halos, su3* orig, int halos) {

	for (int t = -halos ; t < LT+halos ; t++ ) {
		for (int z = 0 ; z < LZ ; z++ ) {
			for (int y = 0 ; y < LY ; y++ ) {
				for (int x = 0 ; x < LX/2 ; x++ ) {
					int tt = (t+LT)%LT;

					for (int i=0 ; i < 8 ; i++ ){
						with_halos[ (((( (t+halos) * LZ + z) * LY + y) * LX/2 ) + x) * 8 + i] =
								orig[ ((((tt*LZ + z) * LY + y) * LX/2 ) + x) * 8 + i];
					}
				}
			}
		}
	}

}


void add_4d_halos_spinor(spinor* with_halos, spinor* orig, int halos, int eo) {
	int lhx, lhy, lhz;
	lhx = LX+2*halos;
	lhy = LY+2*halos;
	lhz = LZ/2+halos;

	eo = eo ^ 1;

	for (int t = -halos ; t < LT+halos ; t++ ) {
		for (int x = -halos ; x < LX+halos ; x++ ) {
			for (int y = -halos ; y < LY+halos ; y++ ) {
				for (int z = -(halos/2) ; z < LZ/2+(halos+1)/2 ; z++ ) {
					int tt = (t+LT)%LT;
					int xx = (x+LX)%LX;
					int yy = (y+LY)%LY;
					int isOddRow = (tt & 1) ^ (xx & 1) ^ (yy & 1) ^ eo;
					int zz;
					if (halos%2 == 0) {
						zz = (z+LZ/2)%(LZ/2);
					} else {
						zz = (z+LZ*2-isOddRow*halos)%(LZ/2);
					}

					with_halos[ (((( (t+halos) * lhx +
					                 (x+halos)) * lhy +
					                 (y+halos)) * lhz ) +
					                 (z+halos/2))] =
							orig[ (((tt*LX + xx) * LY + yy) * LZ/2 ) + zz];
				}
			}
		}
	}
}

void add_4d_halos_gauge(su3* with_halos, su3* orig, int halos, int eo) {
	int lhx, lhy, lhz;
	lhx = LX+2*halos;
	lhy = LY+2*halos;
	lhz = LZ/2+halos;

	eo = eo ^ 1;

	for (int t = -halos ; t < LT+halos ; t++ ) {
		for (int x = -halos ; x < LX+halos ; x++ ) {
			for (int y = -halos ; y < LY+halos ; y++ ) {
				for (int z = -(halos/2) ; z < LZ/2+(halos+1)/2 ; z++ ) {
					int tt = (t+LT)%LT;
					int xx = (x+LX)%LX;
					int yy = (y+LY)%LY;
					int isOddRow = (tt & 1) ^ (xx & 1) ^ (yy & 1) ^ eo;
					int zz;
					if (halos%2 == 0) {
						zz = (z+LZ/2)%(LZ/2);
					} else {
						zz = (z+LZ*2-isOddRow*halos)%(LZ/2);
					}

					for (int i=0 ; i < 4 ; i++ ){
						with_halos[ ((((( (t+halos) * lhx +
				                          (x+halos)) * lhy +
				                          (y+halos)) * lhz ) +
				                          (z+halos/2)) * 4) + i] =
						      orig[ ((((tt*LX + xx) * LY + yy) * LZ/2 ) + zz) * 4 + i];
					}
				}
			}
		}
	}
}

