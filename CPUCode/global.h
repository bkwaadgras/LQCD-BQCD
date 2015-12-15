#ifndef GLOBAL_H
#define GLOBAL_H

#ifdef TEST_MAIN
#  define EXTERN
#else
#  define EXTERN extern
#endif

EXTERN int LT, VOLUME;
EXTERN int LX, LY, LZ;
EXTERN int NUM_PIPES, LOOP_OFFSET;


typedef struct {
	uint64_t in_spinor0_address;
	uint64_t in_spinor1_address;
	uint64_t out_spinor0_address;
	uint64_t out_spinor1_address;
	uint64_t gauge_address;
	uint64_t clover_address;

	unsigned char enable_second_input;
	unsigned char apply_axpy_pre;
	unsigned char apply_axpy_post;

	unsigned char apply_dslash;

	unsigned char enable_clover_input;
	unsigned char apply_clover;

	unsigned char enable_second_output;

	unsigned char even_odd;
	unsigned char dagger;
	unsigned char misc_mode;

	double axpy_pre_mult;
	double axpy_post_mult;

	uint64_t out_sq;
	uint64_t out_sq_shift;
	uint64_t clover_diag_shift;
	uint64_t clover_offdiag_shift;
	uint64_t gauge_shift;
} dirac_run_parameters;


#endif
