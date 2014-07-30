#ifndef LOW_COVERAGE_CONCATENATE_H
#define LOW_COVERAGE_CONCATENATE_H


#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "parser.h"
#include "util.h"
#include "low_coverage.h"



typedef struct Sorter
{
	unsigned int pos;
	unsigned int aln_start;
	unsigned int off_set;
} Sorter;



typedef struct Irreg_2DSorter
{
	unsigned int dim1;
	unsigned int max;
 	unsigned int step_size;
	unsigned int *dim2;
	struct Sorter **data;
} Irreg_2DSorter;



Irreg_2DSorter*
make_irregular_2D(unsigned int dim1);

void
free_Irreg2DSorter(Irreg_2DSorter *irr);



/**
 * \brief Function to produce a t_coffee alignment with concatenated low_coverage sequences.
 *
 * \param aln_f The alignment file
 * \param lib_f The library file
 * \param out_f The file to write the final alignment into
 * \param do_missing_pairs Should missing pairs be filled with lastz
 * \param ref_tree_f The reference tree
 */
int
concatenate_lib(char *aln_f, char *lib_f, char *out_f, int do_missing_pairs, char *ref_tree_f, int n_cores);





#endif //LOW_COVERAGE_CONCATENATE_H