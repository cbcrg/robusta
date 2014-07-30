#ifndef LOW_COVERAGE_H
#define LOW_COVERAGE_H


#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "parser.h"
#include "util.h"
#include "pw_aligners.h"
#include "low_coverage_concatenate.h"

/*! \file low_coverage.h
 * \brief Functions needed to produce low coverage alignments.
 *
 */



/**
 * \brief Produces a tree for a low coverage Alinment
 *
 * To prevent T-Coffee from aligning the short low coverage sequence pieces first it is neccessary to give a tree which alignes them at the end. The tree produced in this function is completeley unbalenced and just adds one sequence after another.
 * \param tree_f The file in which the tree should be stored.
 * \param num_seqs The number of sequences. It is supposed that the first sequences are from the reference alignment.
 */
void produce_tree(char *tree_f, unsigned int num_seqs);


/**
 * \brief Computes the library for low coverage genomes.
 *
 *\param aln_file The alignment file in XMFA format. The first block has to be the reference alignment. The rest pairwise alignments.
 \param lib_file The file into which the library will be written.
 \param result_seqs A pointer to which all Sequences found are written (contains only Meta information, not the sequences themselves).
 \param do_missing_blasts Starts a lastz call for each pair of low coverage genomes and reference sequence which is not found in the pair alignments.
 \return The number of sequences found.
 */
int produce_low_cow_library(char *aln_f, char *lib_f, Sequence ***result_seqs, int do_missing_blasts, char *ref_tree_f);


/**
 * \brief Computes a library for low coverage alignments.
 *
 * It produces the library, runs t_coffee and returns a Alignment in XMFA format.
 * \param aln_f The alignment file.
 * \param out_f The output file.
 */
int
prepare_low_cov(char *aln_f, char *out_f, int do_missing_blasts, char *ref_tree_f, int concatenate_species, int n_cores);



#endif /* LOW_COVERAGE_H */