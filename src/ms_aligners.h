

#ifndef MS_ALIGNERS_H
#define MS_ALIGNERS_H

#include "unistd.h"
#include "string.h"
#include "stdlib.h"
#include "system_caller.h"
#include "ctype.h"

#ifdef PARALLEL
#include <omp.h>
#endif /* PARALLEL */



#include "util.h"


/*! \file ms_aligners.h
\brief Header file for functions calling different multiple genome aligners.
*/



/**
 * \brief Function to call the different multiple sequence alignment programs.
 *
 * This function produces alignments with the methods specified in \a methods. If this function is compiled in parallel mode the different methods are executed at the same time.
 * \param seq_files The list of sequnce files.
 * \param num_seq_files The number of sequence files.
 * \param methods The alignment program to use.
 * \param num_methods The number of different alignment programs.
 * \param tree_file_name The name of the file containing the tree.
 * \return An array which contains at the first position the name of the temporary directory and at the following positions the names of the alignment files produced.
 **/
char **
call_ms_aligners(char **seq_files, int num_seq_files, char **methods, int num_methods, char *tree_file_name);


/** @name Single multiple genome aligner calls
*
* This group contains functions to make single calls to whole genome alignment programs.
*/
//@{


/**
* \brief Produces a whole genome alignment using progressive Mauve.
*
* \param seq_files The names of the sequence files.
* \param num_seq_files The number of sequence files.
* \param out_file_name The name of the output file.
*/
void progressive_mauve_alignment(char **seq_files, int num_seq_files, char *out_file_name);


/**
* \brief Produces a whole genome alignment using Mauve.
*
* \param seq_files The names of the sequence files.
* \param num_seq_files The number of sequence files.
* \param out_file_name The name of the output file.
*/
void mauve_alignment(char **seq_files, int num_seq_files, char *out_file_name);


/**
* \brief Produces a genome alignment using the tba alignment program.
*
* \param tree_file_name A tree in tba format (simplified newick format).
* \param out_file_name The name of the output file.
*/
void
tba_alignment(char * tree_file_name, char *out_file_name);


/**
* \brief Produces a genome alignment using the pecan alignment program.
*
* \param tree_file_name A tree in tba format (simplified newick format).
* \param out_file_name The name of the output file.
*/
void
pecan_alignment(char * tree_file_name, char *out_file_name);



//@}



#endif /* MS_ALIGNERS_H */