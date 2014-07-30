#ifndef PW_ALIGNERS_H
#define PW_ALIGNERS_H

#include "string.h"
#include "stdlib.h"
#include "util.h"
#include "system_caller.h"
#include "stdio.h"
#include "ctype.h"

/*! \file pw_aligners.h
\brief Header file for functions calling different multiple genome aligners.
*/



/**
* \brief Function to call the different pairwise sequence alignment programs.
*
* This function produces alignments with the methods specified in \a methods. If this function is compiled in parallel mode the different methods are executed at the same time.
* \param seq_files The list of sequnce files.
* \param num_seq_files The number of sequence files.
* \param methods The alignment program to use.
* \param num_methods The number of different alignment programs.
* \return An array which contains at the first position the name of the temporary directory and at the following positions the names of the alignment files produced.
**/
char **
call_pw_aligners(char **seq_files, int num_seq_files, char **methods, int num_methods);

/*
* @name Single pairwise aligner calls
*
* This group contains functions to make single calls to pairwise whole genome alignment programs.
*



**
* \brief Produces a whole genome alignment using lastz.
*
* \param seq_files The names of the sequence files.
* \param num_seq_files The number of sequence files.
* \param out_file_name The name of the output file.
*
void lastz(char *seq_file1, char *seq_file2, char *out_file_name);



//@}*/

void
lastz_pw(char *file1, char *file2, char *aln_file, int threshold);







#endif /* PW_ALIGNERS_H */