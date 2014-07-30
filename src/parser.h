#ifndef PARSER_H
#define PARSER_H


#include "stdio.h"
#include "malloc.h"
#include "string.h"
#include "stdlib.h"
#include "util.h"

#include "classes.h"

/*! \file parser.h
 * \brief Functions to parse alignment formats.
 *
 */


/**
 * \brief reads a single block of a file in xmfa format.
 *
 * \param The file containing the alignment. The position has to be at the alignment to read. Function stops after reading the block.
 * \return The read alignment in a Alignment struct.
 *
 */
Alignment* _read_xmfa_block(FILE *xmfa_F);





/**
 * \brief Function turning a set of alignments into a T-Coffee library.
 *
 * \param seq_files The sequence files.
 * \param num_seq_files The number of sequence files.
 * \param aln_files The nubmer of alignment files.
 * \param num_aln_files The number of alignments.
 * \param out_file The output file.
 * \param region_file The file containing the regions to align. If regionfile is NULL the whole sequence is used.
 */
void
make_lib(char **seq_files, int num_seq_files, char **aln_files, int num_aln_files, char *out_file/*, char *region_file*/);





/**
 * \brief Determines the format of a given alignemnt file.
 *
 * \param aln_file The alignemnt file.
 * \return a number depending on the format detected. 0 = maf, 1 = xmfa
 */
int format_agent(char *aln_file);





/**
 * \brief makes a header file with the given sequence files.
 *
 * \param seq_files The sequence files.
 * \param num_seq_files The number of sequence files.
 * \param lib_file_name The output file.
 */
void
make_lib_header(char **seq_files, int num_seq_files, char *lib_file_name);


/**
* \brief Adds the last line to a library file.
*
* \param lib_file_name The output file.
*/
void
make_lib_end(char *lib_file_name);



/** @name File parsing into library.
*
* This group contains functions to read an alignemnt file and turn it into a library.
*/
//@{


/**
* \brief Parses alignments in maf format.
*
* \param input_file_name The alignment file.
* \param lib_file_name The name of the ouput file.
* \param seq_names The names of the sequences.
* \param num_sequences The number of sequences.
* \returns 0 if succesfully else 1;
*/
int maf2lib_body(char *input_file_name, char *lib_file_name, char **seq_names, int num_sequences);




/**
* \brief Parses alignments in xmfa format.
*
* \param input_file_name The alignment file.
* \param lib_file_name The name of the ouput file.
* \param seq_names The names of the sequences.
* \param num_sequences The number of sequences.
*/
int xmfa2lib_body(char *input_file_name, char *lib_file_name,/* char **seq_names,*/ int num_sequences);


/**
* \brief Parses alignments in fasta format.
*
* \param input_file_name The alignment file.
* \param lib_file_name The name of the ouput file.
* \param seq_names The names of the sequences.
* \param num_sequences The number of sequences.
*/
void
fasta2lib_body(char *input_file_name, char *lib_file_name, char **seq_names, int num_sequencess);

void
lastz2lib_body_2(char *input_file_name, FILE *lib_F, unsigned int off_set1, unsigned int off_set2);

//@}





#endif /* PARSER_H */