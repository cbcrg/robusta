#ifndef CLASSES_H
#define CLASSES_H


#include "stdlib.h"
#include "string.h"
#include "stdio.h"

/*! \file classes.h
 * \brief Contains serveral structs and functions for Sequences, Alignments and T-Coffee like libraries.
 *
 *
 */




/**
 * \brief A structure representing a sequence.
 *
 */
typedef struct
{
	char *name;		/**< \brief Name of the sequence */
	char *chr;		/**< \brief Chromosome/scaffold....*/
	char *comment;  /**< \brief Comment */
	unsigned int number; /**< \brief Sequence number */
	unsigned int start; /**< \brief Start position of sequence in chromosome (start is 0) */
	unsigned int end;   /**< \brief Last position belonging to the sequence.*/
	int strand;			/**< \brief The strand +1/-1 */
	char *seq;			/**< \brief The sequence it self. */
	unsigned int length; /**< \brief Length of the sequence */
	unsigned int gap_less_length;
} Sequence;



void free_Sequence(Sequence *seq);




/**
 * \brief A structure representing a sequence
 */
typedef struct
{
	Sequence **seqs; /**< \brief Array of Sequence objects */
	unsigned int num_seqs; /**< \brief Number of sequences in this Alignment */
} Alignment;

void free_Alignment(Alignment *aln);
char * get_gap_free_seq(Alignment *aln, unsigned int in);

void
aln_2_lib_body(const Alignment *aln, int index1, int index2, int off_set1, int off_set2, FILE *lib_F);

//Library struct
/**
 * \brief Structure to save the matches of two sequences.
 *
 * It contains all information of a T-Coffee library "body" for two sequences.
 */
typedef struct
{
	int seq1; /**< \brief Index of sequence 1*/
	int seq2; /**< \brief Index of sequence 2*/
	int num_entries; /**< \brief Number of entries in the Library */
	unsigned int *pairs; /**< \brief The matches.

		Matches are saved one after another. First one index of seq1 the second index of seq2. */
} Library_pair;





//This structure stores a sequence the ensembl style:

void free_Library(Library_pair *lib);
void free_Library_content(Library_pair *lib);

Library_pair *aln_2_pairs(const Alignment *aln, int index1, int index2, int off_set1, int off_set2);






#endif /* CLASSES_H */