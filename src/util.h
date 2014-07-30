
#ifndef UTIL_H
#define UTIL_H


#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "errno.h"
#include "unistd.h"

/*! \file util.h
\brief General useful functions.

Contains functions for memory management and temporary files/directories.
*/


// Temporary files/dirs functions

/** @name Temporary file/directory managment.
*
* This group contains functions which make temporary files or direcoties. If for some reason this is not possible the 'exit' function is calles and an error message displayed. The returned error code in this case is -1.
*/
//@{


/**
 * \brief Makes a temporary file.
 * \param template The template name. Has to end with XXXXXX.
 * \param function The function from which this function is called.
 * \param file The name of the file containing the function.
 */
char *
my_make_temp_file(char *template, char *function, char *file);



/**
* \brief Makes a temporary directory.
* \param template The template name. Has to end with XXXXXX.
* \param function The function from which this function is called.
* \param file The name of the file containing the function.
*/
char *
my_make_temp_dir(char *template, char *function, char *file);

//@}


// Memory Function


/** @name Memory management functions
*
* This group contains functions which deal with 2D arrays.
*/
//@{

/**
* \brief Constructs a two dimensional matrix of doubles.
* \param d1 The size of dimension 1.
* \param d2 The size of dimension 2.
* \return Pointer to the matrix.
*/
double ** make_2D_double(int d1, int d2);

/**
* \brief Frees the memory occupied by a two dimensional array if doubles.
* \param matrix A pointer to a two dimensional array.
* \param d1 Length of the first dimension.
*/
void del_2D_double(double **matrix, int d1);



/**
* \brief Constructs a two dimensional matrix of integers.
* \param d1 The size of dimension 1.
* \param d2 The size of dimension 2.
* \return Pointer to the matrix.
*/
int** make_2D_int(int d1, int d2);

/**
* \brief Frees the memory occupied by a two dimensional array of integers.
* \param matrix A pointer to a two dimensional array.
* \param d1 Length of the first dimension.
*/
void del_2D_int(int **matrix, int d1);




/**
* \brief Constructs a two dimensional matrix of chars.
* \param d1 The size of dimension 1.
* \param d2 The size of dimension 2.
* \return Pointer to the matrix.
*/
char** make_2D_char(int d1, int d2);

/**
* \brief Frees the memory occupied by a two dimensional array of chars.
* \param matrix A pointer to a two dimensional array.
* \param d1 Length of the first dimension.
*/
void del_2D_char(char **matrix, int d1);




typedef struct
{
	unsigned int v1;
	unsigned int v2;
} Pair;





//@}

#endif /* UTIL_H */
