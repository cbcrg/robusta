
#ifndef TEMP_UTILS_H
#define TEMP_UTILS_H

#include "unistd.h"
#include "string.h"
#include "stdlib.h"


/*! \file util.h
\brief General useful functions.

Contains functions to handle temporary files/directories.
*/


char * my_make_temp_file(char *template, char *function, char *file);
char * my_make_temp_dir(char *template, char *function, char *file);

#endif /* TEMP_UTILS_H */