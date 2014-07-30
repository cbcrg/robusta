#include "util.h"




// Temporary files/dirs functions


char *
my_make_temp_file(char *template, char *function, char *file)
{
// 	char *tmp_dir = P_tmpdir;
	char *temp_file_name = malloc(500 * sizeof(char*));
	sprintf(temp_file_name, "%s", template);


	int fd = mkstemp(temp_file_name);
	if (fd == -1)
	{
		int errsv = errno;
		fprintf(stderr, "ERROR! A temporary file could not be created: %s\n",strerror(errsv));
		fprintf(stderr, "This error was caused by '%s' in '%s'\n", function, file);
		exit(-1);
	}
	else
	{
		FILE *f = fdopen(fd, "w+");
		fclose(f);
	}
	return temp_file_name;
}



char *
my_make_temp_dir(char *template, char *function, char *file)
{

	char *temp_dir_name = malloc(20 * sizeof(char*));
	sprintf(temp_dir_name, "%s", template);
	char hostname[50];
	gethostname(hostname, 50);
	printf("%s %s\n",hostname, temp_dir_name);
	if ((temp_dir_name = mkdtemp(temp_dir_name))==NULL)
	{
		int errsv = errno;
		fprintf(stderr, "ERROR! A temporary directory could not be created: %s\n",strerror(errsv));
		fprintf(stderr, "This error was caused by '%s' in '%s'\n", function, file);

		exit(-1);
	}
	return temp_dir_name;
}




// Memory Functions

double **
make_2D_double(int d1, int d2)
{
	double **matrix = (double**)malloc(d1 * sizeof(double*));
	int i;
	for (i = 0; i < d1; ++i)
	{
		matrix[i] = (double*)malloc(d2 * sizeof(double));
	}
	return matrix;
}



int**
make_2D_int(int d1, int d2)
{
	int **matrix = (int**)malloc(d1 * sizeof(int*));
	int i;
	for (i = 0; i < d1; ++i)
	{
		matrix[i] = (int*)malloc(d2 * sizeof(int));
	}
	return matrix;
}


char**
make_2D_char(int d1, int d2)
{
	char **matrix = malloc(d1 * sizeof(char*));
	int i;
	for (i = 0; i < d1; ++i)
	{
		matrix[i] = malloc(d2 * sizeof(char));
	}
	return matrix;
}

void
del_2D_double(double **matrix, int d1)
{
	int i;
	for (i = 0; i < d1; ++i)
	{
		free(matrix[i]);
	}
	free(matrix);
}



void
del_2D_int(int **matrix, int d1)
{
	int i;
	for (i = 0; i < d1; ++i)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void
del_2D_char(char **matrix, int d1)
{
	int i;
	for (i = 0; i < d1; ++i)
	{
		free(matrix[i]);
	}
	free(matrix);
}


int**
make_triangle_int(int dim)
{
	int **matrix = (int**)malloc(dim * sizeof(int*));
	int i;
	for (i = 0; i < dim; ++i)
	{
		matrix[i] = (int*)calloc(dim-i , sizeof(int));
	}
	return matrix;
}

void
del_triangle_int(int **matrix, int dim)
{
	del_2D_int(matrix, dim);
}

