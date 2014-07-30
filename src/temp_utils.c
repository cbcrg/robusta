#include "temp_utils.h"


char *
my_make_temp_file(char *template, char *function, char *file)
{
	char *temp_file_name = malloc(20 * sizeof(char*));
	sprintf(temp_file_name, "%s", template);


	int fd = mkstemp(temp_file_name);
	if (fd == -1)
	{
		fprintf(stderr, "ERROR! A temporary file could not be created.\n");
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

	temp_dir_name = mkdtemp(temp_dir_name);
	if (temp_dir_name == NULL)
	{
		fprintf(stderr, "ERROR! A temporary directory could not be created.");
		fprintf(stderr, "This error was caused by '%s' in '%s'\n", function, file);
		exit(-1);
	}
	return temp_dir_name;
}
