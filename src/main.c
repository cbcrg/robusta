


#include "stdio.h"
#include "string.h"
#include "signal.h"
#include <unistd.h>


#include "ms_aligners.h"
#include "pw_aligners.h"
#include "parser.h"
#include "system_caller.h"
#include "low_coverage.h"

static char *tmp_dir = NULL;
static char *cwd =NULL;


char *
readlink_malloc (const char *filename)
{
	int size = 500;

	char *buffer = (char *) malloc (size);
	int nchars = readlink (filename, buffer, size);
	if (nchars < size)
	{
		buffer[nchars] = '\0';
		return buffer;
	}
	return NULL;
}



/**
 * \brief Splits the input file into several sperate files with one sequence each.
 * \param seq_files The sequence file.
 * \param new_seq_files The new files.
 * \return The number of sequence files generated.
 */
int
preprocess_inputfiles(char *seq_file, char ***new_seq_files)
{

	char **seq_files = *new_seq_files;
	int max = 10;
	seq_files = malloc(10*sizeof(char*));
	const int LINE_LENGTH = 500;
	char line[LINE_LENGTH];
	FILE *input_F = fopen(seq_file, "r");
	FILE *output_F = NULL;
	int num_files = -1;
	int i;


	while (fgets(line, LINE_LENGTH, input_F) != NULL)
	{
		if (line[0] == '>')
		{
			if (++num_files == max)
			{
				max+=10;
				seq_files = realloc(seq_files, sizeof(char*)*max);
			}
			seq_files[num_files] = malloc(100*sizeof(char));
			i = 0;
			while (line[++i] != '\n')
			{
				if (line[i] == ' ')
					seq_files[num_files][i-1] = '_';
				else
					seq_files[num_files][i-1] = line[i];
			}
			seq_files[num_files][i-1] = '\0';

			if (output_F != NULL)
				fclose(output_F);

			output_F = fopen(seq_files[num_files], "w");
			fprintf(output_F, "%s", line);
		}
		else
			fprintf(output_F, "%s", line);

	}

	fclose(input_F);
	fclose(output_F);
	*new_seq_files = seq_files;
	return num_files+1;
}


/**
 * \brief This function does everything, from reading the files, producing the alignments, turning the alignments into a library.
 *
 * \param seq_files The sequence files.
 * \param num_seq_files The number of sequence files.
 * \param methods The methods to use.
 * \param num_methods The number of methods.
 * \param tree_file_name The name of the file tree.
 * \param out_file The output file.
 */
void
total(char **seq_files, int num_seq_files, char **methods, int num_methods, char *tree_file_name, char *out_f, int type)
{

	int i;
	char **aln_files;
	int num_aln_files = 0;
	if (type)
	{
		aln_files = call_ms_aligners(seq_files, num_seq_files, methods, num_methods,tree_file_name);
		num_aln_files=num_methods;
		make_lib(seq_files, num_seq_files, aln_files, num_methods, out_f);
	}
	else
	{
		aln_files = call_pw_aligners(seq_files, num_seq_files, methods, num_methods);
		num_aln_files=((num_seq_files * (num_seq_files-1)) / 2)*num_methods;
		make_lib(seq_files, num_seq_files, aln_files, num_aln_files, out_f);
	}


	for(i = 0; i < num_aln_files; ++i)
	{
		free(aln_files[i]);
	}
	free(aln_files);
}




/**
* \brief Prints out the help message of the program.
*/
void
print_help()
{
	printf("Program: G-Coffee\n");
	printf("-i FILES  --in FILES           A list of files containing the sequences to align.\n");
	printf("                               If only one file is given the sequences in this file are aligned with each other.\n");
	printf("-m NAMES  --m_method NAMES     A list of multiple sequence alignment programs to be used for the alignment\n");
	printf("-p NAMES  --pw_method NAMES    A list of pairwise sequence alignment programs to be used for the alignment\n");
	printf("-t FILE   --tree FILE          This tree is needed to produce alignments.\n");
	printf("-o FILE   --out FILES          File prefix to be used to save the libraries in.\n");
	printf("-n INT    --max_prog           The maximium number of processors to use.\n");
	printf("-f        --final_aln          If set T-Coffee is used to produce the final alignment.\n");
	printf("--temp_dir DIRECTORY           Specifies the temporary folder to use\n");
// 	printf("-c        --code               Sequence_names will be shortened and a file produced containing the connection\n");
	printf("-l        --low_cov            Enables low coverage alignments. Uses a XMFA file of -i assumes the first alignment to be a reference alignment and uses following pairwise alignment to produe a single alignment\n");
	printf("-b        --blast              In connection with \"low_cov\" this option uses lastz to fill missing pairwise alignments\n");
}



void
method_check(char **methods, int num_methods, int type)
{
	char *method_string;
	if (type)
		method_string = "mauve$pmauve$tba$mavid$pecan";
	else
		method_string = "pecan_pw$mavid_pw$lastz$cgaln$tblastx$pmauve_pw";

	int i;
	for (i = 0; i < num_methods; ++i)
	{
		if (strstr(method_string, methods[i]) == NULL)
		{
			fprintf(stderr, "Unknown method: %s\n", methods[i]);
			exit(1);
		}
	}
}



//  The signal handler function */
void thandler( int signal_type ) {
	if (tmp_dir != NULL)
	{
		killpg(getsid(0), SIGTERM);
		printf("Deleting temporary files\n");
		chdir(cwd);
		char command[500];
		sprintf(command, "rm -r %s", tmp_dir);
		system(command);
	}
	exit(signal_type);
}



/**
* \brief The main function of the program.
*
* Contains the parameter parsing.
* \param argc Number of input parameters.
* \param argv The parameters.
*/
int
main(int argc, char *argv[])
{
	setpgid(0, 0);
	signal( SIGINT, thandler);
	if (argc == 1)
	{
		print_help();
		return 0;
	}
	char *methods[10];
	int num_methods = 0;

	int type = 0;
	int type_counter = 0;

	int do_missing_blast = 0;
	char *seq_files[10];
	int num_seq_files = 0;
	char *use_as_temp = P_tmpdir;

	char *aln_files[10];
	int num_aln_files = 0;
	char *out_file = NULL;
	int max_num_processes = 8;
	int final_aln = 0;
	int i = 0;
	int low_cov = 0;
	char *tree_file_name = NULL;
	int nd = 0;

	int concatenate_species = 0;

	while (i < argc)
	{
		if (argv[i][0] == '-')
		{
			if ((!strcmp(argv[i],"-i")) || (!strcmp(argv[i], "--in")))
			{
				char *seq_str = argv[++i];
				char *seq = strtok (seq_str, ", ");
				while (seq != NULL)
				{
					seq_files[num_seq_files] = seq;
					++num_seq_files;
					seq = strtok (NULL, " ,");
				}
			}
			else if ((!strcmp(argv[i],"-m")) || (!strcmp(argv[i], "--m_methods")))
			{
				type = 1;
				++ type_counter;
				char *method_str = argv[++i];
				char *method = strtok (method_str, ", ");
				while (method != NULL)
				{
					methods[num_methods] = method;
					++num_methods;
					method = strtok (NULL, " ,");
				}
			}
			else if ((!strcmp(argv[i],"-p")) || (!strcmp(argv[i], "--pw_methods")))
			{
				type = 0;
				++ type_counter;
				char *method_str = argv[++i];
				char *method = strtok (method_str, ", ");
				while (method != NULL)
				{
					methods[num_methods] = method;
					++num_methods;
					method = strtok (NULL, " ,");
				}
			}
			else if ((!strcmp(argv[i],"-a")) || (!strcmp(argv[i], "--aln")))
			{
				char *aln_str = argv[++i];
				char *aln = strtok (aln_str, ", ");
				while (aln != NULL)
				{
					aln_files[num_aln_files] = aln;
					++num_aln_files;
					aln = strtok (NULL, " ,");
				}
			}
			else if ((!strcmp(argv[i], "-o")) || (!strcmp(argv[i], "--out")))
			{
				out_file = argv[++i];
			}
			else if ((!strcmp(argv[i], "-t")) || (!strcmp(argv[i], "--tree")))
			{
				tree_file_name = argv[++i];
			}
			else if (!strcmp(argv[i], "--temp_dir"))
			{
				use_as_temp = argv[++i];
			}
			else if ((!strcmp(argv[i], "-n")) || (!strcmp(argv[i], "--max_prog")))
			{
				max_num_processes = atoi(argv[++i]);
			}
			else if ((!strcmp(argv[i], "-f")) || (!strcmp(argv[i], "--final_aln")))
			{
				final_aln = 1;
			}
			else if (!strcmp(argv[i], "--nd"))
			{
				nd = 1;
			}
			else if ((!strcmp(argv[i], "-l")) || (!strcmp(argv[i], "--low_cov")))
			{
				low_cov = 1;
			}
			else if ((!strcmp(argv[i], "-c")) || (!strcmp(argv[i], "--concatenate_species")))
			{
				concatenate_species = 1;
			}
			else if ((!strcmp(argv[i], "-b")) || (!strcmp(argv[i], "--blast")))
			{
				do_missing_blast = 1;
			}
			else if ((!strcmp(argv[i], "-h")) || (!strcmp(argv[i], "--help")))
			{
				print_help();
				return 0;
			}
		}
		++i;
	}


	//check for neccesarry options
	if (num_seq_files == 0)
	{
		fprintf(stderr, "At least one input file is needed! Use -h to get all options.\n");
		exit(1);
	}
	if (out_file == NULL)
	{
		fprintf(stderr, "A out_file needs to given. Use -h to get all options.\n");
		exit(1);
	}


	//check methods
	if (!low_cov)
	{
		if (type_counter < 1)
		{
			fprintf(stderr, "At least one method needs to be given. Use -h to get all options.\n");
			exit(2);
		} else if (type_counter > 1)
		{
			fprintf(stderr, "You can only use pairwise methods or multiple methods, but not both at the same time\n");
			exit(2);
		}
		method_check(methods, num_methods, type);
	}


	// make temporary working directory and make links to the sequence files and tree_file
	char template[2000];

	sprintf(template, "%s/g_coffee_tmp_XXXXXX", use_as_temp);
	char *temp_dir_name = my_make_temp_dir(template, "total", "main.c");
	if (nd)
		printf("Creating temporary directory: %s\n", temp_dir_name);
	tmp_dir = temp_dir_name;

	char tmp_str[FILENAME_MAX];
	cwd = getcwd(tmp_str, FILENAME_MAX);

	char command[500];

	if (tree_file_name != NULL)
	{
		sprintf(command, "t_coffee -other_pg seq_reformat -action +tree_prune -in2 %s -in %s > %s/tree_tmp", seq_files[0], tree_file_name, tmp_dir);
		printf("C:%s\n", command);
		system(command);
	}
	chdir(temp_dir_name);


	char out_f[500];
	if (out_file[0] != '/')
		sprintf(out_f, "%s/%s", cwd, out_file);
	else
		sprintf(out_f, "%s", out_file);




	char *tmp, **new_seq_files;

	new_seq_files =  make_2D_char(num_seq_files, 40);

	if (num_seq_files == 1)
	{




		tmp = strrchr(seq_files[0], '/');
		if (tmp == NULL)
			tmp = seq_files[0];
		else
			++tmp;
		if (seq_files[0][0] != '/')
			sprintf(command, "ln -s %s/%s %s", cwd, seq_files[0], tmp);
		else
			sprintf(command, "ln -s %s %s", seq_files[0], tmp);
		sprintf(new_seq_files[0], "%s", tmp);
		system(command);
		if (!low_cov)
			num_seq_files = preprocess_inputfiles(tmp, &new_seq_files);
	}
	else
	{
		for (i = 0; i < num_seq_files; ++i)
		{
			tmp = strrchr(seq_files[i], '/');
			if (tmp == NULL)
				tmp = seq_files[i];
			else
				++tmp;
			if (seq_files[i][0] != '/')
				sprintf(command, "ln -s %s/%s %s", cwd, seq_files[i], tmp);
			else
				sprintf(command, "ln -s %s %s", seq_files[i], tmp);

			sprintf(new_seq_files[i], "%s", tmp);
			system(command);
		}
	}

	set_max_num_p(max_num_processes);
// 	char *new_tree_file_name = NULL;
// 	if (tree_file_name != NULL)
// 	{
// 		new_tree_file_name = strrchr(tree_file_name, '/');
// 		if (new_tree_file_name == NULL)
// 			new_tree_file_name = tree_file_name;
// 		else
// 			++new_tree_file_name;

// 		printf("%s\n", new_tree_file_name);
// 		if (tree_file_name[0] != '/')
// 			sprintf(command, "ln -s %s/%s %s", cwd, tree_file_name, new_tree_file_name);
// 		else
// 			sprintf(command, "ln -s %s %s", tree_file_name, new_tree_file_name);
// 		system(command);

// 	}


	int status = 0;
	if (low_cov)
	{
		status = prepare_low_cov(new_seq_files[0], out_f, do_missing_blast, "tree_tmp", concatenate_species, max_num_processes);
	}
	else
	{
		if (num_aln_files > 0)
		{
			make_lib(new_seq_files, num_seq_files, aln_files, num_aln_files, out_f);
			exit(0);
		}
		total(new_seq_files, num_seq_files, methods, num_methods,"tree_tmp", out_f, type);
	}

	chdir(cwd);
	del_2D_char(new_seq_files, num_seq_files);

	if (!nd)
	{
		sprintf(command, "rm -r %s", temp_dir_name);
		if (system(command) == 0)
			tmp_dir = NULL;
	}
	free(temp_dir_name);

	if (final_aln)
	{
		sprintf(command, "t_coffee -quiet -n_core %i -lib %s -outfile %s.fa -output=fasta_aln", max_num_processes, out_file, out_file);
		status = system(command);
	}

	return status;
}
