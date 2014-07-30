#include "ms_aligners.h"



// RENAMING OF OUTFORMAT NEEDED TO BE COMPATIBLE WITH PAIRWISE SOLUTION!
void
progressive_mauve_alignment(char **seq_files, int num_seq_files, char *out_file_name)
{
	char command[501];
	sprintf(command, "progressiveMauve --collinear --output %s", out_file_name);
	int i;
	for (i = 0; i < num_seq_files; ++i)
	{
		sprintf(command, "%s %s", command, seq_files[i]);
	}
	sprintf(command, "%s 2> /dev/null >/dev/null", command);
// 	printf("%s\n", command);
	system(command);
}


void
mauve_alignment(char **seq_files, int num_seq_files, char *out_file_name)
{

	char command[501];
	sprintf(command, "mauveAligner --output-alignment=%s", out_file_name);
	int i;
	char *tmp_f_name;
	for (i = 0; i < num_seq_files; ++i)
	{
		tmp_f_name = my_make_temp_file("tmp_XXXXXX", "call aligners", "aligners.c");
		sprintf(command, "%s %s %s", command, seq_files[i], tmp_f_name);
	}
	sprintf(command, "%s 2> /dev/null >/dev/null", command);
	system(command);
}


void
concatenate_files(char **seq_files, int num_seq_files, char *out_file_name)
{
	FILE *concat_out_f = fopen(out_file_name, "w");
	FILE *in_f;
	int i;
	const int LINE_LENGTH=201;
	char line[LINE_LENGTH];
	for (i = 0; i < num_seq_files; ++i)
	{
		in_f = fopen(seq_files[i], "r");
		while (fgets(line, LINE_LENGTH, in_f) != NULL)
		{
			if (line[0] == '>')
				fprintf(concat_out_f, ">%s\n", seq_files[i]);
			else
				fprintf(concat_out_f, line);
		}
		fclose(in_f);
	}
	fclose(concat_out_f);
}


void
soft2hardmasked2(char *soft_masked_f, char *hard_masked_f)
{
	// 	printf("HERE\n");
	FILE *soft_masked_F = fopen(soft_masked_f, "r");
	FILE *hard_masked_F = fopen(hard_masked_f, "w");
	const int LINE_LENGTH = 801;
	char line[LINE_LENGTH];
	unsigned int i;
	while (fgets(line, LINE_LENGTH, soft_masked_F) != NULL)
	{
		if (line[0] != '>')
		{
			i = 0;
			while ((line[i] != '\n')  && (line[i] != '\0'))
			{
				if (islower(line[i]))
					line[i] = 'N';
				++i;
			}
		}
		fprintf(hard_masked_F, "%s", line);
	}
	fclose(soft_masked_F);
	fclose(hard_masked_F);
}


void
mavid_alignment(char **seq_files, int num_seq_files, char *tree_file_name, char *aln_file)
{
	concatenate_files(seq_files, num_seq_files, "all_fasta_mavid");

 	char command[501];
	sprintf(command, "all_fasta_mavid.masked");
	soft2hardmasked2("all_fasta_mavid", command);
// 	sprintf(command, "RepeatMasker all_fasta_mavid >/dev/null 2>/dev/null");
// 	system(command);

	sprintf(command, "mavid %s all_fasta_mavid >/dev/null 2>/dev/null", tree_file_name);
	system(command);

	sprintf(command, "scp mavid.mfa %s", aln_file);
	system(command);
}



void
tba_alignment(char * tree_file_name, char *out_file_name)
{
	//read tree
	FILE *tree_f = fopen(tree_file_name, "r");
	const int LINE_LENGTH = 501;
	char tree[LINE_LENGTH];
	fgets(tree, LINE_LENGTH, tree_f);
	fclose(tree_f);

	int i = -1;
	while ((tree[++i] != '\n') && (tree[i] != '\0'))
	{
		if ((tree[i] == ',') || (tree[i] == ';'))
		{
			tree[i] = ' ';
		}
	}

	char tree_string[200];
	char *command[4];

	command[0] = "sh"; // call of sh needed!!! else problems with all_bz
	command[1] = "-c";
	sprintf(tree_string, "all_bz \"%s\"", tree);
	printf("%s\n", tree_string);
	command[2] = tree_string;
	command[3] = NULL;
	my_system_redirect(command, "/dev/null", "/dev/null");

	sprintf(tree_string, "tba \"%s\" *.*.maf %s", tree, out_file_name); // call of sh needed!!! else problems with tba
	command[0] = "sh";
	command[1] = "-c";
	command[2] = tree_string;
	command[3] = NULL;
	my_system_redirect(command, "/dev/null", "/dev/null");
}



/**
 * \brief Puts the the filenames into the order in which they appear in the tree.
 *
 * \param tree The tree in newick format (or a simplified version without ',' and ';'.
 * \return The filenames in a string.
 */
char *
pecan_file_order(char *tree)
{
	char *file_string = malloc(1000 * sizeof(char));
	file_string[0] = '\0';
	int i = 0;
	while ((tree[i] != '\n') && (tree[i] != '\0'))
	{
		if (tree[i] == ':')
			while ((tree[i] != ')') && (tree[i] != ',') )
			{
				if (tree[i] == '\0')
					break;
				tree[i] = ' ';
				++i;
			}
		if (tree[i] == '\0')
			break;
		++i;
	}
	char *delims = " )(,;\n";
	char *part = strtok (tree, delims);
	while (part != NULL)
	{
		sprintf(file_string, "%s %s", file_string, part);
		part = strtok (NULL, delims);
	}

	return file_string;
}


void
pecan_alignment(char * tree_file_name, char *out_file_name)
{
	FILE *tree_f = fopen(tree_file_name, "r");
	const int LINE_LENGTH = 1001;
	char tree[LINE_LENGTH];
	fgets(tree, LINE_LENGTH, tree_f);
	fclose(tree_f);

	char tree2[LINE_LENGTH];
	sprintf(tree2, "%s\n", tree);

//  	printf("T: %s\n", tree2);
	char *file_order = pecan_file_order(tree2);

	char *pecan_path = getenv ( "PECAN_PATH" );
	char command[10001];
	sprintf(command, "java -cp %s bp.pecan.Pecan -E '%s' -F %s -G %s\n", pecan_path, tree, file_order, out_file_name);
// 	printf("C: %s\n", command);
	system(command);
}




char**
call_ms_aligners(char **seq_files, int num_seq_files, char **methods, int num_methods, char *tree_file_name)
{

	// 	stores the output files
	char **aln_files = malloc((num_methods+1) * sizeof(char *));
	int i;

	//	calling the different alignment methods with the different alignment programs.
	#ifdef PARALLEL
	#pragma omp parallel shared(aln_files, seq_files, num_seq_files, tree_file_name, num_methods, methods) private(i)
	{
		#pragma omp for schedule(dynamic) nowait
	#endif

		for (i = 0; i < num_methods; ++i)
		{
			if (!strcmp(methods[i], "pmauve"))
			{
				aln_files[i] = malloc(20*sizeof(char));
				sprintf(aln_files[i], "3_msa_pmauve.xmfa");
				progressive_mauve_alignment(seq_files, num_seq_files, aln_files[i]);
			}
			else if (!strcmp(methods[i], "mauve"))
			{
				aln_files[i] = malloc(20*sizeof(char));
				sprintf(aln_files[i], "3_msa_mauve.xmfa");
				mauve_alignment(seq_files, num_seq_files, aln_files[i]);
			}
			else if (!strcmp(methods[i], "tba"))
			{
				aln_files[i] = malloc(20*sizeof(char));
				sprintf(aln_files[i], "7_msa_tba.maf");
				tba_alignment(tree_file_name, aln_files[i]);
			}
			else if (!strcmp(methods[i], "pecan"))
			{
				aln_files[i] = malloc(20*sizeof(char));
				sprintf(aln_files[i], "4_msa_pecan.fa");
				pecan_alignment(tree_file_name, aln_files[i]);
			}
			else if (!strcmp(methods[i], "mavid"))
			{
  				aln_files[i] = malloc(20*sizeof(char));
				sprintf(aln_files[i], "5_msa_mavid.fa");
				mavid_alignment(seq_files, num_seq_files, tree_file_name, aln_files[i]);
			}
		}
	#ifdef PARALLEL
	}
	#endif

	return aln_files;
}