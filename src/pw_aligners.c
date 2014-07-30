#include "pw_aligners.h"




//Format_numbers:
// cgaln = 1;
// lastz = 2;
// xmfa = 3;
// fasta = 4;
// mavid = 5;
// tblastx = 6;
void
msa_cgaln(char **seq_files, int num_seq_files, char **aln_files)
{
	int i, j, k = -1;
	char **command = malloc(7*sizeof(char*));
	command[0] = "maketable";

	command[2] = NULL;
	for (i = 0; i < num_seq_files; ++i)
	{
		command[1] = seq_files[i];
		my_system_no_wait_redirect(command,"/dev/null", "/dev/null");
	};


	wait_for_all();
	command[0] = "Cgaln";
	command[1] = "-Ens";
	command[2] = "-ia";
	command[5] = NULL;

	for (i = 0; i < num_seq_files; ++i)
	{
		for (j = i+1; j < num_seq_files; ++j)
		{
			aln_files[++k] = malloc(20 * sizeof(char));
			sprintf(aln_files[k], "1_%i_%i", i, j);
			command[3] = seq_files[i];
			command[4] = seq_files[j];
			my_system_no_wait_redirect(command, aln_files[k], "/dev/null");
		}
	}
	free(command);
}



void
msa_lastz(char **seq_files, int num_seq_files, char **aln_files)
{
	int i, j, k = -1;
	char **command = malloc(6*sizeof(char*));
	char* output = malloc(100 * sizeof(char));
	command[0] = "lastz";
	command[1] = "--strand=plus";
	command[2] = "--chain";
	command[5] = output;
	command[6] = NULL;
	for (i = 0; i < num_seq_files; ++i)
	{
		command[3] = seq_files[i];
		for (j = i+1; j < num_seq_files; ++j)
		{
			aln_files[++k] = malloc(20 * sizeof(char));
			sprintf(aln_files[k], "2_%i_%i", i, j);
			command[4] = seq_files[j];
			sprintf(output,"--output=%s", aln_files[k]);
			my_system_no_wait(command);
		}
	}
	free(output);
	free(command);
}



void
msa_mauve_prog(char **seq_files, int num_seq_files, char **aln_files)
{
	char **command = malloc(4*sizeof(char*));
	char command_part[500];
	char tmp_name[20];

	command[0] = "sh";
	command[1] = "-c";
	command[2] = command_part;
	command[3] = NULL;
	int i, j, k = -1;
	for (i = 0; i < num_seq_files; ++i)
	{
		for (j = i+1; j < num_seq_files; ++j)
		{
			sprintf(tmp_name, "%i_%i_mauve_pw", i, j);
			sprintf(command_part, "cat %s %s > %s", seq_files[i], seq_files[j], tmp_name);
			system(command_part);

			aln_files[++k] = malloc(20 * sizeof(char));
			sprintf(aln_files[k], "3_%i_%i", i, j);
			sprintf(command_part, "progressiveMauve --disable-backbone --collinear --output %s %s", aln_files[k], tmp_name);
			my_system_no_wait_redirect(command, "/dev/null", "/dev/null");
		}
	}
	free(command);
}




void
msa_pecan_pw(char **seq_files, int num_seq_files, char **aln_files)
{
	char **command = malloc(9*sizeof(char*));
	char command_part[1001];

	command[0] = "sh";
	command[1] = "-c";
	command[2] = command_part;
	command[3] = NULL;
	char *pecan_path = getenv ( "PECAN_PATH" );
	int i, j, k = -1;
	for (i = 0; i < num_seq_files; ++i)
	{
		for (j = i+1; j < num_seq_files; ++j)
		{
			aln_files[++k] = malloc(20 * sizeof(char));
			sprintf(aln_files[k], "4_%i_%i", i, j);
			sprintf(command_part, "java -cp %s bp.pecan.Pecan -E '(A , B);' -F %s %s -G %s ", pecan_path, seq_files[i], seq_files[j], aln_files[k]);
			my_system_no_wait(command);
		}
	}
	free(command);

}

void
soft2hardmasked(char *soft_masked_f, char *hard_masked_f)
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
msa_mavid_pw(char **seq_files, int num_seq_files, char **aln_files)
{
	char *command = malloc(1000*sizeof(char));

	FILE *tree_F;
	int i, j, k = -1;
	for (i = 0; i < num_seq_files; ++i)
	{
		for (j = i+1; j < num_seq_files; ++j)
		{
			tree_F = fopen("tmp_mavid_tree","w");
			fprintf(tree_F, "(%s:0.1, %s:0.1):0.1;", seq_files[i], seq_files[j]);
			fclose(tree_F);
			sprintf(command, "cat %s %s > pair_file",  seq_files[i], seq_files[j]);
			system(command);
// 			sprintf(command, "cp pair_file pair_file.masked");
// 			system(command);
			soft2hardmasked("pair_file", "pair_file.masked");
			aln_files[++k] = malloc(20 * sizeof(char));
			sprintf(aln_files[k], "5_%i_%i", i, j);

			sprintf(command, "mavid tmp_mavid_tree pair_file 2>/dev/null >/dev/null");
			system(command);
			sprintf(command, "mv mavid.mfa %s", aln_files[k]);
			system(command);
		}
	}
	free(command);
}




void
msa_tblastx_pw(char **seq_files, int num_seq_files, char **aln_files)
{
	char command[100];
	int i;

	for (i = 0; i < num_seq_files; ++i)
	{
		sprintf(command, "formatdb -i %s -p F", seq_files[i]);
		system(command);
	}

	char **commands = malloc(10 * sizeof(char*));
	char out_file[20], in_file[20], db_file[20];
	commands[0] = "blastall";
	commands[1] = "-ptblastx";
	commands[2] = &(in_file[0]);
	commands[3] = &(db_file[0]);
	commands[4] = "-m7";
	commands[5] = "-S1";
	commands[6] = "-e5";
	commands[7] = &(out_file[0]);
	commands[8] = NULL;


	//run blast
	int j, k = -1;
	for (i = 0; i < num_seq_files; ++i)
	{
		sprintf(db_file, "-d%s", seq_files[i]);
		for (j = i+1; j < num_seq_files; ++j)
		{
			sprintf(in_file, "-i%s", seq_files[j]);
			aln_files[++k] = malloc(20*sizeof(char));
			sprintf(aln_files[k], "6_%i_%i_all_blast.out", j, i);
			sprintf(out_file, "-o%s", aln_files[k]);
			my_system_no_wait(commands);
		}
	}
	free(commands);
}




int
check_cgaln_aln(char *aln_f)
{
	FILE *aln_F = fopen(aln_f, "r");
	fseek(aln_F, 0L, SEEK_END);
	unsigned int sz = ftell(aln_F);
	fclose(aln_F);
	if (sz == 0)
		return 0;
	else
		return 1;
}



void
cgaln_pw(char *file1, char *file2, char *aln_file)
{
	char **command = malloc(7*sizeof(char*));
	command[0] = "Cgaln";
	command[1] = "-Ens";
	command[2] = "-ia";
	command[3] = "-X100";
	command[4] = file1;
	command[5] = file2;
	command[6] = NULL;
	my_system_no_wait_redirect(command, aln_file, "/dev/null");
	free(command);
}




/**
 * Checks if the library contains match information.
 *
 * \param aln_f Alignemnt file to check.
 * return 1 if at least one match information is found.
 */
int
check_lastz_aln(char *aln_f)
{
	FILE *aln_F = fopen(aln_f, "r");
	const int LINE_LENGTH = 501;
	char line[LINE_LENGTH];

	int found = 0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == 's')
		{
			found = 1;

			break;
		}
	}
	fclose(aln_F);
	return found;
}



void
lastz_pw(char *file1, char *file2, char *aln_file, int threshold)
{
	char **command = malloc(8*sizeof(char*));
	char* output = malloc(100 * sizeof(char));
	char* thres_param = malloc(20*sizeof(char));
	sprintf(thres_param,"--hspthresh=%i", threshold);
	command[0] = "lastz";
	command[1] = "--strand=plus";
	command[2] = "--chain";
	command[3] = thres_param;
	command[4] = file1;
	command[5] = file2;
	command[6] = output;
	command[7] = NULL;
	sprintf(output,"--output=%s", aln_file);
	my_system_no_wait(command);
	free(thres_param);
	free(output);
	free(command);
}




char **
call_pw_aligners(char **seq_files, int num_seq_files, char **methods, int num_methods)
{
	// 	stores the output files
	int num_files = (num_seq_files * (num_seq_files-1)) /2;
	int num_aln_files =  num_methods * num_files;
	char **aln_files = malloc( (num_aln_files+1) * sizeof(char *));
	int i;
	int problematic_method = 0;
	for (i = 0; i < num_methods; ++i)
	{
		if (!strcmp(methods[i], "lastz"))
		{
			problematic_method = 1;
			msa_lastz(seq_files, num_seq_files, &aln_files[i*num_files]);

		}
		else if (!strcmp(methods[i], "cgaln"))
		{
			problematic_method = 1;
			msa_cgaln(seq_files, num_seq_files, &aln_files[i*num_files]);
		}
		else if (!strcmp(methods[i], "pecan_pw"))
		{
			msa_pecan_pw(seq_files, num_seq_files, &aln_files[i*num_files]);
		}
		else if (!strcmp(methods[i], "mavid_pw"))
		{
			msa_mavid_pw(seq_files, num_seq_files, &aln_files[i*num_files]);
		}
		else if (!strcmp(methods[i], "pmauve_pw"))
		{
			msa_mauve_prog(seq_files, num_seq_files, &aln_files[i*num_files]);
		}
		else if (!strcmp(methods[i], "tblastx"))
		{
			msa_tblastx_pw(seq_files, num_seq_files, &aln_files[i*num_files]);
		}
	}

	if ((num_methods != 1) || (strcmp(methods[0], "mavid_pw")))
		if (wait_for_all() >0)
		{
			fprintf(stderr,"Error occured\n");
			exit(EXIT_FAILURE);
		}


	//check_lastz_aln
	if (problematic_method)
	{
		int threshold = 1000;
		char tmp[50];
		unsigned int f1, f2;
		int found_empty=1;
		int cgaln_done = 0;
		while (found_empty)
		{
			found_empty = 0;
			for (i = 0; i < num_aln_files; ++i)
			{
				if (aln_files[i] != NULL)
				{
					if ((aln_files[i][0] == '2') && (!check_lastz_aln(aln_files[i])))
					{
						found_empty = 1;
						strcpy(tmp, aln_files[i]);
						strtok(tmp, "_");
						f1 = atoi(strtok(NULL, "_"));
						f2 = atoi(strtok(NULL, "_"));
						lastz_pw(seq_files[f1], seq_files[f2], aln_files[i], threshold);
					}
					else if ((!cgaln_done) && (aln_files[i][0] == '1') && (!check_cgaln_aln(aln_files[i])))
					{
						found_empty = 1;
						strcpy(tmp, aln_files[i]);
						strtok(tmp, "_");
						f1 = atoi(strtok(NULL, "_"));
						f2 = atoi(strtok(NULL, "_"));
						cgaln_pw(seq_files[f1], seq_files[f2], aln_files[i]);
					}
				}
			}
			cgaln_done = 1;
			threshold /= 2;
			if (threshold < 500)
				break;
			wait_for_all();
		}
	}
	aln_files[num_aln_files] = NULL;
	return aln_files;
}


