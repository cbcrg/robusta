#include "parser.h"



/**
* \brief Writes the matches blockwise into the library.
*
* \param aln_block The sequences.
* \param aln_block_len Length of the block.
* \param starts The starting positions of the block.
* \param num_sequences The number of sequences in the alignment (not in the block).
* \param lib_file_name The file where the matches will be appended.
*/
void
sequence_pair2lib(char **aln_block, int aln_block_len, int *starts, char *lib_file_name, int seq_num1, int seq_num2)
{

	FILE *lib_f = fopen(lib_file_name, "a");
	int k;
	int block_len = 0;
	int start_seq1 = -1, start_seq2 = -1;
	int seq_pos1, seq_pos2;

	seq_pos1 = starts[0];
	seq_pos2 = starts[1];
	fprintf(lib_f, "#%i %i\n", seq_num1, seq_num2);

	for (k = 0; k < aln_block_len; ++k)
	{
		if ((aln_block[0][k] != '-') && (aln_block[1][k] != '-'))
		{
			if (block_len == 0)
			{
				start_seq1 = seq_pos1;
				start_seq2 = seq_pos2;
			}
			++seq_pos1;
			++seq_pos2;
			++block_len;
		}
	else
		{
			if (aln_block[0][k] != '-')
				++seq_pos1;
			if (aln_block[1][k] != '-')
				++seq_pos2;
			if (block_len > 0)
			{
				fprintf(lib_f, "++ %i %i %i 100\n", block_len, start_seq1, start_seq2 );
				block_len = 0;
			}
		}
	}
	if (block_len > 0)
	{
		fprintf(lib_f, "++ %i %i %i 100\n", block_len, start_seq1, start_seq2 );
		block_len = 0;
	}

	fclose(lib_f);
}



void
make_lib_header(char **seq_files, int num_seq_files, char *lib_file_name)
{
	FILE *lib_f = fopen(lib_file_name, "w");
	FILE *seq_f;
	fprintf(lib_f, "! TC_LIB_FORMAT_01\n%i\n", num_seq_files);
	long int position_marker;
	int i, j;
	char seq_len_str[10];
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];

	unsigned long int seq_length;
	for (i = 0; i < num_seq_files; ++i)
	{
		seq_length = 0;
		fprintf(lib_f, "%s ", seq_files[i]);
		position_marker = ftell(lib_f);
		fprintf(lib_f, "           ");
		seq_f = fopen(seq_files[i], "r");

		while (fgets(line, LINE_LENGTH, seq_f) != NULL)
		{
			if (line[0] == '>')
				continue;
			else
			{
				j = 0;
				while ((line[j] != '\0') && (line[j] != '\n'))
				{
					++seq_length;
					fprintf(lib_f, "%c", line[j]);
					++j;
				}
			}
		}
		fclose(seq_f);
		fprintf(lib_f, "\n");
		fseek (lib_f , position_marker , SEEK_SET );
		sprintf(seq_len_str, "%li", seq_length);
		fputs ( seq_len_str, lib_f );
		fseek (lib_f , 0 , SEEK_END );
	}

	fclose(lib_f);
}


void
make_lib_end(char *lib_file_name)
{
	FILE *lib_f = fopen(lib_file_name, "a");
	fprintf(lib_f, "! CPU 0\n! SEQ_1_TO_N\n");
	fclose(lib_f);
}





/**
 * \brief Writes the matches blockwise into the library.
 *
 * \param aln_block The sequences.
 * \param aln_block_len Length of the block.
 * \param starts The starting positions of the block.
 * \param num_sequences The number of sequences in the alignment (not in the block).
 * \param lib_file_name The file where the matches will be appended.
 */
void
sequences2lib(char **aln_block, unsigned int aln_block_len, int *starts, int num_sequences, char *lib_file_name)
{
	FILE *lib_f = fopen(lib_file_name, "a");
	int i, j;
	unsigned int k;
	int block_len = 0;
	int start_seq1 = -1, start_seq2 = -1;
	int seq_pos1, seq_pos2;
	for ( i = 0; i < num_sequences; ++i)
	{
		if ((aln_block[i][0] == '\0') || (starts[i]<0))
			continue;
		for ( j = i+1; j < num_sequences; ++j)
		{
			if ((aln_block[j][0] == '\0')|| (starts[j]<0))
				continue;
			fprintf(lib_f, "#%i %i\n", i, j);
			seq_pos1 = starts[i];
			seq_pos2 = starts[j];

			for (k = 0; k < aln_block_len; ++k)
			{
				if ((aln_block[i][k] != '-') && (aln_block[j][k] != '-'))
				{
					if (block_len == 0)
					{
						start_seq1 = seq_pos1;
						start_seq2 = seq_pos2;
					}
					++seq_pos1;
					++seq_pos2;
					++block_len;
				}
				else
				{
					if (aln_block[i][k] != '-')
						++seq_pos1;
					if (aln_block[j][k] != '-')
						++seq_pos2;
					if (block_len > 0)
					{
						fprintf(lib_f, "++ %i %i %i 100\n", block_len, start_seq1, start_seq2 );
						block_len = 0;
					}
				}
			}
			if (block_len > 0)
			{
				fprintf(lib_f, "++ %i %i %i 100\n", block_len, start_seq1, start_seq2 );
				block_len = 0;
			}
		}
	}
	fclose(lib_f);
}

/*
int
format_agent(char *aln_file)
{

	FILE *aln_f = fopen(aln_file,"r");

	char line[100];
	fgets(line, 100, aln_f);
	fclose(aln_f);

	if (strstr(line, "##maf") != NULL)
	{
		return 0;
	}
	else if (strstr(line, "#FormatVersion Mauve1") != NULL)
	{
		return 1;
	}
	else if (strstr(line, "#i=") != NULL)
	{
		return 2;
	}
	if (strstr(line, "#:lav") != NULL)
	{
		return 3;
	}
	else
	{
		return 4;
	}
}*/

int
get_seq_num2(char **seq_names, int num_sequences, char *seq_name_line)
{
	int i;
	unsigned int str_length=strlen(seq_name_line) -1;
	if (seq_name_line[str_length] == '\n')
		seq_name_line[str_length] = '\0';

	for (i = 0; i < num_sequences; ++i)
	{
// 		printf("%s-%sE\n", seq_names[i], seq_name_line);
		if (strcmp(seq_name_line, seq_names[i]) == 0 )
			break;
	}
	return i+1;
}



// int
// get_seq_num(char **seq_names, int num_sequences, char *seq_name_line)
// {
// 	int i;
// 	for (i = 0; i < num_sequences; ++i)
// 	{
// 		if (strstr(seq_name_line, seq_names[i]) != NULL )
// 			break;
// 	}
// 	return i+1;
// }



int
maf2lib_body(char *input_file_name, char *lib_file_name, char **seq_names, int num_sequences)
{

	int i;
	++num_sequences;
	int *starts = malloc(num_sequences * sizeof(int));
	int current_size = 1000;
	char **aln_block = make_2D_char(num_sequences, current_size);
	for (i = 0; i < num_sequences; ++i)
		aln_block[i][0] = '\0';

	//skip header if exists
	FILE *aln_f = fopen(input_file_name, "r");
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];

	long int pos = ftell(aln_f);
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] == 'a')
			break;
		pos = ftell(aln_f);
	}
	fseek ( aln_f , pos , SEEK_SET );


	//Parse body
	char *seq_name;
	int j;
	char *seq_start;
	unsigned int seq_num;
	int block_length = -1;
	char strand;
	int found = 1;
	int skip= 0;
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] == 's') // signals start of a new sequence
		{
			block_length = -1;
			seq_name = strtok ( &line[2]," ");
			seq_num = get_seq_num2(seq_names, num_sequences-1, seq_name);
			starts[seq_num] = atoi(strtok ( NULL," "))+1;

			strtok ( NULL, " ");
			strand = strtok ( NULL, " ")[0];
			if (strand == '-')
				skip=1;
			strtok ( NULL, " ");
			seq_start = strtok ( NULL, " ");

			j = 0;
			while ((seq_start[j] != '\n') && (seq_start[j] != '\0'))
			{
				aln_block[seq_num][++block_length] = seq_start[j];
				++j;
			}


			while (fgets(line, LINE_LENGTH, aln_f) != NULL)
			{
				if (current_size - LINE_LENGTH <= block_length)
				{
					current_size += 1000;
					for ( i = 0; i < num_sequences; ++i)
					{
						aln_block[i] = realloc(aln_block[i], current_size * sizeof(char));
					}
				}
				j= 0;
				while ((line[j] != '\n') && (line[j] != '\0'))
				{
					aln_block[seq_num][++block_length] = line[j];
					++j;
				}
				if (line[j] == '\n')
				{
					break;
				}
			}
		}
		else if ((line[0] == ' ') || (line[0] == '\n'))  // signals end of a block
		{
			found = 0;
			if (skip==0)
			{
				sequences2lib(aln_block, block_length+1, starts, num_sequences, lib_file_name);
			}
			skip =0;
			for (i = 0; i < num_sequences; ++i)
				aln_block[i][0] = '\0';
		}

	}

	fclose(aln_f);
	free(starts);
	del_2D_char(aln_block, num_sequences);
	if (block_length > 1)
		return 0;
	else return 1;

}




int
xmfa2lib_body( char *input_file_name, char *lib_file_name, int num_sequences)
{
// 	printf("READING\n");
	int *files= malloc((num_sequences+1)*sizeof(int));
	char tmp_str[200];
	strcpy(&tmp_str[0], input_file_name);
	strtok(&tmp_str[0],"_");
	char *part;
	int f = 0;
	int i;
	part = strtok(NULL, "_");
	if (part[0] != 'm')
	{
		files[++f]= atoi(part)+1;
		while ((part = strtok(NULL, "_")) != NULL)
		{
			files[++f]= atoi(part)+1;
		}
	}
	else
	{
		for (i = 1; i <=num_sequences; ++i)
		{
			files[++f] = i;
		}
	}

	int found = 1;
	++num_sequences;
	int *starts = malloc(num_sequences * sizeof(int));
	int current_size = 1000;
	char **aln_block = make_2D_char(num_sequences, current_size);
	for (i = 0; i < num_sequences; ++i)
		aln_block[i][0] = '\0';


	FILE *aln_F = fopen(input_file_name, "r");
	if (aln_F == NULL)
	{
		fprintf(stderr, "File %s could not be opened", input_file_name);
	}
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];

	long int pos = ftell(aln_F);
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '>')
			break;
		pos = ftell(aln_F);
	}
	fseek ( aln_F , pos , SEEK_SET );


	//Parse body
	int j;
	int seq_start;
	unsigned int seq_num = -1;
	int block_length = -1;
	char *curr_seq_name = NULL;
	char strand;
// 	int skip = 0;
	while (fgets(line, LINE_LENGTH, aln_F) != NULL)
	{
		if (line[0] == '>') // signals start of a new sequence
		{
			unsigned int tmp = 0;
			while (line[++tmp] == ' ');
			seq_num=files[atoi(strtok ( &line[tmp]," -:"))];
			seq_start = atoi(strtok(NULL, " ,-"));
			strtok ( NULL, " ");
			strand = strtok ( NULL ," ")[0];
			curr_seq_name = strtok ( NULL," ");

			if (strand == '-')
				starts[seq_num] = -1;
			else
				starts[seq_num] = seq_start;
			block_length = -1;
		}
		else if (line[0] == '=')  // signals end of a block
		{
			found = 0;
			sequences2lib(aln_block, block_length+1, starts, num_sequences, lib_file_name);
			for (i = 0; i < num_sequences; ++i)
				aln_block[i][0] = '\0';
		}
		else //the sequence
		{
			if (current_size - LINE_LENGTH <= block_length)
			{
				current_size += 1000;
				for ( i = 0; i < num_sequences; ++i)
				{
					aln_block[i] = realloc(aln_block[i], current_size * sizeof(char));
				}
			}
			j = 0;
			while ((line[j] != '\n') && (line[j] != '\0'))
			{
				aln_block[seq_num][++block_length] = line[j];
				++j;
			}
		}

	}
	free(files);
	fclose(aln_F);
	free(starts);
	del_2D_char(aln_block, num_sequences);
	return found;
}




void
fasta2lib_body(char *input_file_name, char *lib_file_name, char **seq_names, int num_sequences)
{
	int i;
	++num_sequences;
	int *starts = malloc((num_sequences) * sizeof(int));
	int current_size = 1000;
	char **aln_block = make_2D_char(num_sequences, current_size);
	for (i = 0; i < num_sequences; ++i)
		aln_block[i][0] = '\0';


	FILE *aln_f = fopen(input_file_name, "r");
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];

	//Parse body
	int j;
	unsigned int seq_num = -1;
	int block_length = -1;
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] == '>') // signals start of a new sequence
		{
			seq_num = get_seq_num2(seq_names, num_sequences, &line[1]);
			starts[seq_num] = 1;

			block_length = -1;
		}
		else //the sequence
		{
			if (current_size - LINE_LENGTH <= block_length)
			{
				current_size += 1000;
				for ( i = 1; i < num_sequences; ++i)
				{
					aln_block[i] = realloc(aln_block[i], current_size * sizeof(char));
				}
			}
			j = 0;
			while ((line[j] != '\n') && (line[j] != '\0'))
			{
				aln_block[seq_num][++block_length] = line[j];
				++j;
			}
		}
	}

	sequences2lib(aln_block, block_length+1, starts, num_sequences, lib_file_name);
	fclose(aln_f);
	free(starts);
	del_2D_char(aln_block, num_sequences);

}


void
cgaln2lib_body(char *input_file_name, char *lib_file_name/*, char **seq_names, int num_sequences*/)
{
	FILE *lib_f = fopen(lib_file_name, "a");
	FILE *aln_f = fopen(input_file_name, "r");

	if (aln_f == NULL) {
		fprintf(stderr, "cgaln2lib_body gives following error:");
		perror(input_file_name);
		exit(EXIT_FAILURE);
	}

	strtok (input_file_name, "_");
	int seq1 = atoi(strtok (NULL, "_"))+1;
	int seq2 = atoi(strtok (NULL, "_"))+1;
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];
	int start1, start2, end1;
	fprintf(lib_f, "#%i %i\n", seq1, seq2);
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] != '#')
		{
			strtok ( line,":-,");
			start1 = atoi(strtok ( NULL,":-,"))+1;
			end1 = atoi(strtok ( NULL,":-,"))+1;
			strtok ( NULL,":-,");
			start2 = atoi(strtok ( NULL,":-,"))+1;
			fprintf(lib_f, "+BLOCK+ %i %i %i 100\n",  end1-start1+1, start1, start2);

		}
	}
	fclose(aln_f);
	fclose(lib_f);
}



int
findandgetvalue(char *query, int query_len, FILE *target)
{
	const int LINE_LENGTH = 200;
	char line[LINE_LENGTH];
	char *tmp;
	while (fgets(line, LINE_LENGTH, target) != NULL)
	{

		tmp = strstr(line, query);
		if (tmp != NULL)
		{
			return atoi(tmp + query_len +1);
		}
	}
	return -1;
}


void
tblastx2lib_body(char *input_f, char *lib_f)
{
	FILE *blast_F = fopen(input_f, "r");
	strtok (input_f, "_");
	int seq1 = atoi(strtok (NULL, "_"))+1;
	int seq2 = atoi(strtok (NULL, "_"))+1;


	FILE *lib_F = fopen(lib_f, "a");
	fprintf(lib_F, "#%i %i\n", seq1, seq2);

	int end = 0;
	int start_q, start_t, aln_len, frame_q, frame_t;

	int len_q = findandgetvalue("BlastOutput_query-len",21,  blast_F);
	int len_t = findandgetvalue("Hit_len",7,  blast_F);


	while (!end)
	{
		start_q = findandgetvalue("Hsp_query-from",14,  blast_F);
		if (start_q == -1)
		{
			end=1;
			continue;
		}

		start_t = findandgetvalue("Hsp_hit-from", 12, blast_F);

		frame_q= findandgetvalue("Hsp_query-frame", 15, blast_F);
		frame_t= findandgetvalue("Hsp_hit-frame", 13, blast_F);

		aln_len = findandgetvalue("Hsp_align-len", 13, blast_F);

		if (frame_q < 0)
			start_q = len_q - start_q -aln_len +1;
		if (frame_t < 0)
			start_t = len_t - start_t -aln_len +1;

		fprintf(lib_F, "++ %i %i %i 100\n", aln_len*3, start_q, start_t);
	}
	fclose(lib_F);
	fclose(blast_F);
}




void
lastz2lib_body(char *input_file_name, char *lib_file_name, char **seq_names, int num_sequences)
{
	FILE *lib_f = fopen(lib_file_name, "a");
	FILE *aln_f = fopen(input_file_name, "r");
	char *name;
	int seq_num1 = -1, seq_num2 = -1;
	char *tmp;
	const int LINE_LENGTH = 100;
	char line[LINE_LENGTH];
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] == 's')
		{
			fgets(line, LINE_LENGTH, aln_f);
			name = strtok ( line,"\" ");
// 			printf("%s\n", name);
			seq_num1 = get_seq_num2(seq_names, num_sequences, name);
			fgets(line, LINE_LENGTH, aln_f);
			name = strtok ( line,"\" ");
			seq_num2 = get_seq_num2(seq_names, num_sequences, name);
			break;
		}
	}
	int start1, start2, end1, length;

	fprintf(lib_f, "#%i %i\n", seq_num1, seq_num2);
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] != '#')
		{
			if (line[2] == 'l')
			{
				tmp=&line[3];
				start1 = atoi(strtok ( tmp," "));
				start2 = atoi(strtok ( NULL," "));
				end1 = atoi(strtok ( NULL," "));
				length = end1 - start1 +1;
				fprintf(lib_f, "++ %i %i %i 100\n",  length, start1, start2);
			}
		}
		else
		{
			break;
		}
	}

	fclose(aln_f);
	fclose(lib_f);
}



void
lastz2lib_body_2(char *input_file_name, FILE *lib_F, unsigned int off_set1, unsigned int off_set2)
{
	FILE *aln_f = fopen(input_file_name, "r");

	int seq_num1 = -1, seq_num2 = -1;
	char *tmp;
	const int LINE_LENGTH = 100;
	char line[LINE_LENGTH];

	strtok(input_file_name, "_");
	seq_num1 = atoi(strtok(NULL, "_"));
	seq_num2 = atoi(strtok(NULL, "_"));
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
		if (line[0] == 's')
			break;

	int start1, start2, end1, length;

	fprintf(lib_F, "#%i %i\n", seq_num1, seq_num2);
	while (fgets(line, LINE_LENGTH, aln_f) != NULL)
	{
		if (line[0] != '#')
		{
			if (line[2] == 'l')
			{
				tmp=&line[3];
				start1 = atoi(strtok ( tmp," "));
				start2 = atoi(strtok ( NULL," "));
				end1 = atoi(strtok ( NULL," "));
				length = end1 - start1 +1;
				fprintf(lib_F, "++ %i %i %i 100\n",  length, start1+off_set1, start2+off_set2);
			}
		}
		else
		{
			break;
		}
	}
	fclose(aln_f);
}







void
make_lib(char **seq_files, int num_seq_files, char **aln_files, int num_aln_files, char *out_file)
{

	make_lib_header(seq_files, num_seq_files, out_file);

	int i;
	int problems_occured = 0;
	for (i = 0; i < num_aln_files; ++i)
	{
		if (aln_files[i] != NULL)
		{
			if (aln_files[i][0] == '1')
			{
				cgaln2lib_body(aln_files[i], out_file/*, seq_files, num_seq_files*/);
			}
			else if (aln_files[i][0] == '2')
			{
				lastz2lib_body(aln_files[i], out_file, seq_files, num_seq_files);
			}
			else if (aln_files[i][0] == '3')
			{
				if (xmfa2lib_body(aln_files[i], out_file, /*seq_files,*/ num_seq_files))
				{
					problems_occured = 1;
				}
			}
			else if ((aln_files[i][0] == '4') || (aln_files[i][0] == '5'))
			{
				fasta2lib_body(aln_files[i], out_file, seq_files, num_seq_files);
			}
			else if (aln_files[i][0] == '6')
			{
				tblastx2lib_body(aln_files[i], out_file);
			}
			else if (aln_files[i][0] == 7)
			{
				maf2lib_body(aln_files[i], out_file, seq_files, num_seq_files);
			}
			else
			{
				printf("PROBLEM\n");
			}
		}
	}

	make_lib_end(out_file);
	if (problems_occured)
	{
		fprintf(stderr, "Problems occured while producing library!\n");
	}
}



Alignment*
_read_xmfa_block(FILE *xmfa_F)
{
	const int LINE_LENGTH = 301;
	char line[LINE_LENGTH];

	Alignment *block = malloc(sizeof(Alignment));
	unsigned int aln_space = 10;
	block->seqs = malloc(aln_space * sizeof(Sequence*));
	unsigned int num_seqs = 0;
	Sequence *seq = NULL;
	char *sequence = NULL;
	unsigned int i;
	const unsigned int STEP_SIZE = 500;
	unsigned int reserved = STEP_SIZE,  used = 0;
	char *tmp;
	unsigned int gap_less_length = 0;
	while (fgets(line, LINE_LENGTH, xmfa_F) != NULL)
	{
		if (line[0] == '>')
		{
			if (num_seqs > 0)
			{
				sequence[used]='\0';
				seq->seq = sequence;
				seq->length = used;
				seq->gap_less_length = gap_less_length;
			}
			if (num_seqs >= aln_space)
			{
				aln_space+=10;
				block->seqs = realloc(block->seqs, aln_space * sizeof(Sequence*));
			}
			block->seqs[num_seqs] = malloc(sizeof(Sequence));
			seq = block->seqs[num_seqs];

			// XMFA sequence header
			// > seq_num:start1-end1 ± comments (sequence name, etc.)
			seq->number = atoi(strtok (&line[1], ":"));
			seq->name = NULL;
			seq->chr = NULL;
			seq->start=atoi(strtok (NULL, "-"))-1;
			seq->end=atoi(strtok (NULL, " "))-1;

			if ((strtok (NULL, " ")[0]) == '+')
				seq->strand = 1;
			else
				seq->strand = -1;

			tmp=strtok (NULL, "\n");
			if (tmp != NULL)
				seq->comment = malloc((strlen(tmp)+1) * sizeof(char));
			else
				seq->comment = NULL;
			strcpy(seq->comment, tmp);
			++num_seqs;
			reserved = STEP_SIZE;
			sequence = malloc(reserved * sizeof(char));
			used=0;
			gap_less_length = 0;
		}
		else if (line[0] == '=')
		{
			sequence[used]='\0';
			seq->seq = sequence;
			seq->length = used;
			seq->gap_less_length = gap_less_length;
			block->num_seqs = num_seqs;
			block->seqs = realloc(block->seqs, num_seqs * sizeof(Sequence));
			return block;
		}
		else
		{
			if (used + LINE_LENGTH > reserved)
			{
				reserved += STEP_SIZE;
				sequence = realloc(sequence, reserved* sizeof(char));
			}
			i = 0;
			while ((line[i] != '\n') && (line[i] != '\0'))
			{
				if (line[i] != '-')
					++gap_less_length;
				sequence[used++] = line[i];
				++i;
			}
		}
	}
	free(block->seqs);
	free(block);
	return NULL;
}

