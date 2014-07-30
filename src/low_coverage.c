#include "low_coverage.h"



void
decode(Sequence **seqs, char *tmp_out_f, char *out_f)
{
	FILE *tmp_out_F = fopen(tmp_out_f, "r");
	FILE *out_F = fopen(out_f, "w");
	char strand;
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];
	unsigned int seq_num;
	while (fgets(line, LINE_LENGTH, tmp_out_F) != NULL)
	{
		if (line[0] == '>')
		{
			seq_num = atoi(&line[1]);
			if (seqs[seq_num]->strand==1)
				strand = '+';
			else
				strand = '-';
			fprintf(out_F, ">%i:%i-%i %c %s\n", seq_num, seqs[seq_num]->start+1, seqs[seq_num]->end+1, strand, seqs[seq_num]->comment);
		}
		else
			fprintf(out_F, "%s", line);
	}
	fprintf(out_F, "=\n");
	fclose(out_F);
	fclose(tmp_out_F);
}



int
prepare_low_cov(char *aln_f, char *out_f, int do_missing_blasts, char *ref_tree_f, int concatenate_species, int n_cores)
{
	char *lib_f = "library";
	char *tree_f = "guide_tree";
	char *tmp_out_f = "t_coffee_out_coded";



	unsigned int num_seqs;
	if (!concatenate_species)
	{
		Sequence **seqs = NULL;
		num_seqs = produce_low_cow_library(aln_f, lib_f, &seqs, do_missing_blasts, ref_tree_f);
		if (ref_tree_f == NULL)
			produce_tree(tree_f, num_seqs);
		char command[500];
		sprintf(command, "t_coffee -n_core %i -quiet -lib %s -usetree %s -outfile %s -output fasta_aln >/dev/null 2>/dev/null", n_cores, lib_f, tree_f, tmp_out_f);
		int status = system(command);
		decode(seqs, tmp_out_f, out_f);

		unsigned int i;
		for ( i = 1; i <= num_seqs; ++i)
		{
			free_Sequence(seqs[i]);
			free(seqs[i]);
		}
		free(seqs);
		return status;
	}
	else
		return concatenate_lib(aln_f, lib_f, out_f, do_missing_blasts, ref_tree_f, n_cores);
}



void
produce_tree(char *tree_f, unsigned int num_seqs)
{
	FILE *tree_F = fopen(tree_f, "w");
	unsigned int i;
	for (i = 0; i < num_seqs-1; ++i)
	{
		fprintf(tree_F, "(");
	}
	fprintf(tree_F, "1, 2)");
	for (i = 3; i <= num_seqs; ++i)
	{
		fprintf(tree_F, ",%i)", i);
	}
	fprintf(tree_F, ";\n");
	fclose(tree_F);
}



int compare (const void * a, const void * b)
{
	int x = ((Pair*)a)->v1 - ((Pair*)b)->v1;
	if (x == 0)
		return (((Pair*)a)->v2 - ((Pair*)b)->v2);
	else
		return x;
}




/**
 * \brief Produces a guide tree incorporating a guide tree for the reference alignment.
 *
 * \param seqs The sequences
 * \param num_ref_seqs The number of reference sequences
 * \param num_seqs The allover number of sequences
 * \param ref_tree_f File to save the guide tree
 *
 */
void
make_reference_guide_tree(Sequence **seqs, unsigned int num_ref_seqs, unsigned int num_seqs, char *ref_tree_f)
{
	FILE *ref_tree_F = fopen(ref_tree_f, "r");
	FILE *new_tree_F = fopen("guide_tree", "w");
	char tree_string[1000];
	fgets(tree_string, 1000, ref_tree_F);

	fclose(ref_tree_F);
	unsigned int tree_len = strlen(tree_string)-1;
	unsigned int i = 0;
	unsigned int j, k;
	char c;

	for (i= num_ref_seqs; i < num_seqs; ++i)
	{
		fprintf(new_tree_F, "(");
	}

	i =0;
	while (i < tree_len)
	{
		if ((tree_string[i] == '(') || (tree_string[i] == ')') || (tree_string[i] == ' ') || (tree_string[i] == ',') || (tree_string[i] == ';'))
		{
			fprintf(new_tree_F, "%c", tree_string[i]);
			++i;
		}
		else
		{
			j = i;
			while ((tree_string[i] != '(') && (tree_string[i] != ')') && (tree_string[i] != ' ') && (tree_string[i] != ',') && (tree_string[i] != ';'))
			{
				++i;
			}
			c = tree_string[i];
			tree_string[i] = '\0';
			for (k = 1; k <= num_ref_seqs; ++k)
			{
				if (!strncmp(seqs[k]->comment, &tree_string[j], strlen(&tree_string[j])))
				{
					fprintf(new_tree_F, "%i", seqs[k]->number);
					break;
				}

			}
			tree_string[i] = c;
		}
	}


	for (i = num_ref_seqs+1; i <= num_seqs; ++i)
	{
		fprintf(new_tree_F, ",%i)", i);
	}
	fprintf(new_tree_F, ";\n");
	fclose(new_tree_F);
}



int
produce_low_cow_library(char *aln_f, char *lib_f, Sequence ***result_seq, int do_missing_blasts, char *ref_tree_f)
{
	Alignment *ref_aln, *pw_aln;
	FILE *lib_F = fopen(lib_f, "w");
	if (lib_F == NULL)
	{
		fprintf(stderr, "ERROR: Could not open file: %s\n", lib_f);
		exit(1);
	}
	FILE *aln_F = fopen(aln_f, "r");
	if (aln_F == NULL)
	{
		fprintf(stderr, "ERROR: Could not open file: %s\n", aln_f);
		exit(1);
	}
	fprintf(lib_F, "! TC_LIB_FORMAT_01\n");
	unsigned int num_pos = ftell(lib_F);
	fprintf(lib_F, "          \n"); //This keeps empty space to be replaced later with the number of sequneces

	//turning reference alignment into library
	ref_aln = _read_xmfa_block(aln_F);



	char *seq;
	FILE *single_seq_F;
	char single_seq_f[8];
	unsigned int num_ref_seqs = ref_aln->num_seqs;
	unsigned int i, j;
	for (i=0; i < num_ref_seqs; ++i)
	{
		seq=get_gap_free_seq(ref_aln, i);
		if (do_missing_blasts)
		{
			sprintf(single_seq_f, "%i", ref_aln->seqs[i]->number);
			single_seq_F = fopen(single_seq_f, "w");
			fprintf(single_seq_F, ">%i\n%s\n", ref_aln->seqs[i]->number, seq);
			fclose(single_seq_F);
		}
		fprintf(lib_F, "%i %lu %s\n", ref_aln->seqs[i]->number, strlen(seq), seq);
		free(seq);
	}


	unsigned int lib_reserved = (num_ref_seqs*(num_ref_seqs-1)/2) + 100;
	unsigned int in_use = 0;
	Library_pair **pairs= malloc(lib_reserved*sizeof(Library_pair*));
	unsigned int seq_reserved = num_ref_seqs + 50;
	unsigned int seqs_use = 1;
	Sequence **seqs = malloc(seq_reserved * sizeof(Sequence*));

	// turning reference alignment into library
	for (i = 0; i < num_ref_seqs; ++i)
	{
		for (j = i+1; j < num_ref_seqs; ++j)
		{
			pairs[in_use++] = aln_2_pairs(ref_aln, i , j, ref_aln->seqs[i]->start, ref_aln->seqs[j]->start);
		}
		free(ref_aln->seqs[i]->seq);
		seqs[seqs_use++] = ref_aln->seqs[i];
		ref_aln->seqs[i]->seq = NULL;
	}



	free(ref_aln->seqs);
	free(ref_aln);
	unsigned int ref_seq_num;
	unsigned int ref_index, low_index;


	unsigned int storage_reserve = 100;
	Pair *pair_storage = NULL;
	if (do_missing_blasts)
		pair_storage = malloc(storage_reserve * sizeof(Pair));
	unsigned int storage_use =0;





	// turn pair low coverage libraries into library
	while ((pw_aln = _read_xmfa_block(aln_F)) != NULL)
	{

		//allocate more memory if needed
		if (storage_use >= storage_reserve)
		{
			storage_reserve += 100;
			pair_storage = realloc(pair_storage, storage_reserve*sizeof(Pair));
		}
		if (in_use >= lib_reserved)
		{
			lib_reserved += 100;
			pairs = realloc(pairs, lib_reserved*sizeof(Library_pair*));
		}
		if (seqs_use >= seq_reserved)
		{
			seq_reserved += 50;
			seqs = realloc(seqs, seq_reserved*sizeof(Sequence*));
		}


		ref_seq_num = pw_aln->seqs[0]->number;
		if (ref_seq_num < num_ref_seqs)
		{
			ref_index = 0;
			low_index = 1;
		}
		else
		{
			ref_index = 1;
			low_index = 0;
		}


		pairs[in_use++] = aln_2_pairs(pw_aln, ref_index , low_index, seqs[ref_seq_num]->start, pw_aln->seqs[low_index]->start);
		seq=get_gap_free_seq(pw_aln, low_index);
		fprintf(lib_F, "%i %lu %s\n", pw_aln->seqs[low_index]->number, strlen(seq), seq);
		seqs[seqs_use++] = pw_aln->seqs[low_index];
		if (do_missing_blasts)
		{
			sprintf(single_seq_f, "%i", pw_aln->seqs[low_index]->number);
			single_seq_F = fopen(single_seq_f, "w");
			fprintf(single_seq_F, ">%i\n%s\n", pw_aln->seqs[low_index]->number, seq);
			fclose(single_seq_F);
			pair_storage[storage_use].v1=pw_aln->seqs[ref_index]->number;
			pair_storage[storage_use++].v2=pw_aln->seqs[low_index]->number;
		}

		//free alignment memory but not "Sequence meta_data", this is stored in array.
		free(seq);
		free_Sequence(pw_aln->seqs[ref_index]);
		free(pw_aln->seqs[ref_index]);
		free(pw_aln->seqs[low_index]->seq);
		pw_aln->seqs[low_index]->seq = NULL;
		free(pw_aln->seqs);
		free(pw_aln);
	}
	fclose(aln_F);
	--seqs_use;

	// write read_in_pairs
	Library_pair *pair = NULL;
	unsigned int num_pairs;
	for (i = 0; i < in_use; ++i)
	{
		pair = pairs[i];
		fprintf(lib_F, "#%i %i\n", pair->seq1, pair->seq2);
		num_pairs = pair->num_entries;
		j = 0;
		while (j < num_pairs)
		{
			fprintf(lib_F, "%i %i 100\n", pair->pairs[j], pair->pairs[j+1]);
			j+=2;
		}
		free_Library(pair);
	}
	free(pairs);

	// produce missing lastz_libs
	// unsigned int num_low_cov_seqs = seqs_use - num_ref_seqs;
	char file1[10], file2[10], lastz_out_f[30];

	if (do_missing_blasts)
	{
		unsigned int current_pos =0;
		unsigned int low_cov_pos;
		qsort (pair_storage, storage_use, sizeof(Pair), compare);

		for (i = 1; i <= num_ref_seqs; ++i)
		{
			low_cov_pos= num_ref_seqs+1;
			for (j = low_cov_pos; j <= seqs_use; ++j)
			{
				if ((current_pos < storage_use) && (pair_storage[current_pos].v1 == i) && (pair_storage[current_pos].v2 == j))
				{
					++current_pos;
					continue;
				}
				sprintf(file1, "%i", i);
				sprintf(file2, "%i", j);
				sprintf(lastz_out_f, "2_%i_%i", i, j);
				lastz_pw(file1, file2, lastz_out_f, 3000);
			}
		}

		//wait for all lastz finished
		wait_for_all();


		//read libraries
		current_pos = 0;
		for (i = 1; i <= num_ref_seqs; ++i)
		{
			low_cov_pos= num_ref_seqs+1;
			for (j = low_cov_pos; j <= seqs_use; ++j)
			{
				if ((current_pos < storage_use) && (pair_storage[current_pos].v1 == i) && (pair_storage[current_pos].v2 == j))
				{
					++current_pos;
					continue;
				}
				sprintf(lastz_out_f, "2_%i_%i", i, j);
				lastz2lib_body_2(lastz_out_f, lib_F, 0, 0);
			}
		}

		free(pair_storage);
	}

	if (ref_tree_f != NULL)
		make_reference_guide_tree(seqs, num_ref_seqs, seqs_use, ref_tree_f);

	fprintf(lib_F, "! SEQ_1_TO_N\n");
	fseek(lib_F, num_pos, SEEK_SET);
	fprintf(lib_F, "%i", seqs_use);
	fclose(lib_F);
	result_seq[0] = seqs;
	return seqs_use;
}