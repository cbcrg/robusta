#include "low_coverage_concatenate.h"



Irreg_2DSorter *
make_irregular_2D(unsigned int dim1)
{
	unsigned int i;
	Irreg_2DSorter *irr = malloc(sizeof(Irreg_2DSorter));
	irr->dim1 = dim1;
	irr->dim2 = malloc(2*dim1 * sizeof(int));
	irr->data = malloc(dim1 *sizeof(Sorter*));
	irr->step_size = 20;
	for (i = 0; i < dim1; ++i)
	{
		irr->dim2[i*2] = 0;
		irr->dim2[i*2+1] = irr->step_size;
		irr->data[i] = malloc(irr->step_size * sizeof(Sorter));
	}
	irr->max = 0;
	return irr;

}


void
free_Irreg2DSorter(Irreg_2DSorter *irr)
{
	unsigned int i = 0;
	for (i = 0; i < irr->dim1; ++i)
		free(irr->data[i]);
	free(irr->data);
	free(irr->dim2);
	free(irr);
}



void
push_sorter(Irreg_2DSorter *manager, unsigned int pos, unsigned int ending, unsigned int id)
{
// 	printf("id %i\n", id);
	unsigned int i;
	if (id >= manager->max)
		manager->max = id+1;
	if (id >= manager->dim1)
	{
		unsigned int new_dim = manager->dim1 + 20;
		manager->data = realloc(manager->data,new_dim * sizeof(Irreg_2DSorter*));
		manager->dim2 = realloc(manager->dim2, 2*new_dim * sizeof(int));
		for (i = manager->dim1; i < new_dim; ++i)
		{
			manager->dim2[i*2] = 0;
			manager->dim2[i*2+1] = manager->step_size;
			manager->data[i] = malloc(manager->step_size * sizeof(Sorter));
		}
		manager->dim1 = new_dim;
	}


	if (manager->dim2[2*id] == manager->dim2[2*id+1])
	{
		manager->dim2[2*id+1]+= manager->step_size;
		manager->data[id] = realloc(manager->data[id], manager->dim2[2*id+1] * sizeof(Sorter));
	}

	Sorter *tmp = &(manager->data[id][manager->dim2[2*id]]);
	tmp->pos=pos;
	tmp->aln_start = ending;
	++manager->dim2[2*id];
}





int **
seq_2_aln(Alignment *aln)
{
	unsigned int num_seqs = aln->num_seqs;
	int **hash = malloc(num_seqs * sizeof(int*));
	int *hash_tmp;
	unsigned int i, j, k, len;
	char *seq;
	for (i = 0; i < num_seqs; ++i)
	{
		len = aln->seqs[i]->length;
		hash[i] = malloc(aln->seqs[i]->gap_less_length * sizeof(int));
		hash_tmp = hash[i];
		seq= aln->seqs[i]->seq;
		k = 0;
		for (j = 0; j < len; ++j)
		{
			if (seq[j] != '-')
				hash_tmp[k++] = j;
		}
	}
	return hash;
}



int sequence_sort_compare (const void * a, const void * b)
{
	return ((Sorter*)a)->aln_start - ((Sorter*)b)->aln_start;
}


void
sort_sequences(Irreg_2DSorter *manager)
{
	unsigned int i;
	for (i = 0; i < manager->max; ++i)
	{
		qsort(manager->data[i], manager->dim2[i*2], sizeof(Sorter), sequence_sort_compare);
	}
}




int
concatenate_seqs(Sorter *data, unsigned int size, Alignment **alns, char **seq_p, unsigned int *reserved)
{
	char *seq = *seq_p;
	unsigned int pos = 0;
	unsigned int i = 0;
	unsigned int len;
	unsigned int max_len = *reserved;

	char *gap_free_seq = NULL;
	Sequence *current_seq = NULL;
	Sequence *next_seq = alns[data[0].pos]->seqs[0];
	for (i = 0; i < size-1; ++i)
	{
		current_seq = next_seq;
		next_seq = alns[data[i+1].pos]->seqs[0];
		if (current_seq->end  <  next_seq->start)
		{
			data[i].off_set = pos;
			gap_free_seq = get_gap_free_seq(alns[data[i].pos], 1);
			len = alns[data[i].pos]->seqs[1]->gap_less_length; //strlen(gap_free_seq);
			if (pos + len >= max_len)
			{
				max_len += len*2;
				seq = realloc(seq, max_len * sizeof(char));
			}
			strcpy(&seq[pos], gap_free_seq);
			pos += len;
			free(gap_free_seq);
		}
		else
		{
			printf("%i %i %i %i\n", current_seq->start, current_seq->end, next_seq->start, next_seq->end);
			fprintf(stderr, "WARNING - Concatenation problem\n");
		}

	}
	gap_free_seq = get_gap_free_seq(alns[data[i].pos], 1);
	len = alns[data[i].pos]->seqs[1]->gap_less_length;
	data[i].off_set = pos;
	if (pos + len >= max_len)
	{
		max_len += len;
		seq = realloc(seq, max_len * sizeof(char));
	}


	strcpy(&seq[pos], gap_free_seq);
	pos += strlen(gap_free_seq);
	free(gap_free_seq);


	*reserved = max_len;
	seq_p[0] = seq;
	seq[pos] = '\0';
	return pos;
}



void
split_sequence(Sorter *sort, unsigned int num_pieces, char *aln_seq, unsigned int aln_length, Alignment **pw_alns, FILE *out_F)
{
	unsigned int tmp_len, i, j, aln_pos = 0;
	unsigned int seq_pos;
	Sequence *tmp_seq;
	char strand;
	for (i = 0; i < num_pieces; ++i)
	{
		seq_pos = 0;
		tmp_seq = pw_alns[sort[i].pos]->seqs[1];
		if (tmp_seq->strand==1)
			strand = '+';
		else
			strand = '-';
		tmp_len = tmp_seq->end - tmp_seq->start +1;

		fprintf(out_F, ">%i:%i-%i %c %s\n", tmp_seq->number, tmp_seq->start+1, tmp_seq->end+1, strand, tmp_seq->comment);

		//first gaps
		for (j = 0; j < aln_pos; ++j)
			fprintf(out_F, "-");

		//sequence (with gaps)
		while (seq_pos < tmp_len)
		{
			if (aln_seq[aln_pos] != '-')
				++seq_pos;
			fprintf(out_F, "%c", aln_seq[aln_pos++]);
		}
		j=aln_pos;
		//end gaps
		for (; j < aln_length; ++j)
			fprintf(out_F, "-");
		fprintf(out_F, "\n");
	}

}



void
decode_and_split(Alignment *ref_aln, Irreg_2DSorter *manager, Alignment **pw_alns, char *tmp_out_f, char *out_f)
{
	FILE *tmp_out_F = fopen(tmp_out_f, "r");
	FILE *out_F = fopen(out_f, "w");
	char strand;
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];
	unsigned int num_ref_seqs = ref_aln->num_seqs;
	unsigned int num =0;
	Sequence *tmp_seq;

	//Decode reference alignment
	while (fgets(line, LINE_LENGTH, tmp_out_F) != NULL)
	{
		if (line[0] == '>')
		{
			if (num == num_ref_seqs)
				break;
			tmp_seq = ref_aln->seqs[num];
			if (tmp_seq->strand==1)
				strand = '+';
			else
				strand = '-';
			fprintf(out_F, ">%i:%i-%i %c %s\n", ++num, tmp_seq->start+1, tmp_seq->end+1, strand, tmp_seq->comment);
		}
		else
			fprintf(out_F, "%s", line);
	}

	unsigned int reserved = 1000;
	char *aln_seq = malloc(reserved * sizeof(char));
	unsigned int end, used = 0;
	num = 0;


	//read low coverage sequence and align it to the output
	while (fgets(line, LINE_LENGTH, tmp_out_F) != NULL)
	{
		if (line[0] == '>')
		{
			split_sequence(manager->data[num], manager->dim2[num*2], aln_seq, used, pw_alns, out_F);

			++num;
			used=0;
		}
		else
		{
			if (used + LINE_LENGTH > reserved)
			{
				reserved += 1000;
				aln_seq = realloc(aln_seq, reserved * sizeof(char));
			}
			end = strlen(line)-1;
			if (line[end] == '\n')
				line[end] = '\0';
			else
				++end;
			strcpy(&aln_seq[used], line);
			used += end;
		}
	}
	split_sequence(manager->data[num], manager->dim2[num*2], aln_seq, used, pw_alns, out_F);

	free(aln_seq);
	fprintf(out_F, "=\n");
	fclose(out_F);
	fclose(tmp_out_F);
}




void
make_reference_guide_tree2(Alignment *ref_aln, unsigned int num_seqs, char *ref_tree_f, char *tree_out)
{
	FILE *ref_tree_F = fopen(ref_tree_f, "r");
	FILE *new_tree_F = fopen(tree_out, "w");
	char tree_string[1000];
	fgets(tree_string, 1000, ref_tree_F);

	fclose(ref_tree_F);
	unsigned int tree_len = strlen(tree_string)-1;
	unsigned int i = 0;
	unsigned int j, k;
	char c;
	unsigned int num_ref_seqs = ref_aln->num_seqs;
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
			for (k = 0; k < num_ref_seqs; ++k)
			{
				if (!strncmp(ref_aln->seqs[k]->comment, &tree_string[j], strlen(&tree_string[j])))
				{
					fprintf(new_tree_F, "%i", ref_aln->seqs[k]->number);
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
concatenate_lib(char *aln_f, char *lib_f, char *out_f, int do_missing_pairs, char *ref_tree_f, int n_cores)
{

	//open neede files
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


	//get reference alignment and write sequences into library
	Alignment *ref_aln = _read_xmfa_block(aln_F);
	unsigned int num_ref_seqs = ref_aln->num_seqs;
	unsigned int i,j;
	char *seq;
	char single_seq_f[10];
	FILE *single_seq_F;
	for (i=0; i < num_ref_seqs; ++i)
	{
		seq=get_gap_free_seq(ref_aln, i);
		if (do_missing_pairs)
		{
			sprintf(single_seq_f, "%i", ref_aln->seqs[i]->number);
			single_seq_F = fopen(single_seq_f, "w");
			fprintf(single_seq_F, ">%i\n%s\n", ref_aln->seqs[i]->number, seq);
			fclose(single_seq_F);
		}
		fprintf(lib_F, "%i %lu %s\n", ref_aln->seqs[i]->number, strlen(seq), seq);
		free(seq);
	}


	int **seq_2_aln_h = seq_2_aln(ref_aln);
	Alignment *pw_aln;
	Sequence *ref_seq, *low_cov_seq;
	unsigned int max_aln_storage = 100;
	unsigned int aln_storage_usage = 0;
	Alignment **pw_alns  = malloc(max_aln_storage * sizeof(Alignment*));
	unsigned int ref_seq_number;
	Sequence *tmp;


	//Read pairwise alignment
	Irreg_2DSorter *manager =  make_irregular_2D(50);
	while ((pw_aln = _read_xmfa_block(aln_F)) != NULL)
	{
		if (aln_storage_usage == max_aln_storage)
		{
			max_aln_storage += 100;
			pw_alns = realloc(pw_alns, max_aln_storage * sizeof(Alignment*));
		}
		pw_alns[aln_storage_usage] = pw_aln;

		//make ref_seq first sequence
		if (pw_aln->seqs[0]->number > num_ref_seqs)
		{
			tmp = pw_aln->seqs[0];
			pw_aln->seqs[0] = pw_aln->seqs[1];
			pw_aln->seqs[1] = tmp;
		}

		ref_seq = pw_aln->seqs[0];
		low_cov_seq = pw_aln->seqs[1];

		ref_seq_number = ref_seq->number-1;

		push_sorter(manager, aln_storage_usage, seq_2_aln_h[ref_seq_number][ref_seq->start-ref_aln->seqs[ref_seq_number]->start], low_cov_seq->number-num_ref_seqs-1);

		++aln_storage_usage;
	}
	fclose(aln_F);

	for (i = 0; i < num_ref_seqs; ++i)
		free(seq_2_aln_h[i]);
	free(seq_2_aln_h);

	sort_sequences(manager);

	unsigned int concat_length;
	unsigned int reserved = 1000;
	char **concat_seq = malloc(sizeof(char*));
	concat_seq[0] = malloc(reserved*sizeof(char*));

	for (i = 0; i < manager->max; ++i)
	{
		concat_length = concatenate_seqs(manager->data[i], manager->dim2[i*2], pw_alns, concat_seq, &reserved);
		fprintf(lib_F, "%i %i %s\n", pw_alns[manager->data[i]->pos]->seqs[1]->number, concat_length, concat_seq[0]);
	}
	free(concat_seq[0]);
	free(concat_seq);
	//write reference matches into library
	for (i =0; i < num_ref_seqs; ++i)
	{
		for (j = i+1; j < num_ref_seqs; ++j)
		{
			aln_2_lib_body(ref_aln, i , j, ref_aln->seqs[i]->start, ref_aln->seqs[j]->start, lib_F);
		}
	}


	//write low cov matches into library
	unsigned int dim, low_cov_id;
	Alignment *tmp_aln;
	Sorter *tmp_sorter;
	unsigned int ref_num, k;
	char *gap_free_seq;
	char file1[10], file2[10], lastz_out_f[20];
	for (i = 0; i < manager->max; ++i)
	{
		low_cov_id = pw_alns[manager->data[i]->pos]->seqs[1]->number;
		ref_num = pw_alns[manager->data[i]->pos]->seqs[0]->number;

		dim = manager->dim2[i*2];
		for (j = 0; j < dim; ++j)
		{
			tmp_sorter=&manager->data[i][j];
			tmp_aln = pw_alns[tmp_sorter->pos];
			aln_2_lib_body(tmp_aln, 0 , 1, ref_aln->seqs[tmp_aln->seqs[0]->number-1]->start, (int)tmp_aln->seqs[1]->start-(int)tmp_sorter->off_set, lib_F);


			//Produce missing lastz libraries if wanted
			if (do_missing_pairs)
			{
// 				low_cov_id = tmp_aln->seqs[1]->number;
				gap_free_seq = get_gap_free_seq(tmp_aln, 1);
				sprintf(file1, "%i", low_cov_id );
				FILE *lastz_F = fopen(file1, "w");
				fprintf(lastz_F, ">%i\n%s\n",tmp_aln->seqs[1]->number, gap_free_seq);
				fclose(lastz_F);
				free(gap_free_seq);
				for (k = 1; k <= num_ref_seqs; ++k)
				{
					if (k == ref_num)
						continue;
					sprintf(file2, "%i", k );
					sprintf(lastz_out_f, "2_%i_%i", low_cov_id, k);
					lastz_pw(file1, file2, lastz_out_f, 3000);
				}

				wait_for_all();

				for (k = 1; k <= num_ref_seqs; ++k)
				{
					if (k == ref_num)
						continue;
					sprintf(lastz_out_f, "2_%i_%i", low_cov_id, k);
					lastz2lib_body_2(lastz_out_f, lib_F, tmp_sorter->off_set, 0);
				}

			}


		}
	}



	fprintf(lib_F, "! SEQ_1_TO_N\n");
	fseek(lib_F, num_pos, SEEK_SET);
	unsigned int all_seq_num=manager->max+num_ref_seqs;
	fprintf(lib_F, "%i", manager->max+num_ref_seqs);
	fclose(lib_F);




	//Call T-Coffee
	char *tree_f = "guided_tree";

	if (ref_tree_f == NULL)
		produce_tree(tree_f, all_seq_num);
	else
		make_reference_guide_tree2(ref_aln, all_seq_num, ref_tree_f, tree_f);

	char *tmp_out_f = "t_coffee_out_coded";
	char command[500];
	sprintf(command, "t_coffee -quiet -n_core %i -lib %s -usetree %s -outfile %s -output fasta_aln >/dev/null 2>/dev/null", n_cores, lib_f, tree_f, tmp_out_f);
	int status = system(command);

	//split and decode
// 	decode_and_split(Alignment *aln_ref, Irreg_2DSorter *manager, Alignment *pw_alns, char *tmp_out_f, char *out_f)
	decode_and_split(ref_aln, manager, pw_alns, tmp_out_f, out_f);

	free_Alignment(ref_aln);
	free(ref_aln);
	for (i = 0; i < aln_storage_usage; ++i)
	{
		free_Alignment(pw_alns[i]);
		free(pw_alns[i]);
	}
	free(pw_alns);

	free_Irreg2DSorter(manager);
	return status;
}


