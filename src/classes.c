#include "classes.h"


//Struct functions
void
free_Sequence(Sequence *seq)
{
	if (seq->name != NULL)
		free(seq->name);
	if (seq->chr !=NULL)
		free(seq->chr);
	if (seq->seq !=NULL)
		free(seq->seq);
	if (seq->comment !=NULL)
		free(seq->comment);
}


void
print_seq(Sequence *seq)
{
	printf("%s %s %i %i %i\n %s\n", seq->name, seq->chr, seq->start, seq->end, seq->strand, seq->seq);
}


// Alignment functions
void
free_Alignment(Alignment *aln)
{
	unsigned int i;
	unsigned int num_seqs = aln->num_seqs;
	for (i = 0; i < num_seqs; ++i)
	{
		free_Sequence(aln->seqs[i]);
		free(aln->seqs[i]);
	}
	free(aln->seqs);
}


char*
get_gap_free_seq(Alignment *aln, unsigned int in)
{
	unsigned int i, pos = 0;
	char *seq = aln->seqs[in]->seq;
	unsigned int l = aln->seqs[in]->length;
	char *gap_free_seq = malloc((l+1) * sizeof(char));
	for (i =0; i < l; ++i)
	{
		if (seq[i] != '-')
		{
			gap_free_seq[pos++] = seq[i];
		}
	}
	gap_free_seq[pos]='\0';
	gap_free_seq = realloc(gap_free_seq, (++pos)*sizeof(char));
	return gap_free_seq;
}

// Library_pair functions
void
free_Library_content(Library_pair *lib)
{
	if (lib->pairs != NULL)
		free(lib->pairs);
}


void
free_Library(Library_pair *lib)
{
	if (lib->pairs != NULL)
		free(lib->pairs);
	free(lib);
}


void
aln_2_lib_body(const Alignment *aln, int index1, int index2, int off_set1, int off_set2, FILE *lib_F)
{
	char *seq1=aln->seqs[index1]->seq;
	char *seq2=aln->seqs[index2]->seq;

	unsigned int len= aln->seqs[index1]->length;
	unsigned int s1_pos = aln->seqs[index1]->start-off_set1+1;
	unsigned int s2_pos = aln->seqs[index2]->start-off_set2+1;
	unsigned int i;
	fprintf(lib_F, "#%i %i\n", aln->seqs[index1]->number, aln->seqs[index2]->number);
	for (i = 0; i < len; ++i)
	{
		if (( seq1[i] != '-' ) && ( seq2[i] != '-' ))
			fprintf(lib_F, "%i %i 100\n", s1_pos, s2_pos);
		if ( seq1[i] != '-' )
			++s1_pos;
		if ( seq2[i] != '-' )
			++s2_pos;
	}


}




Library_pair*
aln_2_pairs(const Alignment *aln, int index1, int index2, int off_set1, int off_set2)
{
	char *seq1=aln->seqs[index1]->seq;
	char *seq2=aln->seqs[index2]->seq;

	Library_pair *lib = malloc(sizeof(Library_pair));
	lib->seq1 = aln->seqs[index1]->number;
	lib->seq2 = aln->seqs[index2]->number;
	unsigned int len= aln->seqs[index1]->length;
	unsigned int *pairs = malloc(len*2 * sizeof(unsigned int));
	unsigned int s1_pos = aln->seqs[index1]->start-off_set1+1;
	unsigned int s2_pos = aln->seqs[index2]->start-off_set2+1;
	unsigned int i, lib_pos = 0;
	for (i = 0; i < len; ++i)
	{
		if (( seq1[i] != '-' ) && ( seq2[i] != '-' ))
		{
			pairs[lib_pos++] = s1_pos;
			pairs[lib_pos++] = s2_pos;
		}
		if ( seq1[i] != '-' )
			++s1_pos;
		if ( seq2[i] != '-' )
			++s2_pos;
	}
	lib->pairs = realloc(pairs, lib_pos * sizeof(unsigned int)); //shrinking memory
	lib->num_entries = lib_pos;

	return lib;
}






void
write_pairs(const Library_pair *pair, FILE *lib_F)
{
	fprintf(lib_F, "#%i %i\n", pair->seq1, pair->seq2);
	unsigned int i = 0;
	unsigned int len = pair->num_entries;
	unsigned int *pairs = pair->pairs;
	while (i < len)
	{
		fprintf(lib_F, "%i %i\n", pairs[i], pairs[i+1]);
		i+=2;
	}
}
