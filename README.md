Robusta
=======


Robusta is a meta multiple genome alignment package. 

It does not produce any alignment by itself but can combine alignments of different genome alignment programs 
into a single final alignment. It is able to call serveral different pairwise and multiple genome alignment programs. 
The results of these programs are parsed and turned into a T-Coffee library. 

Using T-Coffee the alignments are then turned into the final alignment. Currently two pairwise genome aligner 
and three multiple genome alignment programs are included into this framework. A typical usage would be to cut 
the genomes to align into a set of collinear homologous blocks using Mercator. 

These blocks are than aligned and turned into libraries using Robusta. As a final step T-Coffee is called to 
produce an alignment for each block.


Read more: http://www.tcoffee.org/Projects/robusta/
