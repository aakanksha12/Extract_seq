# Extract_seq

	This pipeline is to extract loci from given taxa using nhmmer and blast

	This program requires folder named as "genomes","queries" and "hmm_prof_align" in current working directory 
	It requires 3 command line arguments passed while using it.

	 perl Extract_seq.pl <current working directory path> <text_file_database_list> <text_file_queries_list>
		
	<current working directory path> : It is the path of current working directory which contains genomes and queries folder
	<text_file_database_list> : It is text file that contains list of genomes in genomes folder without ".fasta". Each line should 		contain one genome name. It should be located in current working directory
	<text_file_queries_list>  : It is a text file that contains list of queries in the query folder without ".fasta". Each line should 		contain one query name.It should be located in current working directory
	genomes folder : It contains all the genome files after processing them using blast makedb  
	queries folder : It contains all the query files in .fasta format
	hmm_prof_align folder : It will store the alignment files in ".fasta" format(used by nhmmer) for all the loci
	Following directories are created in CWD as part of the program:
	gallus_blast: stores blast result with the given loci and gallus genome
	other_blast: stores blast result of given loci and all genomes
	extract_seq: stores extracted sequence and scaffold of best hit
	blast_again: stores blast result of extracted sequence and gallus genome
