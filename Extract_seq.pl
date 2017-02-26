#!/usr/local/bin/perl
##################################### USAGE ###############################

# This pipeline is to extract loci from given taxa using nhmmer and blast

#This program requires folder named as "genomes","queries" and "hmm_prof_align" in current working directory 
#It requires 3 command line arguments passed while using it
#		perl Extract_seq.pl <current working directory path> <text_file_database_list> <text_file_queries_list>
#	<current working directory path> : It is the path of current working directory which contains genomes and queries folder
#	<text_file_database_list> : It is text file that contains list of genomes in genomes folder without ".fasta". Each line should contain one genome name. It should be located in current working directory
#	<text_file_queries_list>  : It is a text file that contains list of queries in the query folder without ".fasta". Each line should contain one query name.It should be located in current working directory
#	genomes folder : It contains all the genome files after processing them using blast makedb  
#	queries folder : It contains all the query files in .fasta format
#	hmm_prof_align folder : It will store the alignment files in ".fasta" format(used by nhmmer) for all the loci
#	Following directories are created in CWD as part of the program:
#	gallus_blast: stores blast result with the given loci and gallus genome
#	other_blast: stores blast result of given loci and all genomes
#	extract_seq: stores extracted sequence and scaffold of best hit
#	blast_again: stores blast result of extracted sequence and gallus genome


use strict;
use warnings;
use List::Util qw[min max];

############# Declaration of CWD , GENOME and LOCI arrays###############

if ($#ARGV != 2 ) {
	print "usage: perl Extract_seq.pl <current working directory path> <text_file_database_list> <text_file_queries_list> \n";
	exit;
}

my($dir,$input_db,$input_query,$flank) = @ARGV;
print "The arguments passed are :\n";
foreach my $temp(@ARGV) #prints arguments
{
print"$temp\n";
}

my @db =();
my @query =();

open(my $fh1,"<","$dir/$input_db") or die "Can't open $input_db\n"; #creates db array from file provided
{
while(my $line = <$fh1>)
{
chomp($line);
push(@db,$line);
}
}
open(my $fh2,"<","$dir/$input_query") or die "Can't open $input_query\n"; #creates query array from the file provided
{
while(my $line = <$fh2>)
{
chomp($line);
push(@query,$line);
}
}

print" The list of genomes contain following names:\n"; #prints names of the genomes in the file
foreach my $genome(@db)
{
print "$genome\n";
}

print" The list of queries contain following names:\n"; #prints names of the queries in the file
foreach my $qu1(@query)
{
print"$qu1\n";
}


##############CREATING DIRECTORIES TO STORE RESULTS############################

mkdir("$dir/gallus_blast") or die "Couldn't create the gallus_blast directory"; # stores blast result with the given loci and gallus genome
mkdir("$dir/other_blast") or die "Couldn't create the other_blast directory"; # stores blast result of given loci and all genomes
mkdir("$dir/extract_seq") or die "Couldn't create the extract_seq directory";# stores extracted sequence and scaffold of best hit
mkdir("$dir/blast_again") or die "Couldn't create the blast_again directory"; # stores blast result of extracted sequence and gallus genome 
 
 ############### MAIN PROGRAM##################################################
foreach my $query(@query)
{
print"Performing Blast on $query with Gallus_gallus database and retrieving the files\n";
system("blastn -task blastn -db $dir/genome/Gallus_gallus.fasta -query $dir/queries/$query.fasta -outfmt 6 -out $dir/gallus_blast/$query.gallus.txt");
system("sed \"s/#/_/g\" $dir/gallus_blast/$query.gallus.txt > $dir/gallus_blast/temp");
unlink("$dir/gallus_blast/$query.gallus.txt") or die "couldn't delete this file!";
rename("$dir/gallus_blast/temp" ,"$dir/gallus_blast/$query.gallus.txt"); # "$query.gallus.txt" this file stores blast result of gallus genome and given loci

foreach my $db(@db)
	{
	print"Performing Blast on $query with $db database and retrieving the files\n";
	system("blastn -task blastn -db $dir/genome/$db.fasta -query $dir/queries/$query.fasta -outfmt 6 -out $dir/other_blast/$query.$db.txt");
	system("sed \"s/#/_/g\" $dir/other_blast/$query.$db.txt > $dir/other_blast/temp1");
	unlink("$dir/other_blast/$query.$db.txt") or die "couldn't delete this file!";
	rename("$dir/other_blast/temp1" ,"$dir/other_blast/$query.$db.txt"); # "$query.$db.txt"this file stores blast between given loci and other genomes ($db here)
	
	my $file = "$query.$db.txt";
	open(my $IN,"$dir/other_blast/$file") or die "Couldn't open file $file\n";
	my @IN = <$IN>;
	close($IN);
	chomp(@IN);
	my @out_db_subject = ();
	
	foreach my $lines (@IN)
	{
	  my @blast_res = split_line($lines);         # calling a function split_line that will help in parsing the blast output into different values
	  push(@out_db_subject,$blast_res[0]);
	}
	
	my $scaffold = $out_db_subject[0];    # the scaffold of the top hit is used here as the main scaffold.
	print"The scaffold from $db used is :$scaffold\n";
	print"The query file used is : $query.fasta\n";

	system("sed -n \"/^>$scaffold/,/^>/p\"  $dir/genome/$db.fasta > $dir/extract_seq/$query.$db.scaffold.txt");    #the main scaffold is then extracted from genome database file and stored as separate file
	system("egrep -v \">|>\" $dir/extract_seq/$query.$db.scaffold.txt > $dir/extract_seq/$query.$db.scaffold01.txt"); # It is to get that scaffold sequence without fasta description
	my $file1 = "$query.$db.scaffold01.txt";
	my $file2 = "$query.$db.scaffold_final.txt";
	
	open(my $out,"$dir/extract_seq/$file1") or die "Couldn't open file $file1\n"; 
	my @out = <$out>;
	chomp(@out);
	close($out);
	
	my $sub = join('',@out);
	my $seq;
	my $rev_seq;
	my $len_new;
	
	open(my $new_scaff ,">>","$dir/extract_seq/$file2") or die "Couldn't open file $file2\n";  #Created a new database or scaffold for HMMER, to use HMMER needs a scaffold with description label 
	print $new_scaff ">$db.$scaffold\n";
	print $new_scaff "$sub\n";
	close($new_scaff);
	
	##################Use HMMER to build a profile from alignment and then use that profile to get a tabular output#############################
	system("hmmbuild $dir/hmm_prof_align/$query.hmm $dir/hmm_prof_align/$query.fasta");
	system("nhmmer --dfamtblout $dir/extract_seq/$query.$db.tbl $dir/hmm_prof_align/$query.hmm  $dir/extract_seq/$query.$db.scaffold_final.txt");
	system("egrep -v \"\#\" $dir/extract_seq/$query.$db.tbl > $dir/extract_seq/test"); #This step is to remove unnecessary comment lines from the tabular output and store it in a different place 
	system("sed -n '1 p' $dir/extract_seq/test > $dir/extract_seq/test1");
	unlink("$dir/extract_seq/test") or die "couldn't delete this file!";
	unlink("$dir/extract_seq/$query.$db.tbl") or die "couldn't delete this file!";
	rename("$dir/extract_seq/test1" ,"$dir/extract_seq/$query.$db.tbl");
	
	open(my $hmm_tab,"$dir/extract_seq/$query.$db.tbl") or die "Couldn't open file $query.$db.tbl\n";
	my @hmm_res = <$hmm_tab>;
	chomp(@hmm_res);
	close($hmm_tab);
	
	my @start;
	my @end;
	my @subject_name;
	my @e_val;
	my $len_scaffold;
	my $len_new;
	my $use_start;
	my $use_end;
	my $length;
	
	foreach my $line_input(@hmm_res)
	{
		print"Maine line :$line_input\n";
		chomp($line_input);
		my @out_tab = Split_Hmm_Out($line_input);
		if ($out_tab[3] < (1e-10))
		{
		push(@start,$out_tab[0]);
		push(@end,$out_tab[1]);
		push(@e_val,$out_tab[3]);
		}
		else
		{
		next;
		}
	}
	
		if($start[0] < $end[0])
		{
		$use_start = min(@start);
		$use_end = max(@end);

		print"Best start : $use_start\n";
		print"Best end : $use_end\n";
		$len_new = ($use_end - $use_start);
		$seq = substr($sub,$use_start,$len_new);
		
		$len_scaffold = length($sub);
		$length = length($seq);
		}
		elsif($start[0] > $end[0])
		{
		print"This is a reverse hit\n";
		$use_start = max(@start);
		$use_end = min(@end);
		print"Before swapping start : $use_start \t end: $use_end\n";
	
		#swapping start and stop
		($use_start,$use_end) = ($use_end,$use_start);
		print"After start : $use_start\n";
		print"After end : $use_end\n";
	
		$len_new = ($use_end - $use_start);
		my $rev_seq = substr($sub,$use_start,$len_new);
		$seq  = reverse($rev_seq); #reverse compliment of the strand
		$seq =~ tr/ACGT/TGCA/;
		print"Reverse complimented the given seq\n";
		$len_scaffold = length($sub);
		$length = length($seq);
		}
	
	
	print"Total Length of Scaffold = $len_scaffold\n";
	print"Seq starting $use_start\n";
	print"Seq ending $use_end\n";
	print"length of extracted sequence : $length\n";

	my $extract_seq = "$dir/extract_seq/$db.$query.extract.fasta"; # "$db.$query.extract.fasta" this file stores the extracted sequence from the genome
	my $label = $db . "_genome.$query";
	open(my $two,">>","$extract_seq");
	print $two ">$label\n";
	print $two "$seq\n" ;
	close($two);
	
	open(my $fh,">>","$dir/extract_seq/log.txt"); # file "log.txt" act as  log and store the results for all hits
		print $fh "The query file opened is : $file\n";
		print $fh "The query name is :$query.fasta\n";
		print $fh "The database used is : $db\n";
		print $fh "The start position is : $use_start\n";
		print $fh "The end position is : $use_end\n";
		print $fh "The Total len of scafflod  used is  = $len_scaffold\n";
		print $fh "The total length of extracted query is : $length\n";
		print $fh "\n";

print" Performing blast on extracted sequnces with gallus_gallus database\n";
system("blastn -task blastn -db $dir/genome/Gallus_gallus.fasta -query $dir/extract_seq/$db.$query.extract.fasta -outfmt 7 -out $dir/blast_again/$query.$db.again.txt");
system("sed -n '6 p' $dir/blast_again/$query.$db.again.txt > $dir/blast_again/temp");
unlink("$dir/blast_again/$query.$db.again.txt") or die "couldn't delete this file!";
rename("$dir/blast_again/temp" ,"$dir/blast_again/$query.$db.again.txt");

my $loc = "$dir/blast_again/$query.$db.again.txt";
my @data_again_blast =location($loc);
my $gal_loc = "$dir/gallus_blast/$query.gallus.txt";
my @data_gallus = location($gal_loc);
if ("$data_again_blast[0]" eq "$data_gallus[0]")
{
	open(my $fh1,">>","compile_results.txt");
		{
			print $fh1  "$db \t $query \t $data_gallus[0] \t $data_gallus[1] \t $data_again_blast[1] \t $data_gallus[2] \t $data_again_blast[2]\n";
			close($fh1);
		}
}
else
{
		open(my $fh2,">>","mismatched_results.txt");
		{
			print $fh2 "$db \t $query \t $data_again_blast[0] \t $data_gallus[0]\n";
			close($fh2);
		}
}
}
}

sub split_line
{
	my $line = shift;
	chomp($line);
	$line =~ tr/\n/\t/;
	my @parse = split(/\t/,$line);
    my $len= scalar(@parse);
	my @query_name = ();
	my @subject_name = ();
	my @per_idn = ();
	my @align_len = ();
	my @num_mismatch = ();
	my @num_gap_pos = ();
	my @query_start = ();
	my @query_end = ();
	my @subject_start = ();
	my @subject_end = ();
	my @e_val = ();
	my @bit_score = ();
	

	for(my $i=0;$i<=12;$i++)
	{
		for(my $j=$i;$j<=$len;($j=($j+12)))
		{
		if($i==0)
		{
		push(@query_name,$parse[$j]);
		}
		elsif($i==1)
		{
		push(@subject_name,$parse[$j]);
		}
		elsif($i==2)
		{
		push(@per_idn,$parse[$j]);
		}
		elsif($i==3)
		{
		push(@align_len,$parse[$j]);
		}
		elsif($i==4)
		{
		push(@num_mismatch,$parse[$j]);
		}
		elsif($i==5)
		{
		push(@num_gap_pos,$parse[$j]);
		}
		elsif($i==6)
		{
		push(@query_start,$parse[$j]);
		}
		elsif($i==7)
		{
		push(@query_end,$parse[$j]);
		}
		elsif($i==8)
		{
		push(@subject_start,$parse[$j]);
		}
		elsif($i==9)
		{
		push(@subject_end,$parse[$j]);
		}
		elsif($i==10)
		{
		push(@e_val,$parse[$j]);
		}
		elsif($i==11)
		{
		push(@bit_score,$parse[$j]);
		}
		}
	}
	my $db_sub_name = $subject_name[0];
	my $db_start = $subject_start[0];
	my $db_end = $subject_end[0];
	my $e_val_given = $e_val[0];
	
	return $db_sub_name, $db_start, $db_end,$e_val_given;
}

sub location
{
my $cur_file = shift;
my $fp;
	open($fp,"$cur_file");
	my @fp = <$fp>;
	my $line = join('',@fp);
	$line =~ tr/\n/\t/;
	my @parse = split(/\t/,$line);
    my $len= scalar(@parse);
	my @query_name = ();
	my @subject_name = ();
	my @per_idn = ();
	my @align_len = ();
	my @num_mismatch = ();
	my @num_gap_pos = ();
	my @query_start = ();
	my @query_end = ();
	my @subject_start = ();
	my @subject_end = ();
	my @e_val = ();
	my @bit_score = ();
	

	for(my $i=0;$i<=12;$i++)
	{
		for(my $j=$i;$j<=$len;($j=($j+12)))
		{
		if($i==0)
		{
		push(@query_name,$parse[$j]);
		}
		elsif($i==1)
		{
		push(@subject_name,$parse[$j]);
		}
		elsif($i==2)
		{
		push(@per_idn,$parse[$j]);
		}
		elsif($i==3)
		{
		push(@align_len,$parse[$j]);
		}
		elsif($i==4)
		{
		push(@num_mismatch,$parse[$j]);
		}
		elsif($i==5)
		{
		push(@num_gap_pos,$parse[$j]);
		}
		elsif($i==6)
		{
		push(@query_start,$parse[$j]);
		}
		elsif($i==7)
		{
		push(@query_end,$parse[$j]);
		}
		elsif($i==8)
		{
		push(@subject_start,$parse[$j]);
		}
		elsif($i==9)
		{
		push(@subject_end,$parse[$j]);
		}
		elsif($i==10)
		{
		push(@e_val,$parse[$j]);
		}
		elsif($i==11)
		{
		push(@bit_score,$parse[$j]);
		}
		}
	}
	my $db_sub_name = $subject_name[0];
	my $db_start = $subject_start[0];
	my $db_end = $subject_end[0];
	
	return $db_sub_name, $db_start, $db_end;
	close($fp);
}
sub Split_Hmm_Out
{
	my $line = shift;
	chomp($line);
	my @parse = split(/\s+/,$line);
    my $len= scalar(@parse);
    my $target_name;
	my $ali_st;
	my $ali_end;
	my $e_value;
	
	for(my $i=0;$i<=14;$i++)
	{
		if($i==0)
		{
		$target_name = $parse[$i];
		}
		elsif($i==9)
		{
		$ali_st = $parse[$i];
		}
		elsif($i==10)
		{
		$ali_end = $parse[$i];
		}
		elsif($i == 4)
		{
		$e_value = $parse[$i];
		}
	}
	return $ali_st, $ali_end, $target_name, $e_value;
}
