#!/usr/bin/perl -w

use strict;
use warnings;
use Switch;

die "Usage: perl $0  <aligner_output:blast m6 format>  <nr_length>  <query_length>" if ( @ARGV != 3 );

my ( $align, $nr_length, $query_length ) = @ARGV;

my %hsp_h          = ();
my %hit_identifier = ();
my %hit_h          = ();
my %query_h        = ();

load_query_table( );
read_align( );
load_nr_table( );
print_xml( );


sub load_query_table {

    open( DATA,    $query_length )   || die "$!";

    while (<DATA>) {
        chomp;
        my @its  = split /\t/, $_;
        $query_h{ $its[0] } = \@its;
    }

    close DATA;
}

sub read_align{

    open( DATA,    $align )   || die "$!";
    
    my %hsp_flags = ();
    my %hit_flags = ();

    while (<DATA>) {
        
        chomp;
        next if(/^@/);
       
        my @its  = split /\t/, $_;
       
        $hit_flags{ $its[0] }++;
        $hsp_flags{ $its[0] }{ $its[1] } ++;

		    $hit_identifier{$its[1]} = ();

        #Hsp_num bit-score evalue query-from query-to hit-from hit-to identity positive gaps align-len
        
        my $identity =  int($its[2] * $its[3]/100);
        my $positive =  $its[3] - $its[4] - $its[5]; 

        my @metrics  = ($hsp_flags{ $its[0] }{ $its[1] }, $its[11], $its[10], $its[6], $its[7], $its[8], $its[9], $identity, $identity, $its[5], $its[3]); 

        $hsp_h{ $its[0] }{ $hit_flags{ $its[0] } }{hit} = $its[1];
        $hsp_h{ $its[0] }{ $hit_flags{ $its[0] } }{hsp}{ $hsp_flags{ $its[0] }{ $its[1] } } =  \@metrics;

    }
    close DATA;
}


sub load_nr_table {
    
    open( DATA,    $nr_length )   || die "$!";

    while (<DATA>) {
        chomp;
        my @its  = split /\t/, $_;
        next if(!exists $hit_identifier{ $its[0] });
        
        my $hit_accession = hit_accession($its[0]);

        my @tuple = ($its[0], $its[1], $hit_accession, $its[2]);

        $hit_h{ $its[0] } = \@tuple;
    }

    close DATA;

}


sub hit_accession{

	my $val  = shift;
	my @vals = split /\|/, $val;
	return $vals[-1];

}


sub print_xml {

    &print_head();

    foreach my $x (keys %hsp_h) {
        print_query_length( $query_h{$x} );
        
        flush_hits($x);

        print_query_tail();
    }

    &print_tail();

}


sub print_query_length{
  
  my $val = shift;

  printf qq{<Iteration>
  <Iteration_iter-num>-1</Iteration_iter-num>
  <Iteration_query-ID>%s</Iteration_query-ID>
  <Iteration_query-def>%s</Iteration_query-def>
  <Iteration_query-len>%d</Iteration_query-len>
<Iteration_hits>
} , @{$val};

}

sub print_query_tail{

  print qq{</Iteration_hits>\n</Iteration>\n};

}


sub flush_hits {

    my $query  = shift;
  
    foreach my $hit (sort {$a <=> $b} keys %{$hsp_h{$query}}) {

        print_hit_head($hit, $hit_h{ $hsp_h{$query}{$hit}{hit}});

        foreach my $hsp ( sort {$a <=> $b} keys %{$hsp_h{$query}{$hit}{hsp} } ) {
            print print_hsp( $hsp_h{$query}{$hit}{hsp}{$hsp} ); 
        }

        print_hit_tail();
    }
}


sub print_hit_head{

 my ($ind, $val) = @_;
 printf qq{<Hit>
  <Hit_num>%d</Hit_num>
  <Hit_id>%s</Hit_id>
  <Hit_def>%s</Hit_def>
  <Hit_accession>%s</Hit_accession>
  <Hit_len>%d</Hit_len>
  <Hit_hsps>}, ($ind, @{$val});

}


sub print_hit_tail{

  print qq{  </Hit_hsps>\n</Hit>\n};

}


sub print_hsp{

  my $val = shift;
  printf qq{
    <Hsp>
      <Hsp_num>%s</Hsp_num>
      <Hsp_bit-score>%s</Hsp_bit-score>
      <Hsp_score>-1</Hsp_score>
      <Hsp_evalue>%s</Hsp_evalue>
      <Hsp_query-from>%s</Hsp_query-from>
      <Hsp_query-to>%s</Hsp_query-to>
      <Hsp_hit-from>%s</Hsp_hit-from>
      <Hsp_hit-to>%s</Hsp_hit-to>
      <Hsp_query-frame>-1</Hsp_query-frame>
      <Hsp_hit-frame>-1</Hsp_hit-frame>
      <Hsp_identity>%s</Hsp_identity>
      <Hsp_positive>%s</Hsp_positive>
      <Hsp_gaps>%s</Hsp_gaps>
      <Hsp_align-len>%s</Hsp_align-len>
      <Hsp_qseq>na</Hsp_qseq>
      <Hsp_hseq>na</Hsp_hseq>
      <Hsp_midline>na</Hsp_midline>
    </Hsp>}, @{$val};

}


sub print_head{

  print qq{<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.2.30+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>1</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
<BlastOutput_iterations>
}; 

}

sub print_tail{

   print qq{</BlastOutput_iterations>\n</BlastOutput>\n}

}