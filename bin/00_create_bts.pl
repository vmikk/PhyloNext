#!/usr/bin/env perl

## Script to prepare phylogenetic tree for Biodiverse
## Parameters:
# --input_tree_file  = Phylogenetic tree in NEXUS format
# --out_file         = Output file (bts)

## Based on `create_bts.pl` (988881a, Aug 6, 2014)
# https://github.com/NunzioKnerr/biodiverse_pipeline/blob/master/perl/create_bds.pl


use strict;
use warnings;
use Carp;        #  warnings and dropouts
use File::Spec;  #  for the cat_file sub
use English qw ( -no_match_vars );

use Biodiverse::BaseData;
use Biodiverse::ElementProperties;  #  for remaps
use Biodiverse::ReadNexus;
use Biodiverse::Tree;


#############################################
######       SET PARAMETERS HERE       ######


use Getopt::Long::Descriptive;

$| = 1;

my ($opt, $usage) = describe_options(
  '%c <arguments>',
  [ 'input_tree_file=s',   'The input tree file .nex', { required => 1 } ],
  [ 'out_file=s',  'The output biodiverse tree file .bts', { required => 1 }],
  [],
  [ 'help',       "print usage message and exit" ],
);

 
if ($opt->help) {
    print($usage->text);
    exit;
}

my $nexus_file    = $opt->input_tree_file;
my $tree_out_file = $opt->out_file;



#############################################
######          LOAD THE DATA          ######

###  read in the trees from the nexus file

#  but first specify the remap (element properties) table to use to match the tree names to the basedata object
#my $remap = Biodiverse::ElementProperties->new;
#$remap->import_data (
#    #file => $nex_remap,
#    input_sep_char        => 'guess',   #  or specify yourself, eg ',' for a comma
#    input_quote_char      => 'guess', #  or specify yourself, eg:'"' or "'" (for double or single quotes)
#    #input_element_cols    => @remap_input_element_cols,  # defined above
#    #remapped_element_cols => @remapped_element_cols,  #  defined above
#    #include_cols          => undef,  #  undef or empty array [] if none are to be solely included
#    #exclude_cols          => undef,  #  undef or empty array [] if none are to be excluded
#);


#  read the nexus file
my $read_nex = Biodiverse::ReadNexus->new;
$read_nex->import_data (
    file                   => $nexus_file,
    use_element_properties => 0,  #  set to zero or undef if you don't have a remap
    #element_properties     => $remap,
);

#  get an array of the trees contained in the nexus file
my @trees = $read_nex->get_tree_array;

#  just a little feedback
my $tree_count = scalar @trees;
print "$tree_count trees parsed from $nexus_file\nNames are: ";
my @names;
foreach my $tree (@trees) {
    push @names, $tree->get_param ('NAME');
}
print join (q{, }, @names), "\n";

#  only one tree, so get the first
my $tree = shift @trees;

$tree->save (filename => $tree_out_file);
