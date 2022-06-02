#!/usr/bin/env perl

## Script to prepare phylogenetic tree for Biodiverse
## Parameters:
# --input_bds_file   # Basedata file (bds) produced by `00_create_bds.pl`
# --input_bts_file   # Tree file (bts) produced by `00_create_bts.pl`
# --output_bds_file  # Output trimmed bds
# --output_bts_file  # Output trimmed bts

## based on `trim_bds_and_bts.pl` script (988881a, Aug 6, 2014)
# https://github.com/NunzioKnerr/biodiverse_pipeline/blob/master/perl/trim_bds_and_bts.pl


use strict;
use warnings;
use Carp;
use English qw /-no_match_vars/;

use Data::Dumper;

use Biodiverse::BaseData;
use Biodiverse::Common;
use Biodiverse::ElementProperties;

use Data::Dumper qw/Dumper/;

use Getopt::Long::Descriptive;

$| = 1;

my ($opt, $usage) = describe_options(
  '%c <arguments>',
  [ 'input_bds_file=s',   'The input basedata file .bds', { required => 1 } ],
  [ 'input_bts_file=s',   'The input tree file .bts', { required => 1 } ],
  [ 'output_bds_file=s',  'The output biodiverse basedata file .bds', { required => 1 }],
  [ 'output_bts_file=s',  'The output biodiverse tree file .bts', { required => 1 }],
  [],
  [ 'help',       "print usage message and exit" ],
);


if ($opt->help) {
    print($usage->text);
    exit;
}

my $input_bds_file    = $opt->input_bds_file;
my $input_bts_file = $opt->input_bts_file;
my $output_bds_file  = $opt->output_bds_file;
my $output_bts_file = $opt->output_bts_file;


my $double_quotes = q{"};

my $bd;
eval {
    $bd = Biodiverse::BaseData -> new (file => $input_bds_file);
    1;
} // do {croak $EVAL_ERROR if $EVAL_ERROR};

my $tree;
eval {
    $tree = Biodiverse::Tree -> new (file => $input_bts_file);
    1;
} // do {croak $EVAL_ERROR if $EVAL_ERROR};


$bd->delete_all_outputs;

my $label_hash = $bd->get_labels;
my $trimmed_tree = $tree->trim (keep => $label_hash);


my %summary = $bd->trim (keep => $tree);

say 'BASEDATA DELETIONS:';
say join q{ }, %summary;

$bd->save_to (filename => $output_bds_file);
$tree->save_to (filename => $output_bts_file);

