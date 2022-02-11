#!/usr/bin/env perl


## Script to convert results into csv format.
## Parameters:
# --input_bds_file     = The biodiverse file (.bds) to load and export results from
# --output_csv_prefix  = The output file prefix

## based on `load_bds_and_export_results.pl` script (69b2bd2, May 28, 2017)
# https://github.com/NunzioKnerr/biodiverse_pipeline/blob/master/perl/load_bds_and_export_results.pl


use strict;
use warnings;

use Carp;        #  warnings and dropouts
use File::Spec;  #  for the cat_file sub
use English qw ( -no_match_vars );
use feature qw(say);

use Biodiverse::BaseData;


#  load up the user defined libs
use Biodiverse::Config;



#############################################
######       SET PARAMETERS HERE       ######


use Getopt::Long::Descriptive;

$| = 1;

#############################################
######       SET PARAMETERS HERE       ######

my ($opt, $usage) = describe_options(
  '%c <arguments>',
  [ 'input_bds_file=s',     'The input basedata file .bds', { required => 1 } ],
  [ 'output_csv_prefix:s',  'The output biodiverse results file prefix (.csv is appended)'],
  [],
  [ 'help',       "print usage message and exit" ],
);

 
if ($opt->help) {
    print($usage->text);
    exit;
}

my $input_bds_file    = $opt->input_bds_file;
my $output_csv_prefix = $opt->output_csv_prefix;


###  read in the basedata object
my $bd;
eval {
    $bd = Biodiverse::BaseData->new (file => $input_bds_file);
};
croak $@ if $@;


my @outputs = $bd->get_spatial_output_refs;



foreach my $output (@outputs) {
    my @lists = $output->get_lists_across_elements;
    foreach my $list (@lists) {
        say $list;
        my $list_fname = $list;
        if ($list_fname =~ s/>>/--/g) {  #  CHANGE THIS LINE TO ALTER THE >> REPLACEMENT TEXT IN FILENAME
            say 'Output name changed - now it is ' . $list_fname;
        }
        my $csv_file = sprintf "%s_%s.csv", $output_csv_prefix, $list_fname;
        $output->export (
            file   =>  $csv_file,
            format => 'Delimited text',
            list   => $list,
        );
        
    }
};

my $csv_file = sprintf "%s_%s.csv", $output_csv_prefix, "groups";
$bd->get_groups_ref->export (
    file   =>  $csv_file,
    format => 'Delimited text',
    list => 'SUBELEMENTS',
);
say 'Process completed';

