#!/usr/bin/env perl

## Script to create a BaseData for Biodiverse
## Parameters:
# --csv_file  = File in CSV format
# --out_file  = Output file (bds)

## Based on `create_bds.pl` (988881a, Aug 6, 2014)
# https://github.com/NunzioKnerr/biodiverse_pipeline/blob/master/perl/create_bds.pl

## Input data format:
# Column 1 = species name
# Column 2 = H3 grid cell ID
# Column 5 = Latitude
# Column 6 = Longitude
# Column 7 = Number of records (per species, per grid cell)

use strict;
use warnings;
use Biodiverse::BaseData;
use Biodiverse::Common;
use Biodiverse::ElementProperties;

use Data::Dumper qw/Dumper/;

use Getopt::Long::Descriptive;

$| = 1;


## NB! column numbers are zero-based!
# label = species name
# group = grid cell

my ($opt, $usage) = describe_options(
  '%c <arguments>',
  
  [ 'csv_file=s',              'The input csv file',       { required => 1 } ],
  [ 'out_file|output_bd=s',    'The output basedata file', { required => 1 } ],

  [ 'label_column_number:i',     'Column containing the label name [default= "0"]',    { default => 0 } ],
  [ 'sampcount_column_number:i', 'Column containing the sample counts [default= "6"]', { default => 6 } ],
# [ 'group_column_text:i',       'Column containing the group name [default= "1"]',    { default => 1 } ],
  [ 'group_column_number_x:i',   'Column containing the x-axis values [default= "2"]', { default => 4 } ],
  [ 'group_column_number_y:i',   'Column containing the y-axis values [default= "3"]', { default => 5 } ],
  
  [ 'cell_size_x:f',           'Cell size of x-axis [default= "0"]', { default => 0 } ],
  [ 'cell_size_y:f',           'Cell size of y-axis [default= "0"]', { default => 0 } ],

  [ 'help',       "print usage message and exit" ],
);

 
if ($opt->help) {
    print($usage->text);
    exit;
}

print "Preparing occurrence data for Biodiverse\n";

my $csv_file              = $opt->csv_file;
my $out_file              = $opt->out_file;

my $label_column_number   = $opt->label_column_number;

my $sampcount_column_number = $opt->sampcount_column_number;


# my $group_column_text   = $opt->group_column_text;
my $group_column_number_x = $opt->group_column_number_x;
my $group_column_number_y = $opt->group_column_number_y;

my $cell_size_x           = $opt->cell_size_x;
my $cell_size_y           = $opt->cell_size_y;


my @table;

my $double_quotes = q{"};

my $dist = 'all';

my @names = (
	$csv_file,
);

my %input_files = (
    $names[0] => {
        input_quotes         => $double_quotes,
        label_columns        => [$label_column_number],
        group_columns        => [$group_column_number_x, $group_column_number_y],
        # group_columns      => [$group_column_number_x],
        sample_count_columns => [$sampcount_column_number],
    }
);

my $bd;
eval {
    $bd = Biodiverse::BaseData -> new (file => $out_file . '.bds');
};
if (not defined $bd) {
    $bd = create_bd();
    $bd -> save;
}


sub create_bd {
    my $bd = Biodiverse::BaseData -> new (
        NAME        => $out_file,
        CELL_SIZES  => [$cell_size_x, $cell_size_y],
        # CELL_SIZES  => [$cell_size_x],
    );

    foreach my $file (@names) {
        my %args = %{$input_files{$file}};
        $bd -> load_data (
            input_sep_char        => ",",
            include_columns       => [],
            %args,
            input_quotes          => $double_quotes,
            input_files           => [$file],
            use_label_properties  => 0,
            #label_properties      => $lb_remap,
        );
    }

    return $bd;
}

