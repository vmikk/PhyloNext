#!/usr/bin/env perl

## Script to run Biodiverse analyses
## Parameters:
# --input_bds_file = The full path to the biodiverse file to use (.bds)
# --input_bds_file = The full path to the tree file to use (.bts)
# --calcs          = A comma delimited list of the analyses to run
#                    E.g. 'calc_endemism_whole,calc_pd,calc_pe,calc_phylo_rpd1'

# The full list of `calc` names can be found on the Biodiverse indicies page:
# https://github.com/shawnlaffan/biodiverse/wiki/IndicesDevVersion
# E.g., "Subroutine: calc_endemism_central"

## Based on `run_analyses.pl` (988881a, Aug 6, 2014)
# https://github.com/NunzioKnerr/biodiverse_pipeline/blob/master/perl/run_analyses.pl



use strict;
use warnings;
use Carp;        #  warnings and dropouts
use File::Spec;  #  for the cat_file sub
use English qw ( -no_match_vars );

use Biodiverse::BaseData;
use Biodiverse::Tree;

use rlib;

#  load up the user defined libs
use Biodiverse::Config;

use Getopt::Long::Descriptive;

$| = 1;

#############################################
######       SET PARAMETERS HERE       ######

my ($opt, $usage) = describe_options(
  '%c <arguments>',
  [ 'input_bds_file=s',   'The input basedata file .bds', { required => 1 } ],
  [ 'input_bts_file=s',   'The input tree file .bts', { required => 1 } ],
  [ 'calcs=s',  'The comma delimited list of calculations to run', { required => 1 }],
  [ 'output_bds_file:s',  'The output biodiverse file with anaylsed .bds'],
  [],
  [ 'help',       "print usage message and exit" ],
);

 
if ($opt->help) {
    print($usage->text);
    exit;
}

print "Running diversity analysis with Biodiverse\n";

my $input_bds_file = $opt->input_bds_file;
my $input_bts_file = $opt->input_bts_file;
my $calcs          = $opt->calcs;

my $output_bds_file = $opt->output_bds_file;
if (!defined $output_bds_file) {
	$output_bds_file = $input_bds_file;
	$output_bds_file =~ s/\.bds$/_analysed\.bds/;
}


#say 'ARGS:  ' . join '::', @ARGV;


#  set to 1 to retain outputs after calculation
my $retain_outputs = 1;

#  Set to 1 to save basedata after analysis.
#  This makes sense primarily if the outputs are retained,
#  and can thus be saved with the basedata
my $save_basedata = 1;

### ANALYSIS PARAMETERS
my $do_spatial = 1;  # (do the spatial analysis: 1 for yes, 0 for no)
my $spatial_output_prefix = 'sp_';

# Define neighbour sets 1 and 2
my $spatial_conditions = ['sp_self_only()'];
#my $analyses_to_run    = [qw /calc_pd calc_pe calc_endemism_whole/];


###
$calcs =~ s/\s//g;
my @calc_list = split /,/, $calcs;
my $analyses_to_run    =  [@calc_list];
print "Calculations are " . join (q{ }, @calc_list) . "\n";
###


# (do the matrix (and cluster) analyses: 1 for yes, 0 for no)
my $do_matrix_cluster = 0;
my $matrix_cluster_output_prexif = 'matrix_';

# Define spatial conditions or subsampling for the cluster analysis
my $spatial_conditions_clust = ['sp_select_all ()'];

# (create a spatial index: 1 for yes, 0 for no)
#       this is reccomended for spatial analyses
#       and for cluster analyses with spatial conditions
#       but must not be used with matrix subsampling
my $do_build_index = 1;


######        END OF PARAMETERS        ######
#############################################



#############################################
######          LOAD THE DATA          ######

###  read in the trees from the nexus file

###  read in the basedata object
my $bd;
eval {
	$bd = Biodiverse::BaseData->new (file => $input_bds_file);
};
croak $@ if $@;

my $tree;
eval {
	$tree = Biodiverse::Tree->new(file => $input_bts_file);
};
croak $@ if $@;


#############################################
######         RUN THE ANALYSES        ######

if ($do_build_index) {
    #  build the spatial index
    $bd->build_spatial_index (
        resolutions => $bd->get_groups_ref->get_param ('CELL_SIZES')
    );
};

###### SPATIAL ANALYSIS
if ($do_spatial) {
    #  loop over the trees and add a spatial analysis to the basedata for each of the trees
    #  assuming the trees have unique names, which I think is reasonable

    ###  add a spatial analysis for each tree but use the same spatial params for each
    #  the first will take longer then the rest as they can recycle the neighbourhoods
    #  and thus save search times

    foreach my $tree_ref ($tree) {  #  dirty hack - should just use directly since we have only one tree
        #my $name = $spatial_output_prefix . $tree_ref->get_param ('NAME');
		my $name = $input_bds_file;
        my $output = $bd->add_spatial_output (name => $name);
        my $success = eval {
            $output->run_analysis (
                spatial_conditions => $spatial_conditions,
                calculations       => $analyses_to_run,
                tree_ref           => $tree_ref,
            );
        };
        croak $EVAL_ERROR if $EVAL_ERROR;

        if ($success) {  #  export to CSV using defaults
            $output->export (
                file   => $name . '.csv',
                format => 'Delimited text',
                list   => 'SPATIAL_RESULTS',
            );
			#$output->export ( # export asci grids
			#   file   => $name . '.asc',
			#   format => 'ArcInfo asciigrid files',
			#	list   => 'SPATIAL_RESULTS',
            #);
        }

        if (not $retain_outputs) {
            $bd->delete_output (output => $output);
        }
    }
};

###### MATRIX AND CLUSTER ANALYSES
###  add a cluster analysis for each tree, using the same spatial params for each

if ($do_matrix_cluster) {

    foreach my $tree_ref ($tree) {
        my $name_clust = $matrix_cluster_output_prexif . $tree_ref->get_param ('NAME');
        my $output = $bd->add_cluster_output (name => $name_clust);
        my $success = eval {
            $output->run_analysis (
                spatial_conditions   => $spatial_conditions_clust,
                #spatial_calculations => $analyses_to_run,  #  run these analyses for each node on the tree - comment out if not needed
                tree_ref             => $tree_ref,
                index                => 'S2',  #  this isn't publicly available yet (2011-01-30)
                linkage_function     => 'link_average',
            );
        };

        if ($success) {
             to export the cluster tree
            $output->export (file => $output->get_param ('NAME') . '.nex',
                              format => 'Nexus');  #  export cluster tree to Nexus using defaults

            #  export sparse matrix format to CSV
            $output->export (
                file   => $output->get_param ('NAME') . '.csv',
                format => 'Matrices', type => 'sparse',
            );

        }
        if (not $retain_outputs) {
            $bd->delete_output (output => $output);
        }
    }
};

if ($save_basedata) {
    $bd->save (filename => $output_bds_file);
};
