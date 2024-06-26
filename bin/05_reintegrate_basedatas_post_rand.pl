#!/usr/bin/env perl

## Script to reintegrate multiple BDS files after randomizations
## Author: Shawn Laffan
## Based on 4adf0f0 commit (Nov 29, 2019)
## https://raw.githubusercontent.com/shawnlaffan/biodiverse/73522b74e52a5fb77ae5cfc2e90010350a3abf70/bin/reintegrate_basedatas_post_rand.pl

use 5.016;

use rlib;

use Biodiverse::BaseData;

use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options(
  '%c <arguments>',

  [ 'glob=s',               'glob to find files',                     { required => 1 } ],
  [ 'output_prefix|opfx=s', 'The output prefix for exported files',   { required => 1 } ],
  [ 'no_verify|n!',         'do not verify that the basedatas match', { default  => 1 } ],

  [ 'help',                 "print usage message and exit" ],
);


print "Reintegrating randomization chunks\n";

my $glob = $opt->glob;
my $opfx = $opt->output_prefix;
my $no_verify = $opt->no_verify;

print "\nInput file with paths to randomization chunks:\n";
print "$glob\n";

## Read paths to randomization chunks
open my $glob_file, '<', $glob or die $!;
chomp( my @files = <$glob_file> );
close $glob_file;

print "\nList of randomization chunks:\n";
print "@files\n\n";

print "Multiple randomization chunks fond, reinegrating..\n" if !@files;
print "Just a single randomization file found\n" if @files == 1;

my $first_bd = shift @files;
my $bd = Biodiverse::BaseData->new (file => $first_bd);

## Reintegrate chunks
if (scalar(@files) > 1) {
    foreach my $from_file (@files) {
        say "Reintegrating from $from_file";

        my $from_bd = Biodiverse::BaseData->new(file => $from_file);
        $bd->reintegrate_after_parallel_randomisations (
            from => $from_bd,
            no_check_groups_and_labels => $no_verify,
        );
    }
}

$bd->save_to (filename => $opfx, method => 'save_to_sereal');

