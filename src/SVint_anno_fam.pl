# read in family info: parents/probands/other relavtives; genders; affected or not
# read in SV

#! /h/jingqichen/localperl/bin/perl
use warnings;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require SV_SUB;

$dir = "";
$fam_file = ""; # id, parents/probands/relatives, gender, affected
$sv_loc_file = ""; # id, sv file location

($proband, $father, $mother, $id2gender, $id2status, $id2svfile) = SV_SUB::read_fam($fam_file);
($same_a, $diff_a, $lines_a) = SV_SUB::read_sv($proband); # %{$all_a} all sv for proband
($same_b, $diff_b, $lines_b) = SV_SUB::read_sv($mother); # %{$all_b} all sv for mother
($same_c, $diff_c, $lines_c) = SV_SUB::read_sv($father); # %{$all_c} all sv for father


