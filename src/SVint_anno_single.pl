# read in family info: parents/probands/other relavtives; genders; affected or not
# read in SV

#! /h/jingqichen/localperl/bin/perl
use warnings;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require SV_SUB;

$sv_file = $ARGV[0]; # $hg=$ARGV[1];
$gene_anno = "data/gene_anno/Homo_sapiens.GRCh37.87.PC.genes.bed"; # sorted
$re_anno = "data/reg_anno/combined_regulatory_elements_for_SVint.bed"; #sorted
$rg_HIC = "data/target/INTREPID/INTREPID-HIC_RE2target_PCgenes.txt";
$rg_PHY = "data/target/INTREPID/INTREPID-PHY_RE2target_PCgenes.txt";
$rg_IMPET = "data/target/impet/IMPET_RE2target.txt";
$rg_PRESTIGE = "data/target/prestige/preSTIGE_RE2target_PCgenes.txt";
$rg_TARGETFINDER = "data/target/targetfinder/TargetFinder_combined_EPW_GBM_RE2target.txt";
$rg_HICHIP = "data/target/hichip/HiChip_RE2target.txt";
$rg_EQTL = "data/target/eqtl/GTEx_singleTissue_eQTL_RE2target.txt";
$kg_fi_file = "data/kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.hg19.vcf_based_FI.txt";
print "warning: SVint is using $kg_fi_file\n";

print "Read in data...\n";
($same_a, $diff_a, $lines_a) = SV_SUB::read_sv($sv_file);
# same chr, diff chr; in this tool only SVs that don't cross different chr are used for analysis
%same_a = %{$same_a}; %diff_a = %{$diff_a}; %lines_a = %{$lines_a};
# $same_a{$var}{$type}=$gt; $diff_a{"$loc1\t$loc2"}{$type}=$gt; $lines_a{$id}{$line}=$var or "$loc1\t$loc2";

# filter: only use intergenic SVs
print "Filter by gene annotation...\n";
&SV_SUB::write_bed($same_a, "$sv_file\_temp"); # output $sv_file\_temp.sorted.bed
$overlap_gene = SV_SUB::overlap("$sv_file\_temp.sorted.bed", $gene_anno);
%overlap_gene = %{$overlap_gene};
my %see; map{$see{$_}=1}(keys %overlap_gene);
my %to_use;
for $var(keys %same_a){
	for $tp(keys %{$same_a{$var}}){
		next if exists $see{$var}; $to_use{$var}{$tp}=$same_a{$var}{$tp};
	}
}
&SV_SUB::write_bed(\%to_use, "$sv_file\_temp2"); # output $sv_file\_temp2.sorted.bed

# annotate: with regulatory elements
print "Annotate with regulatory elements...\n";
$sv2re = SV_SUB::overlap("$sv_file\_temp2.sorted.bed", $re_anno); #$sv2re{$var}{$re} = overlapped length ratio in RE; $var & $re: "$chr\t$stt\t$end"

# annotate with each regulatory element - target gene file:
print "Annotate with RE-target information...\n";
$sv2target_HIC = SV_SUB::rg_anno($sv2re, $re_anno, $rg_HIC);
$sv2target_PHY = SV_SUB::rg_anno($sv2re, $re_anno, $rg_PHY);
$sv2target_IMPET = SV_SUB::rg_anno($sv2re, $re_anno, $rg_IMPET);
$sv2target_PRESTIGE = SV_SUB::rg_anno($sv2re, $re_anno, $rg_PRESTIGE);
$sv2target_TARGETFINDER = SV_SUB::rg_anno($sv2re, $re_anno, $rg_TARGETFINDER);
$sv2target_HICHIP = SV_SUB::rg_anno($sv2re, $re_anno, $rg_HICHIP);
$sv2target_EQTL = SV_SUB::rg_anno($sv2re, $re_anno, $rg_EQTL);

# annotate with FIR
print "Annotate with Frequency of Impact for REs calculated through 1000genome SV data...\n";
$fi = SV_SUB::read_fi($kg_fi_file);
$sv2fi_HIC = SV_SUB::anno_fi($sv2target_HIC, \%to_use, $fi);
$sv2fi_PHY = SV_SUB::anno_fi($sv2target_PHY, \%to_use, $fi);
$sv2fi_IMPET = SV_SUB::anno_fi($sv2target_IMPET, \%to_use, $fi);
$sv2fi_PRESTIGE = SV_SUB::anno_fi($sv2target_PRESTIGE, \%to_use, $fi);
$sv2fi_TARGETFINDER = SV_SUB::anno_fi($sv2target_TARGETFINDER, \%to_use, $fi);
$sv2fi_HICHIP = SV_SUB::anno_fi($sv2target_HICHIP, \%to_use, $fi);
$sv2fi_EQTL = SV_SUB::anno_fi($sv2target_EQTL, \%to_use, $fi);

# write out new VCF
print "Write out the new VCF...\n";
undef %gs; undef %fs; undef %g2f_names; undef @order;
@order =($sv2target_HIC, $sv2target_PHY, $sv2target_IMPET, $sv2target_PRESTIGE, $sv2target_TARGETFINDER, $sv2target_HICHIP, $sv2target_EQTL);
$gs{$sv2target_HIC}=$sv2target_HIC; $fs{$sv2target_HIC}=$sv2fi_HIC;
$gs{$sv2target_PHY}=$sv2target_PHY; $fs{$sv2target_PHY}=$sv2fi_PHY;
$gs{$sv2target_IMPET}=$sv2target_IMPET; $fs{$sv2target_IMPET}=$sv2fi_IMPET;
$gs{$sv2target_PRESTIGE}=$sv2target_PRESTIGE; $fs{$sv2target_PRESTIGE}=$sv2fi_PRESTIGE;
$gs{$sv2target_TARGETFINDER}=$sv2target_TARGETFINDER; $fs{$sv2target_TARGETFINDER}=$sv2fi_TARGETFINDER;
$gs{$sv2target_HICHIP}=$sv2target_HICHIP; $fs{$sv2target_HICHIP}=$sv2fi_HICHIP;
$gs{$sv2target_EQTL}=$sv2target_EQTL; $fs{$sv2target_EQTL}=$sv2fi_EQTL;
$g2f_names{$sv2target_HIC}="INTREPID-HIC";
$g2f_names{$sv2target_PHY}="INTREPID-PHY";
$g2f_names{$sv2target_IMPET}="IMPET";
$g2f_names{$sv2target_PRESTIGE}="PreSTIGE";
$g2f_names{$sv2target_TARGETFINDER}="TargetFinder";
$g2f_names{$sv2target_HICHIP}="HiChip";
$g2f_names{$sv2target_EQTL}="GTEx-eQTL";

&SV_SUB::write_svint_anno($sv_file."_SVint_anno.vcf", $sv_file, $lines_a, \%to_use, \@order, \%gs, \%fs, \%g2f_names);
print "SVint done.\n";

unlink "$sv_file\_temp.*bed";
