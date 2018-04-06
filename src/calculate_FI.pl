# calculate Frequency of Impact for the RE (FIR)
# for INV and INS, the impact is always 0.5, since the RE dosage is most probably not changed
# if in the same sample one RE was impacted in different levels, the more severe levels should be also counted for the less severe levels: e.g. if in one sample, RE1 has 100%impact for 1 allele, then it should also have 50% and <50% impact for 1 allele, unless it already have >=1 record of other impact levels

use warnings;
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0);

require SV_SUB;

$sv_file = "data/kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.hg19.vcf";
$re_anno = "data/reg_anno/combined_regulatory_elements_for_SVint.bed"; #sorted

open F, $sv_file || die;
while(<F>){
	chomp; $line = $_;
	next if $line =~ /^##/;
	if($line =~ /^#/){
		$sample_no = scalar(split /\t/, $line) - 9; next;
	}
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  samples...
	($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = split /\t/, $_;
	next if $filter ne "PASS";
	$chr="chr$chr"; $stt = $pos - 1;
	next if $info !~ /END=/; next if $info !~ /SVTYPE=/;
	if($info =~ /END=(\d+);/){
		$end = $1;
	}else{
		$info =~ /SVLEN=(\d+);/; $end = $stt + $1;
	}
	$info =~ /SVTYPE=(.*)/; $tp=$1; # ALU;CNV;DEL;DEL_ALU;DEL_HERV;DEL_LINE1;DEL_SVA;DUP;INS;INV;LINE1;SVA
		$var = "$chr\t$stt\t$end";
# DEL
	if($tp =~ /^DEL/){
		$to_use{"DEL"}{$var}=$line;
		$to_use{"all"}{$var}=$line;	
	}elsif($tp =~ /DUP/){
# DUP / CNV
		$to_use{"DUP"}{$var}=$line;
		$to_use{"all"}{$var}=$line;
	}elsif($tp =~ /CNV/){
#CNV
		$to_use{"CNV"}{$var}=$line;
		$to_use{"all"}{$var}=$line;
	}elsif($tp =~ /INV/){
# INV
		$to_use{"INV"}{$var}=$line;
		$to_use{"all"}{$var}=$line;
	}else{
# INS, ALU, LINE1, SVA
		$to_use{"INS"}{"$chr\t$stt\t$pos"}=$line;
		$to_use{"all"}{"$chr\t$stt\t$pos"}=$line;
	}
}
close F || die;

&SV_SUB::write_bed2(\%to_use, "$sv_file\_temp");
$sv2re = SV_SUB::overlap("$sv_file\_temp.sorted.bed", $re_anno);
%sv2re = %{$sv2re};
unlink "$sv_file\_temp.sorted.bed";
unlink "$sv_file\_temp.bed";

# $..{$tp}{$sample_id}{$re}; for each sample, the same type of impact on the same RE will only be count once; for the same sample, different impact will all be taken into consideration
open O, ">$sv_file\_based_FIR.txt" || die;
print O "SVTYPE\tCHR\tSTT\tEND\tFIR(100% IMPACT)\tFIR(50-100% IMPACT)\tFIR(0-50% IMPACT)\n";
undef %count; # $count{$tp}{$re_loc}{$im}
for $tp(keys %to_use){
	print "$tp:\n";
	for $i(1 .. $sample_no){
		undef %re2count;
		print "$i\n";
		for $var(keys %{$to_use{$tp}}){
			$line = $to_use{$tp}{$var};
			@tmp = split /\t/, $line;
			$col = 9 + $i -1; $gt = $tmp[$col];
			next if $gt eq ".";
			$sum_gt = 0;                            
			for $tt(split /\|/, $gt){
				$sum_gt = $sum_gt + 1 if $tt >= 1;
			}
			for $re_loc(keys %{$sv2re{$var}}){
				$ratio = $sv2re{$var}{$re_loc};
				if($ratio >= 1){
					$im = "set3";
				}elsif($ratio < 0.5){
					$im = "set1";
				}else{
					$im = "set2";
				}
				if($tp eq "INV" || $tp eq "INS"){
					$im = "set2";
				}
				if(exists $re2count{$re_loc}{$im}){
					undef $c;
					$c = $re2count{$re_loc}{$im};
					$re2count{$re_loc}{$im}=$c>$sum_gt?$c:$sum_gt;
				}else{
					$re2count{$re_loc}{$im}=$sum_gt;
				}
			}
		}
# 100% impact (set3) will also be counted for 0-0.5 (set1) and 0.5-1 (set2); 0.5-1 (set2) impact will also be counted for 0-0.5 (set1);
# for each case, the same re_loc with each impact should at most be counted twice (diplotype) 
		for $re_loc(keys %re2count){
			$count{$tp}{$re_loc}{"set1"} = 0 if (!exists $count{$tp}{$re_loc}{"set1"});
			$count{$tp}{$re_loc}{"set2"} = 0 if (!exists $count{$tp}{$re_loc}{"set2"});
			$count{$tp}{$re_loc}{"set3"} = 0 if (!exists $count{$tp}{$re_loc}{"set3"});
			$add3 = 0; $add2 = 0; $add1 = 0;
			$add3 = $re2count{$re_loc}{"set3"} if exists $re2count{$re_loc}{"set3"};
			$add2 = $re2count{$re_loc}{"set2"} if exists $re2count{$re_loc}{"set2"}; $add2 = $add2>$add3?$add2:$add3;
			$add1 = $re2count{$re_loc}{"set1"} if exists $re2count{$re_loc}{"set1"}; $add1 = $add1>$add2?$add1:$add2;
			$count{$tp}{$re_loc}{"set3"} += $add3;
			$count{$tp}{$re_loc}{"set2"} += $add2;
			$count{$tp}{$re_loc}{"set1"} += $add1;
		}
	}
	for $re_loc(keys %{$count{$tp}}){
		$c1 = $count{$tp}{$re_loc}{"set1"};
		$c2 = $count{$tp}{$re_loc}{"set2"};
		$c3 = $count{$tp}{$re_loc}{"set3"};
		$fir1 = $c1/($sample_no * 2);
		$fir2 = $c2/($sample_no * 2);
		$fir3 = $c3/($sample_no * 2);
		print O "$tp\t$re_loc\t$fir3\t$fir2\t$fir1\n";
	}
}
close O || die;
