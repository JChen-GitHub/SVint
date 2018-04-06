#! /h/jingqichen/localperl/bin/perl

package SV_SUB;
use warnings;

sub add {
	my ($x, $y) = @_;
	return $x + $y;
}

# read in family information
# ID parents/proband gender affected/unaffected sv-file-location
sub read_fam {
	my ($file) = @_;
	my ($proband, $father, $mother, %id2gender, %id2status, %id2svfile); 
	open I, $file;
	while(<I>){
		next if $_ =~ /^#/;
		chomp; ($id, $nm, $gender, $dis_status, $svfile) = split /\t/;
		$id2gender{$id} = $gender;
		$id2status{$id} = $dis_status;
		$id2svfile{$id} = $svfile;
		if($nm eq "proband"){
			$proband = $id; next;
		}elsif($nm eq "father"){
			$father = $id; next;
		}elsif($nm eq "mother"){
			$mother = $id; next;
		}else{
			push @relatives, $id;
		}
	}
	close I || die;
	return (\$proband, \$father, \$mother, \@relatives, \%id2gender, \%id2status, \%id2svfile);
}


# read in SV files in VCF format (using 10X longranger call format as an example)
# output: same chr SV; cross-chr SV
# use only paired BND; use "0/1" if two BND disagree on GT
sub read_sv {
	my ($file) = @_;
	my (%same, %diff, %mates, %bnd_lines, %lines);
	open I, $file;
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  BC00101_longranger_noloupe_GATK
	my $count=0;	
	while(<I>){
		next if $_ =~ /^#/;
		chomp; $line=$_; ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $tgt) = split /\t/;
		$count++;
		$bnd_lines{$id}{$line} = 1; $id2count{$id}=$count; $lines{$count}{$line}=1;
		$chr=~ s/Chr/chr/; $chr = "chr$chr" if $chr !~ /chr/;
		print "warning for SV call: format not contain \"GT\"!\n" if $format !~ "GT";
		@f = split ":", $format; @tgt = split ":", $tgt;
		for $k(0 .. $#f){
			if($f[$k] eq "GT"){
				$gt = $tgt[$k]; last;
			}
		}
		if($alt =~ /^<(.*)>$/){
			$tp = $1;
			$info =~ /END=(\d+);/; $end=$1;
			@loc = sort{$a<=>$b}($pos, $end); $var = "$chr\t$loc[0]\t$loc[1]";
			$same{$var}{$tp} = $gt; $lines{$count}{$line}=$var;
		}else{
			undef $m_id; $e=0;
			@in = split ";", $info;
			for $in(@in){
				if($in =~ /MATEID=(.*)/){
					$m_id = $1; $e=1;
					last;
				}}
			if($e == 1){
				@ids = sort{$a cmp $b}($id, $m_id);
				$mates{$ids[0]}{$ids[1]}=1;
			}
		}
	}
	close I || die;
	for $id1(keys %mates){
		@line1 = keys %{$bnd_lines{$id1}}; @t1 = split /\t/, $line1[0];$count1=$id2count{$id1};
		for $id2(keys %{$mates{$id1}}){
			@line2 = keys %{$bnd_lines{$id2}}; @t2 = split /\t/, $line2[0]; $count2=$id2count{$id2};
			$chr1 = $t1[0]; $pos1 = $t1[1]; $chr2 = $t2[0]; $pos2 = $t2[1]; $gt1=$t1[$#t1]; $gt2=$t2[$#t2];
			undef $tp; undef $gt;
			if($gt1 ne $gt2){
				$gt="0/1";
			}else{
				$gt=$gt1;
			}
			@in = split ";", $t1[7];
			for $in(@in){
				if($in =~ /SVTYPE2=(.*)/){
					$tp = $1;
				}elsif($in =~ /SVTYPE=(.*)/){
					$tp = $1;
				}
			}
			if($chr1 eq $chr2){
				@loc = sort{$a<=>$b}($pos1, $pos2); $var = "$chr\t$loc[0]\t$loc[1]";
				$same{$var}{$tp}=$gt; 
$lines{$count1}{$line1[0]}=$var; $lines{$count2}{$line2[0]}=$var;
			}else{
				$var1 = "$chr1\t$pos1"; $var2 = "$chr2\t$pos2";
				@loc = sort{$a cmp $b}($var1, $var2); $var = "$loc[0]\t$loc[1]";
				$diff{$var}{$tp}=$gt;
				$lines{$count1}{$line1[0]} = $var; $lines2{$count2}{$line2[0]} = $var;
			}
		}}
	return(\%same, \%diff, \%lines);
}


#
sub write_bed{ # based on the result from read_in()
	my ($list, $bed) = @_;
	my %list = %{$list};
	open O, ">$bed.bed";
	for $var(keys %list){
		($chr, $stt, $end) = split /\t/, $var;
		print O "$chr\t$stt\t$end\n";
	}
	close O || die;
	system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed";
	return()
}

#
sub write_bed2{ # based on the result from read_in()
	my ($list, $bed) = @_;
	my %list = %{$list};
	open O, ">$bed.bed";
	for $tp(keys %list){
		for $var(keys %{$list{$tp}}){
			($chr, $stt, $end) = split /\t/, $var;
			print O "$chr\t$stt\t$end\n";
		}        
	}
	close O || die;
	system "sort -k1,1 -k2,2n $bed.bed > $bed.sorted.bed";
	return()
}


#
sub overlap{ # overlap("query.bed", "feature.bed"); need to be sorted
	my ($query, $feature) = @_;
	my %over;
	open I, "bedtools intersect -wa -wb -a $query -b $feature -sorted |";
	while(<I>){
		chomp;
		($chr1, $stt1, $end1, $chr2, $stt2, $end2) = split /\t/;
		@tmp = sort{$a <=> $b}($stt1, $end1, $stt2, $end2);
		$over{"$chr1\t$stt1\t$end1"}{"$chr2\t$stt2\t$end2"} = ($tmp[2] - $tmp[1])/($end2-$stt2);
	}
	close I;
	return (\%over);
} 

# annotate SV with RE, targets, and impact group of RE
sub rg_anno{
	my ($sv2re, $re_anno, $rg_file) = @_;
	my %sv2re = %{$sv2re};
	my %loc2re; my %sv2target;
	open F, $re_anno;
	while(<F>){
		($chr, $stt, $end, $re) = split /\t/;
		$loc2re{"$chr\t$stt\t$end"}=$re;
	}
	close F || die;
	open I, $rg_file || die;
	while(<I>){
		chomp; ($re, $gn) = split /\t/;
		$re2gn{$re}{$gn}=1;
	}
	close I || die;
	for $var(keys %sv2re){
		my @set1; my @set2; my @set3;
		for $re_loc(keys %{$sv2re{$var}}){
			$ratio = $sv2re{$var}{$re_loc};
			$re = $loc2re{$re_loc}; $loc = join ":", (split "\t", $re_loc);
			if($ratio < 0.5){
				for $gn(keys %{$re2gn{$re}}){
					push @set1, $loc."--".$gn;
				}
			}elsif($ratio >= 1){
				for $gn(keys %{$re2gn{$re}}){
					push @set3, $loc."--".$gn;
				}			
			}else{
				for $gn(keys %{$re2gn{$re}}){
					push @set2, $loc."--".$gn;
				}
			}
		}
my $j3 = join ";", @set3; $j3 = "-" if $j3 eq "";
my $j2 = join ";", @set2; $j2 = "-" if $j2 eq "";
my $j1 = join ";", @set1; $j1 = "-" if $j1 eq "";
		$sv2target{$var}="$j3|$j2|$j1";
}
	return (\%sv2target);
}

# annotate FI for RE based on 1000Genome SV calls
sub read_fi{
	my ($file) = @_;
	my %fi;
	open I, $file || die; 
	$tmp=<I>;
	while(<I>){
		chomp;
		($tp, $chr, $stt, $end, $fi3, $fi2, $fi1)=split /\t/;
		$re_loc = "$chr:$stt:$end";
		$fi{$tp}{$re_loc}{"set3"} = $fi3;
		$fi{$tp}{$re_loc}{"set2"} = $fi3;
		$fi{$tp}{$re_loc}{"set1"} = $fi3;
	}
	close I || die;
	return (\%fi);
}

sub anno_fi{
	my ($sv2target, $to_use, $fi) = @_;
	%fi = %{$fi}; %sv2target = %{$sv2target}; %to_use = %{$to_use};
	my %sv2fi;
	for $var(keys %sv2target){
		next if(!exists $to_use{$var});
		for $tp0(keys %{$to_use{$var}}){
			if($tp0 =~ /DUP/){
				$tp = "DUP";
			}elsif($tp0 =~ /DEL/){
				$tp = "DEL";
			}elsif($tp0 =~ /CNV/){
				$tp = "CNV";
			}elsif($tp0 =~ /INV/){
				$tp = "INV";
			}elsif($tp0 =~ /INS/){
				$tp = "INS";
			}else{
				$tp = "all"; # if cannot match the SV type, just use the FI by all SVs impacting the RE
			}
			undef @ts3; undef @ts2; undef @ts1; my @ts = ("-","-","-");
			($ts[0], $ts[1], $ts[2]) = split /\|/, $sv2target{$var};
			undef @set3; undef @set2; undef @set1;
			if($ts[0] ne "-"){
				@ts3 = split ";", $ts[0];
				for $tmp(@ts3){
					undef $re_loc;
					($re_loc) = split "--", $tmp;
					if(exists $fi{$tp}{$re_loc}{"set3"}){
						push @set3, $fi{$tp}{$re_loc}{"set3"}."($tp)";
# if no FI for the specific SV type, try using the FI by all SVs
					}elsif(exists $fi{"all"}{$re_loc}{"set3"}){
						push @set3, $fi{"all"}{$re_loc}{"set3"}."(all)";
					}else{
						push @set3, "";
					}
				}
			}
			if($ts[1] ne "-"){
				@ts2 = split ";", $ts[1];
				for $tmp(@ts2){
					undef $re_loc;
					($re_loc) = split "--", $tmp;
					if(exists $fi{$tp}{$re_loc}{"set2"}){
						push @set2, $fi{$tp}{$re_loc}{"set2"}."($tp)";
					}elsif(exists $fi{"all"}{$re_loc}{"set2"}){
						push @set2, $fi{"all"}{$re_loc}{"set2"}."(all)";
					}else{
						push @set2, "";
					}
				}
			}
			if($ts[2] ne "-"){
				@ts1 = split ";", $ts[2];
				for $tmp(@ts1){
					undef $re_loc;
					($re_loc) = split "--", $tmp;
					if(exists $fi{$tp}{$re_loc}{"set1"}){
						push @set1, $fi{$tp}{$re_loc}{"set1"}."($tp)";
					}elsif(exists $fi{"all"}{$re_loc}{"set1"}){
						push @set1, $fi{"all"}{$re_loc}{"set1"}."(all)";
					}else{
						push @set1, "";
					}
				}
			}
my $j3 = join ";", @set3; $j3 = "-" if $j3 eq "";
my $j2 = join ";", @set2; $j2 = "-" if $j2 eq "";
my $j1 = join ";", @set1; $j1 = "-" if $j1 eq "";
			$sv2fi{$var}{$tp0} = "$j3|$j2|$j1"; 
		}
	}
	return (\%sv2fi);
}

#
sub write_svint_anno{
	my ($out_file, $sv_file, $lines_a, $to_use, $order, $gs, $fs, $g2f_names) = @_;
	@order = @{$order}; %gs = %{$gs}; %fs = %{$fs};  %g2f_names = %{$g2f_names}; %lines_a = %{$lines_a}; %to_use = %{$to_use};
	undef @gs; undef @fs; undef @gfs_names;
	for $i(0 .. $#order){
		$tar = $order[$i];
		push @gs, $gs{$tar}; push @fs, $fs{$tar}; 
		push @gfs_names, $g2f_names{$tar}."(100% IMPACT)"; 
		push @gfs_names, $g2f_names{$tar}."(50-100% IMPACT)";
		push @gfs_names, $g2f_names{$tar}."(0-50% IMPACT)";
		push @gfs_names, "FI_for_".$g2f_names{$tar}."(100% IMPACT)";
		push @gfs_names, "FI_for_".$g2f_names{$tar}."(50-100% IMPACT)";
		push @gfs_names, "FI_for_".$g2f_names{$tar}."(0-50% IMPACT)";
	}
	open O, ">$out_file" || die;
	print O "##SVint version1.0 annotated targets and FI (Frequency of Impact for Regulatory elements)  for intergenic SVs\n";
	open I, $sv_file || die;
	while(<I>){
		chomp; $tmp = $_; last if $tmp !~ /^#/;
		if($tmp =~ /^##/){
			print O $tmp."\n";
		}elsif($tmp =~ /^#CHR/){
			print O "##SVINT Format: ".(join "|", @gfs_names)."\n";
			print O "$tmp\tSVINT\n";
		}
	}
	close I || die;
	for $id(sort{$a<=>$b}(keys %lines_a)){
		for $line(keys %{$lines_a{$id}}){
			undef @svint; undef @tt; undef $tp0;
			$var = $lines_a{$id}{$line};
			if(scalar keys %{$to_use{$var}}<1){
				print O "$line\t.\n"; next;
			}
			@tt = split /\t/, $line; $alt = $tt[4];
			if($alt =~ /<(.*)>/){
				$tp0 = $1;
			}else{
				undef @in;
				@in = split ";", $tt[7];
				for $in(@in){
					if($in =~ /SVTYPE2=(.*)/){
						$tp0 = $1;
					}elsif($in =~ /SVTYPE=(.*)/){
						$tp0 = $1;
					}
				}
			}
			for $i(0 .. $#gs){
				%tmp1 = %{$gs[$i]}; %tmp2 = %{$fs[$i]};
#print "$tp0\t$var\t$tmp1{$var}\t$tmp2{$var}{$tp0}\n";
				push @svint, "-|-|-"; $svint[$#svint] = $tmp1{$var} if exists $tmp1{$var};
				push @svint, "-|-|-"; $svint[$#svint] = $tmp2{$var}{$tp0} if exists $tmp2{$var}{$tp0};
			}
			print O "$line\t".(join "|", @svint)."\n";
		}
	}
	close O || die;
}



1;
