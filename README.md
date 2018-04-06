# SVint
a light-weight tool for annotating structure variants located outside the coding genome



0. dependants: perl; bedtools;

1. download "SVint.zip", and unzip the whole thing

2. go to your-directory/SVint/, and use this as your working directory
>your-directory is where you unzip your SVint.zip

3.if you want to annotate an SV VCF from one single sample, please first put your data (presumably an SV VCF file) into your working directory (your-directory/SVint/), and then run this line: 
perl src/SVint_anno_single.pl your-data.vcf
>the results will be written in your working directory

