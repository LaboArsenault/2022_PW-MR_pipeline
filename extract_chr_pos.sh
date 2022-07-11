#!/bin/bash

#### INPUT PARAMETERS ####
help(){
echo "==== VARIABLE ARGUMENTS ===="
echo "-p : protein name"
echo "-f : protein file"
echo "-c : chromosome"
echo "-s : start of gene (including window)"
echo "-e : end of gene (including window)"
echo "-t : file type (1=vcf (default); 0=txt,tsv,etc.)"
echo "-n : study name (default deCODE)"
echo "-A : SNP column number (default 3)"
echo "-B : chromosome column number (default 1)"
echo "-C : position column number (default 2)"
echo "-D : vcf lines to skip (default 34)"
}

while getopts p:f:c:s:e:t:n:A:B:C:D:h flag
do
    case "${flag}" in
        p) prot=${OPTARG};;
	f) file=${OPTARG};;
	c) chr=${OPTARG};;
	s) start=${OPTARG};;
	e) end=${OPTARG};;
	t) type=${OPTARG:-1};;
	n) study=${OPTARG:-1};;
	A) snpcol=${OPTARG:-3};;
	B) chrcol=${OPTARG:-4};;
	C) poscol=${OPTARG:-5};;
	D) vcfskip=${OPTARG:-5};;
        h) help; exit;;
    esac
done

if [ -z "$type" ]; then type=1; fi;
if [ -z "$study" ]; then study="deCODE"; fi;
if [ -z "$snpcol" ]; then snpcol=3; fi;
if [ -z "$chrcol" ]; then chrcol=1; fi;
if [ -z "$poscol" ]; then poscol=2; fi;
if [ -z "$vcfskip" ]; then vcfskip=34; fi;

declare prot_name=$prot
declare prot_file=$file
## start and end include 500kb window
declare -i start=$start
declare -i end=$end
declare -l file_type=$type
declare -i snp_col=$snpcol
declare -i chr_col=$chrcol
declare -i pos_col=$poscol

if [ `echo $prot_file | grep gz` ]
then
	com_cat="zcat"
	com_grep="zgrep"
else
	com_cat="cat"
	com_grep="grep"
fi;

if [ $file_type == 1 ]
then
	echo "File is in vcf format"
	$com_cat $prot_file | awk -v chr=$chr -v start=$start -v end=$end -v skip=$vcfskip -v chrc=$chr_col -v pos=$pos_col -v snp=$snp_col '(NR>skip && $chrc==chr && $pos>=start && $pos<=end) {print $snp}' >> ./temp_files/${study}/${prot_name}.clean.1.txt
	$com_cat $prot_file | awk -v skip=$vcfskip '(NR==skip+1)' >> ./temp_files/${study}/${prot_name}.header.txt
else
	echo "File is not in vcf format"
	$com_cat $prot_file | awk -v chr=$chr -v start=$start -v end=$end -v chrc=$chr_col -v pos=$pos_col -v snp=$snp_col '($chrc==chr && $pos>=start && $pos<=end) {print $snp}' >> ./temp_files/${study}/${prot_name}.clean.1.txt
	$com_cat $prot_file | awk '(NR==1)' >> ./temp_files/${study}/${prot_name}.header.txt
fi;

grep -v "NA" ./temp_files/${study}/${prot_name}.clean.1.txt >> ./temp_files/${study}/${prot_name}.clean.2.txt
grep -v "," ./temp_files/${study}/${prot_name}.clean.2.txt >> ./temp_files/${study}/${prot_name}.txt
$com_grep -f ./temp_files/${study}/${prot_name}.txt $prot_file >> ./temp_files/${study}/${prot_name}.region.nohead.txt
cat ./temp_files/${study}/${prot_name}.header.txt ./temp_files/${study}/${prot_name}.region.nohead.txt >> ./temp_files/${study}/${prot_name}.region.txt
rm ./temp_files/${study}/${prot_name}.txt ./temp_files/${study}/${prot_name}.clean.1.txt ./temp_files/${study}/${prot_name}.clean.2.txt ./temp_files/${study}/${prot_name}.header.txt ./temp_files/${study}/${prot_name}.region.nohead.txt
