#!/bin/bash

#### INPUT PARAMETERS ####
help(){
echo "==== VARIABLE ARGUMENTS ===="
echo "-p : protein name"
echo "-f : protein file"
echo "-s : study name"
echo "-r : snp file"
echo "-t : file type (1=vcf (default); 0=txt,tsv,etc.)"
echo "-D : vcf lines to skip (default 34)"
}

while getopts p:f:s:r:t:D:h flag
do
    case "${flag}" in
        p) prot=${OPTARG};;
	f) file=${OPTARG};;
	s) study=${OPTARG};;
	r) snp=${OPTARG};;
	t) type=${OPTARG};;
	D) vcfskip=${OPTARG};;
        h) help; exit;;
    esac
done

declare prot_name=$prot
declare prot_file=$file
declare snp_file=$snp
declare study_name=$study
if [ -z "$type" ]; then type=1; fi;
if [ -z "$vcfskip" ]; then vcfskip=34; fi;
declare file_type=$type


if [ `echo $prot_file | grep gz` ]
then
	com_grep="zgrep"
	com_cat="zcat"
else
	com_grep="grep"
	com_cat="cat"
fi;


if [ $file_type == 1 ]
then
        echo "File is in vcf format"
	vcfskip=$((vcfskip+1))
	$com_cat $prot_file | head -n$vcfskip | tail -n1 >> ./temp_files/$study_name/$prot_name.header.txt
else
        echo "File is not in vcf format"
	$com_cat $prot_file | head -n1 >> ./temp_files/$study_name/$prot_name.header.txt
fi;

#$com_cat $prot_file | head -n1 >> ./temp_files/$study_name/header.txt
$com_grep -f $snp_file $prot_file >> ./temp_files/$study_name/$prot_name.reverse.noheader.txt
cat ./temp_files/$study_name/$prot_name.header.txt ./temp_files/$study_name/$prot_name.reverse.noheader.txt >> ./temp_files/$study_name/$prot_name.reverse.txt
rm ./temp_files/$study_name/$prot_name.header.txt ./temp_files/$study_name/$prot_name.reverse.noheader.txt

