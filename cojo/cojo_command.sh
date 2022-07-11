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
echo "-l : p-value format (1=log10P (default); 0=P)"
echo "-n : study name (default 'deCODE')"
echo "-r : retry analysis if no snp was found in case allele cols are inverted"
echo "-i : does reference file contains SNP rsids as reference instead of chr:pos:a1:a2 ? (1=yes; 0=no) (default 0)"
echo "-A : SNP column number (default 3)"
echo "-B : A1 column number (default 4)"
echo "-C : A2 column number (default 5)"
echo "-D : MAF column number (if file is a vcf, this column must correspond to its order in the format column + 3 (if format column goes like AF:ES:SE:LP:SS, then AF column number should be 4 (1+3)) (default 6)"
echo "-E : Beta column number (if file is a vcf, this column must correspond to its order in the format column + 3 (if format column goes like AF:ES:SE:LP:SS, then ES column number should be 5 (2+3)) (default 7)"
echo "-F : SE column number (if file is a vcf, this column must correspond to its order in the format column + 3 (if format column goes like AF:ES:SE:LP:SS, then SE column number should be 6 (3+3)) (default 8)"
echo "-G : P-value column number (if file is a vcf, this column must correspond to its order in the format column + 3 (if format column goes like AF:ES:SE:LP:SS, then LP column number should be 7 (4+3)) (def9ult 9)"
echo "-H : N column number (if file is a vcf, this column must correspond to its order in the format column + 3 (if format column goes like AF:ES:SE:LP:SS, then N column number should be 8 (1+3)) (default 10)"
echo "-I : If file is a vcf, how long is the header (number of lines excluding column names) (default 34)"
echo "-J : If file is a vcf, format column number (default 12)"
echo "-K : Path to LD reference files in Plink BED format WITHOUT chromosome number, bed extension or '.' at the end (default './ldref/chr1_22.b37/all_phase3')"
echo "-L : maf_threshold (default 0.01)"
echo "-M : chromosome column number (default 1)"
echo "-N : position column number (default 2)"
}

while getopts p:f:c:s:e:t:l:n:r:i:A:B:C:D:E:F:G:H:I:J:K:L:M:N:h flag
do
    case "${flag}" in
        p) prot=${OPTARG};;
	f) file=${OPTARG};;
	c) chr=${OPTARG};;
	s) start=${OPTARG};;
	e) end=${OPTARG};;
	t) type=${OPTARG};;
	l) logp=${OPTARG};;
	n) study=${OPTARG};;
	r) retry=${OPTARG};;
	i) refrsids=${OPTARG};;
	A) snpcol=${OPTARG};;
	B) a1col=${OPTARG};;
	C) a2col=${OPTARG};;
	D) mafcol=${OPTARG};;
	E) betacol=${OPTARG};;
	F) secol=${OPTARG};;
	G) pcol=${OPTARG};;
	H) ncol=${OPTARG};;
	I) vcfskip=${OPTARG};;
	J) formcol=${OPTARG};;
	K) ldref=${OPTARG};;
	L) mafthresh=${OPTARG:-0.01};;
	M) chrcol=${OPTARG:-1};;
	N) poscol=${OPTARG:-2};;
        h) help; exit;;
    esac
done

if [ -z "$type" ]; then type=1; fi;
if [ -z "$logp" ]; then logp=1; fi;
if [ -z "$study" ]; then study="deCODE"; fi;
if [ -z "$retry" ]; then retry=0; fi;
if [ -z "$refrsids" ]; then refrsids=0; fi;
if [ -z "$snpcol" ]; then snpcol=3; fi;
if [ -z "$a1col" ]; then a1col=4; fi;
if [ -z "$a2col" ]; then a2col=5; fi;
if [ -z "$mafcol" ]; then mafcol=6; fi;
if [ -z "$betacol" ]; then betacol=7; fi;
if [ -z "$secol" ]; then secol=8; fi;
if [ -z "$pcol" ]; then pcol=9; fi;
if [ -z "$ncol" ]; then ncol=10; fi;
if [ -z "$vcfskip" ]; then vcfskip=34; fi;
if [ -z "$formcol" ]; then formcol=12; fi;
if [ -z "$ldref" ]; then ldref="./cojo/ldref/chr1_22.b37/all_phase3"; fi;
if [ -z "$mafthresh" ]; then mafthresh=0.01; fi;
if [ -z "$chrcol" ]; then chrcol=1; fi;
if [ -z "$poscol" ]; then poscol=2; fi;

declare prot_name=$prot
declare prot_file=$file
## start and end include 500kb window
# declare -i chr=$chr
declare -i start=$start
declare -i end=$end
declare -l file_type=$type
declare -i log_type=$logp
declare -i retry=$retry
declare -i ref_file_rsids=$refrsids
declare -i snp_col=$snpcol
declare -i a1_col=$a1col
declare -i a2_col=$a2col
declare -i maf_col=$mafcol
declare -i beta_col=$betacol
declare -i se_col=$secol
declare -i p_col=$pcol
declare -i n_col=$ncol
declare -i format_col=$formcol
declare -i vcfskip=$vcfskip
declare -i chr_col=$chrcol
declare -i pos_col=$poscol

if [ `echo $prot_file | grep gz` ]
then com_cat="zcat"
else com_cat="cat"
fi;

cd ./cojo/
if [ $file_type == 1 ]
then
	echo "File is in vcf format"
	if [ $ref_file_rsids == 0 ]
	then
		echo "SNPs reference file contains SNPs as chr:pos:a1:a2"
		$com_cat $prot_file | awk -v chr="${chr}" -v start=$start -v end=$end -v skip=$vcfskip -v chrc=$chr_col -v pos=$pos_col '(NR>skip && $chrc==chr && $pos>=start && $pos<=end)' | awk -v chrc=$chr_col -v pos=$pos_col -v a1=$a1_col -v a2=$a2_col -v form=$format_col '{print $chrc "@" $pos "@" $a1 "@" $a2 "\t" $a1 "\t" $a2 "\t" $form}' | sed 's/:/\t/g' | awk -v maf=$maf_col -v b=$beta_col -v se=$se_col -v p=$p_col -v n=$n_col '{print $1 "\t" $2 "\t" $3 "\t" $maf "\t" $b "\t" $se "\t" $p "\t" $n}' | sed 's/@/:/g' >> ./temp_files/${study}/${prot_name}.txt
	else
		echo "SNPs reference file contains SNPs as rsids"
		$com_cat $prot_file | awk -v chr="${chr}" -v start=$start -v end=$end -v skip=$vcfskip -v chrc=$chr_col -v pos=$pos_col '(NR>skip && $chrc==chr && $pos>=start && $pos<=end)' | awk -v snp=$snp_col -v a1=$a1_col -v a2=$a2_col -v form=$format_col '{print $snp "\t" $a1 "\t" $a2 "\t" $form}' | sed 's/:/\t/g' | awk -v maf=$maf_col -v b=$beta_col -v se=$se_col -v p=$p_col -v n=$n_col '{print $1 "\t" $2 "\t" $3 "\t" $maf "\t" $b "\t" $se "\t" $p "\t" $n}' >> ./temp_files/${study}/${prot_name}.txt
	fi;
	if [ $log_type == 1 ]
	then
		awk -v p=$p_col '{print v=10^-$p}' ./temp_files/${study}/${prot_name}.txt >> ./temp_files/${study}/${prot_name}.pval.txt
		paste ./temp_files/${study}/${prot_name}.txt ./temp_files/${study}/${prot_name}.pval.txt | cut -f1-6,8-9 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $8 "\t" $7}' >> ./temp_files/${study}/${prot_name}.clean.txt
		rm ./temp_files/${study}/${prot_name}.pval.txt
	else
		awk -v snp=$snp_col -v a1=$a1_col -v a2=$a2_col -v maf=$maf_col -v b=$beta_col -v se=$se_col -v n=$n_col '{print $snp "\t" $a1 "\t" $a2 "\t" $maf "\t" $b "\t" $se "\t" $8 "\t" $n}' >> ./temp_files/${study}/${prot_name}.clean.txt
	fi;
	rm ./temp_files/${study}/${prot_name}.txt
	cat header.vcf.txt ./temp_files/${study}/${prot_name}.clean.txt >> ./input/${study}/${prot_name}.txt
	rm ./temp_files/${study}/${prot_name}.clean.txt
else
	echo "File is not in vcf format"
	if [ $ref_file_rsids == 0 ]
	then
		echo "SNPs reference file contains SNPs as chr:pos:a1:a2"
		$com_cat $prot_file | awk -v chr="${chr}" -v start=$start -v end=$end -v chrc=$chr_col -v pos=$pos_col '($chrc==chr && $pos>=start && $pos<=end)' | awk -v chrc=$chr_col -v pos=$pos_col -v a1=$a1_col -v a2=$a2_col -v maf=$maf_col -v b=$beta_col -v se=$se_col -v p=$p_col -v n=$n_col '{print $chrc ":" $pos ":" $a1 ":" $a2 "\t" $a1 "\t" $a2 "\t" $maf "\t" $b "\t" $se "\t" $p "\t" $n}' >> ./temp_files/${study}/${prot_name}.clean.txt
	else
		echo "SNPs reference file contains SNPs as rsids"
		$com_cat $prot_file | awk -v chr="${chr}" -v start=$start -v end=$end -v chrc=$chr_col -v pos=$pos_col '($chrc==chr && $pos>=start && $pos<=end)' | awk -v snp=$snp_col -v a1=$a1_col -v a2=$a2_col -v maf=$maf_col -v b=$beta_col -v se=$se_col -v p=$p_col -v n=$n_col '{print $snp "\t" $a1 "\t" $a2 "\t" $maf "\t" $b "\t" $se "\t" $p "\t" $n}' >> ./temp_files/${study}/${prot_name}.clean.txt
	fi;
        cat header.vcf.txt ./temp_files/${study}/${prot_name}.clean.txt >> ./input/${study}/${prot_name}.txt
        rm ./temp_files/${study}/${prot_name}.clean.txt
fi;

#### COJO ANALYSIS ####
declare ld_ref=$ldref
maf_threshold=$mafthresh

if [ `echo $chr | grep "chr"` ]
then
        chr=`echo $chr | sed 's/chr//g'`
        sed 's/chr//g' ./input/${study}/${prot_name}.txt >> ./input/${study}/${prot_name}.nochr.txt
	rm ./input/${study}/${prot_name}.txt
	mv ./input/${study}/${prot_name}.nochr.txt ./input/${study}/${prot_name}.txt
fi;

if [ $retry == 1 ]
then
	echo "Retrying CoJo analysis by inverting alleles in SNP name from GWAS file..."
	awk '(NR>1)' ./input/${study}/${prot_name}.txt | sed 's/:/\t/g' | cut -f1-4 | awk '{print $1 ":" $2 ":" $4 ":" $3}' | grep -v 'NA' >> ./input/${study}/${prot_name}.newid.txt
	echo -e "SNP" | cat - ./input/${study}/${prot_name}.newid.txt >> ./input/${study}/${prot_name}.newid.clean.txt
	awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' ./input/${study}/${prot_name}.newid.clean.txt ./input/${study}/${prot_name}.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' >> ./input/${study}/${prot_name}.final.txt
	rm ./input/${study}/${prot_name}.txt
	rm ./input/${study}/${prot_name}.newid.txt
	rm ./input/${study}/${prot_name}.newid.clean.txt
	mv ./input/${study}/${prot_name}.final.txt ./input/${study}/${prot_name}.txt
fi;

cd ..

# You can also add the option 'cojo-gc' to adjust p-values by the genomic control method (inflation lambda)
gcta64 --bfile ${ld_ref}.${chr} --chr $chr --maf $maf_threshold --cojo-file ./cojo/input/${study}/${prot_name}.txt --cojo-slct --out ./cojo/output/${study}/${prot_name} --cojo-wind 1000 --cojo-collinear 0.9 --diff-freq 0.2 --cojo-p 5e-8
rm ./cojo/input/${study}/${prot_name}.txt
