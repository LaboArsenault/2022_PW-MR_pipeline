#!/usr/bin/env bash

# ================== INFO =====================	#
#						#
# This is a bash script to loop on deCODE files	#
# splitting them by batches of 500 proteins.	#
# This script can be used with '.vcf' files,	#
# but its main purpose is to be used when	#
# '.txt.gz' files are used, since the 'PWAS.R'	#
# script will run quite faster on '.vcf' files	#
# (~1 day) than on '.txt.gz' files (~4-5 days).	#
#						#
# This script automatically uses 'NOHUP' to run #
# the analyses.					#
#						#
#		*	*	*		#
# MAKE SURE YOU FILL IN THE 'SPECIFY OPTIONS'	#
# PARTS WITH WHAT IS SPECIFIED ABOVE THEM.	#
#		*	*	*		#
#						#
# ============================================= #

# ************ 'SPECIFIY OPTIONS' ************* #
# Specify the number of proteins total that are present in the 'genes_ref_file' (for deCODE : 4549)
n_proteins=''
# Specify the number of batches (for deCODE : 10 is appropriate to run the analysis)
n_batches=''
# ********************************************* #

for i in {1..${n_batches}}
do
	declare -i start=$(((($i-1)*500)+1))
	declare -i end=$(($i*500))
	if [ $i == ${n_batches} ]
	then
		declare -i end=${n_proteins}
	fi;

	# ************ 'SPECIFIY OPTIONS' ************* #
	# Here, fill in the parts between ' ', uncomment the line and run the bash script.
	# nohup Rscript ./PWAS.R '[all arguments to supply]' --from $start --to $end >> ./nohup/'[name of your nohup outfile]'.from${start}.to${end}.out &
	# ********************************************* #
done;
