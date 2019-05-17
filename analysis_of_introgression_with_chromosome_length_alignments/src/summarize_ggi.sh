# m_matschiner Fri Nov 30 13:40:26 CET 2018

# Get the command-line arguments.
alignment_dir=${1}
output_table=${2}

# Specify the output table.
if [ ! -f ${output_table} ]
then
	echo -e "block_id\tbest\tdelta_lik\tlik_t1\tlik_t2\tlik_t3" >> ${output_table}
	for i in ${alignment_dir}/*.nex
	do
		block_id=`basename ${i%.nex}`

		# Test if all three iqtree info files exist.
		t1_file=${i%.nex}.t1.log
		t2_file=${i%.nex}.t2.log
		t3_file=${i%.nex}.t3.log
		if [[ -f ${t1_file} && -f ${t2_file} && -f ${t3_file} ]]
		then

			# Make sure that all files contain likelihoods.
			for t_file in ${t1_file} ${t2_file} ${t3_file}
			do
				lik_included=`cat ${t_file} | grep "Log-likelihood of the tree" | wc -l`
				if [[ ${lik_included} == 0 ]]
				then
					echo "ERROR: No likelihood found in file ${t_file}!"
					exit 1
				fi
			done

			# Get the three likelihoods.
			t1_lik=`cat ${t1_file} | grep "Log-likelihood of the tree" | cut -d ":" -f 2 | cut -d " " -f 2`
			t2_lik=`cat ${t2_file} | grep "Log-likelihood of the tree" | cut -d ":" -f 2 | cut -d " " -f 2`
			t3_lik=`cat ${t3_file} | grep "Log-likelihood of the tree" | cut -d ":" -f 2 | cut -d " " -f 2`

			# Determine the best and second-best likelihood, and their difference.
			best_lik=`echo -e "${t1_lik}\n${t2_lik}\n${t3_lik}" | sort -n -r | head -n 1`
			second_best_lik=`echo -e "${t1_lik}\n${t2_lik}\n${t3_lik}" | sort -n -r | head -n 2 | tail -n 1`
			delta_lik=`echo "${best_lik} - ${second_best_lik}" | bc`
			if (( 1 == `echo "${t1_lik} > ${t2_lik}" | bc` )) && (( 1 == `echo "${t1_lik} > ${t3_lik}" | bc` ))
			then
				echo -e "${block_id}\tt1\t${delta_lik}\t${t1_lik}\t${t2_lik}\t${t3_lik}" >> ${output_table}
			elif (( 1 == `echo "${t2_lik} > ${t1_lik}" | bc` )) && (( 1 == `echo "${t2_lik} > ${t3_lik}" | bc` ))
			then
				echo -e "${block_id}\tt2\t${delta_lik}\t${t1_lik}\t${t2_lik}\t${t3_lik}" >> ${output_table}
			elif (( 1 == `echo "${t3_lik} > ${t1_lik}" | bc` )) && (( 1 == `echo "${t3_lik} > ${t2_lik}" | bc` ))
			then
				echo -e "${block_id}\tt3\t${delta_lik}\t${t1_lik}\t${t2_lik}\t${t3_lik}" >> ${output_table}
			fi
		fi
	done
fi
