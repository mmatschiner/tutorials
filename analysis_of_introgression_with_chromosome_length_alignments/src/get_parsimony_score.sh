# m_matschiner Mon Sep 17 18:49:14 CEST 2018

# Get the command-line arguments.
nex=${1}

# Prepare a temporary paup command block.
echo "begin paup;" > tmp_run_paup.nex
echo "  execute ${nex};" >> tmp_run_paup.nex
echo "  set autoclose=yes warntree=no warnreset=no;" >> tmp_run_paup.nex
echo "  log start file=tmp_paup_log.txt replace;" >> tmp_run_paup.nex 
echo "  set criterion=parsimony;" >> tmp_run_paup.nex
echo "  hsearch addseq=simple nreps=1000 swap=tbr hold=1;" >> tmp_run_paup.nex
echo "  pscores / ci=no ri=no rc=no scoreFile=tmp_paup_out.txt;" >> tmp_run_paup.nex
echo "end;" >> tmp_run_paup.nex

# Execute the paup command block.
paup -n tmp_run_paup.nex &> /dev/null

# Get the parsimony score from the paup pscores output file.
parsimony_score=`cat tmp_paup_out.txt | head -n 2 | tail -n 1 | cut -f 2`

echo ${parsimony_score}

# Clean up.
rm tmp_run_paup.nex
rm tmp_paup_log.txt
rm tmp_paup_out.txt