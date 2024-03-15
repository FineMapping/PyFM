for k in 5; do
  for iteration in 1000; do
    for r in {1..100}; do
      python src/fm.py \
        -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
        -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.ld \
        -o test_output/i"${iteration}"_K"${k}"_r"${r}" \
        --configs-method SSSConfigurations \
        --SSS-iterations ${iteration} \
        -n 471 -c ${k} -t 0 -e 0.1 -a 1.6 -p 0.95
    done
  done
done


#o=../pyfm_results_c3_sss1000_${i}
#mkdir $o
#python src/fm.py \
	#-z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z  \
	#-r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD  \
	#-o $o  \
	#--configs-method SSSConfigurations  \
	#--SSS-iterations 1000  \
	#-n 471 -c 3 -t 0 -e 0.1 -a 1.6 --Random-Seed ${i} > ${o}_pyfm.log



