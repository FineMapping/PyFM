for k in 7; do
  for iteration in 700 ; do
    for r in {1..30}; do
      python src/fm.py \
        -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
        -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD \
        -o sampling_run/i"${iteration}"_K"${k}"_r"${r}" \
        --configs-method SSSConfigurations \
        --SSS-iterations ${iteration} \
        -n 5363 -c ${k} -t 0 -e 0.1 -a 1.6 -p 0.95
    done
  done
done



