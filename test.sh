iteration=100
k=5

for r in {1..3}; do
  python src/fm.py \
    -z example/finemap_examples/finemap_data.standardized.z \
    -r example/finemap_examples/finemap_data.ld \
    -o sampling_run/I"${iteration}"_K"${k}"_run"${r}" \
    --configs-method SSSConfigurations \
    --SSS-iterations ${iteration} \
    -n 5363 -c ${k} -t 0 -e 0.1 -a 1.6
done