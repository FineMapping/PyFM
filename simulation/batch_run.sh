ld_dir=/media/zhaoyang-new/workspace/westonelison_finemap_simulator/CSE_284_Finemapping/simulations/1KGP_Regions/

file=$1
output_dir=$2

echo $file
region=$(echo "file" | cut -d. -f1)
seed=$(echo "file" | awk -F_ '{print $NF}')
file_name=$(basename "$file")
region="${file_name%%.vcf.*}"        # Removes _CV_* from the end
seed="${file_name##*gz_}"        # Removes * from the beginning
seed="${seed%%_CV_*}"
middle="${file_name#*_CV_}"         # Removes shortest match of *_CV_ from the beginning
cv="${middle%%.*}"             # Removes shortest match of .* from the end

iteration=300

for k in 3; do
    for r in {1..5}; do
      python src/fm.py \
        -z "${file}" \
        -r ${ld_dir}/"${region}".vcf.gz.r.maf.order.ld \
        -o "${output_dir}/i${iteration}_K${k}_r${r}_cv${cv}_region${region}_seed${seed}" \
        --configs-method SSSConfigurations \
        --SSS-iterations ${iteration} \
        -n 20000 -c ${k} -t 0 -e 0.1 -a 1.6 -p 1
    done
done




