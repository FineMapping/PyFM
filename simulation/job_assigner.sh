data_dir=/media/zhaoyang-new/workspace/PyFM/example/cv2_data/
output_dir=simulated_data_run_cv4_k3

for f in "${data_dir}"/*
do
  bash batch_run.sh "${f}" ${output_dir} &
done