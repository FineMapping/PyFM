ld_dir=/media/zhaoyang-new/workspace/westonelison_finemap_simulator/CSE_284_Finemapping/simulations/1KGP_Regions/
region=1KGP_hg19_AIRE_1MB

../caviarbf/caviarbf \
  -z example/simulated_data/1KGP_hg19_AIRE_1MB.vcf.gz_104851018_CV_5.caviarbf.z \
  -r ${ld_dir}/"${region}".vcf.gz.r.maf.order.ld \
  -o /media/zhaoyang-new/Photo_Archive/k5/k5.bf \
  -n 20000 -c 3 -t 0 -e 0.1 -a 1.6

../caviarbf/model_search \
  -i /media/zhaoyang-new/Photo_Archive/k5/k5.bf \
  -o /media/zhaoyang-new/Photo_Archive/k5/ \
  -s -m 3192 -p 0 > /media/zhaoyang-new/Photo_Archive/k5/log.txt
#
#../caviarbf/model_search \
#  -i /media/zhaoyang-new/workspace/PyFM/CV4/CV4.bf \
#  -o /media/zhaoyang-new/workspace/PyFM/CV4/CV4 \
#  -s -m 237 -p 0 > /media/zhaoyang-new/workspace/PyFM/CV4/log.txt