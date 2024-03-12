../caviarbf/caviarbf \
  -z example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z \
  -r example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD \
  -o /media/zhaoyang-new/Photo_Archive/k5/k5.bf \
  -n 471 -c 5 -t 0 -e 0.1 -a 1.6

../caviarbf/model_search \
  -i /media/zhaoyang-new/Photo_Archive/k5/k5.bf \
  -o /media/zhaoyang-new/Photo_Archive/k5/ \
  -s -m 237 -p 0 > /media/zhaoyang-new/workspace/PyFM/k5/log.txt

../caviarbf/model_search \
  -i /media/zhaoyang-new/workspace/PyFM/k4/k4.bf \
  -o /media/zhaoyang-new/workspace/PyFM/k4/k4 \
  -s -m 237 -p 0 > /media/zhaoyang-new/workspace/PyFM/k4/log.txt