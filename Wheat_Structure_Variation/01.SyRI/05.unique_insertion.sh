#get all assembly all chromosome lists
ls */syri/*/*.out|grep -v 3504 >list.syri.txt
#get 3504 all chromosome syri out
ls */syri/*/*.out|grep 3504|while read line
  do
  cat $line >>3504_CS21syri.out
  done
#run script unique_ins.py 
python unique_ins.py -F 3504_CS21syri.out -L list.syri.txt -O 3504_CS21syri.unique_insertion.out
