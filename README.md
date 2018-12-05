# voicerecogn-cw

## Training 
* train.sh 

## Testing 
* HVite -T 1 -S lists/list1.txt -d hmms/ -w lib/NET -l results lib/dict lib/words3

## Results
* HResults -p -e "???" Silent -L labels/test lib/words3 results/*.rec