dimusign=true
accW=true
effW=true
isMC=false
HFSel=0

root -l -b <<EOF
.L makeV2Hist.C++
.q
EOF
rm *.so *.d *.pcm

echo make hist... 

for pt in '0,3' '0,50' '3,6' '6,50' 
do
  for cBin in '20,180' '20,100' '100,180' '0,20' 
  do
    for muPt in '3.5' 
    do
      for mBin in '8,14'
      do
        for corrW in '1,1'
        do
          root -l -q -b 'makeV2Hist.C('$cBin','$pt','$yLow','$yHigh','$muPt','$mBin','$isMC','$dimusign','$corrW','$HFSel')'
        done
      done
    done
  done
done

for pt in '0,50' '0,3' '3,6' '6,15'
do
  for cBin in '10,120'
  do
    for rap in '0,2.4' '0,1.2' '1.2,2.4' 
    do
      for mBin in '8,14'
      do
        for corrW in '1,1'
        do
          root -l -q -b 'makeV2Hist.C('$cBin','$pt','$rap',3.5,'$mBin','$isMC','$dimusign','$corrW','$HFSel')'
        done
      done
    done
  done
done



