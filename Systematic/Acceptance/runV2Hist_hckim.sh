dimusign=true
accW=true
effW=true
isMC=false
HFSel=0

root -l -b <<EOF
.L doSimultaneousV2MassFit_SysAcc_AllState.C++
.q
EOF
rm *.so *.d *.pcm

echo make hist... 

for pt in '0,3' '0,50' '3,6' '6,50' 
do
  for cBin in '20,180' '20,100' '100,180' '0,20' 
  do
    for rap in '0,2.4' '0,1.2' '1.2,2.4' 
    do
      for SiMuPtCut in '3.5' 
      do
        for mBin in '8,14'
        do
#          root -l -q -b 'makeV2Hist_hckim.C('$cBin','$pt','$yLow','$yHigh','$muPt','$mBin','$isMC','$dimusign','$corrW','$HFSel')'
#			  root -l -q -b 'makeV2Hist_hckim.C('$cBin','$pt','$rap','$muPt','$mBin','$isMC','$dimusign','$corrW','$HFSel')'
 			  root -l -q -b 'doSimultaneousV2MassFit_SysAcc_AllState.C('$cBin','$pt','$rap','$SiMuPtCut','$mBin',true,fpol2)'
#(int cLow = 20, int cHigh = 180,float ptLow = 0, float ptHigh = 50,yLow = 0, float yHigh=2.4,float SiMuPtCut = 3.5, float massLow = 8, float massHigh =14, bool dimusign=true, int ibkg_vn_sel = fpol2)
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
	   for SiMuPtCut in '3.5'
      do
        for mBin in '8,14'
        do
#           root -l -q -b 'makeV2Hist_hckim.C('$cBin','$pt','$rap',3.5,'$mBin','$isMC','$dimusign','$corrW','$HFSel')'
 			   root -l -q -b 'doSimultaneousV2MassFit_SysAcc_AllState.C('$cBin','$pt','$rap','$SiMuPtCut','$mBin',true,fpol2)'
        done
      done
    done
  done
done



