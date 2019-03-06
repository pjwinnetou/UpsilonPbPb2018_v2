#void makeV2Hist(int cLow = 20, cHigh = 120,
#                float ptLow = 0, float ptHigh = 6,
#                float yLow = 0, float yHigh = 2.4,
#                float SiMuPtCut = 4.0,
#                float massLow = 8, float massHigh = 14,


SiMuPtCut='3.5'
cLow='0'
cHigh='200'
ptLow='0'
ptHigh='30'
yLow='0'
yHigh='2.4'
massLow='7'
massHigh='14'

#root -q -b makeV2Hist.C'('$cLow','$cHigh','$ptLow','$ptHigh','$yLow','$yHigh','$SiMuPtCut','$massLow','$massHigh')'

outputdir='plots/MassV2_190306'
mkdir -p $outputdir
cp commonUtility.h cutsAndBinUpsilonV2.h HiEvtPlaneList.h Style_jaebeom.h tdrstyle.C CMS_lumi_v2mass.* $outputdir 

root -l -b <<EOF
.L makeV2Hist.C++
.q
EOF

echo make hist... 

for pt in '0,30' '0,3' '0,6' '3,6' '6,30' 
do
  for cBin in '20,60' '20,120' '60,100' '100,200'
  do
    for muPt in '3.5'
    do
      for mBin in '8,14'
      do
        ./condor_root.sh $outputdir 'makeV2Hist.C('$cBin','$pt','$yLow','$yHigh','$muPt','$mBin')'
      done
    done
  done
done
