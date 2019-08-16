#void getEfficiency_syst( float ptLow = 0.0, float ptHigh = 50.0, float yLow = 0.0, float yHigh = 2.4, 
#          int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = true, 
#          int fMuId = -1, int fInnTrk=-1, int fTrig = -1, int state=1)

ptLow=0
ptHigh=50
TnP=1
PtWeight=1

root -l -b <<EOF
.L getEfficiency_syst++
.q
EOF

echo getEfficiency_syst...

for state in '1' '2'
do
  for Idx in '-1' '-2' '1' '2'
  do
    root -l -q -b 'getEfficiency_syst.C('$ptLow','$ptHigh',0,2.4,0,180,'$TnP','$PtWeight','$Idx',0,0,'$state')'
    root -l -q -b 'getEfficiency_syst.C('$ptLow','$ptHigh',0,2.4,0,180,'$TnP','$PtWeight',0,'$Idx',0,'$state')'
    root -l -q -b 'getEfficiency_syst.C('$ptLow','$ptHigh',0,2.4,0,180,'$TnP','$PtWeight',0,0,'$Idx','$state')'
  done
done
