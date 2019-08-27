dimusign=true
accW=true
effW=true
isMC=false
HFSel=0

root -l -b <<EOF
.L doSimultaneousV2MassFit_SysEff_AllState_hckim.C++
.q
EOF
rm *.so *.d *.pcm

echo make hist... 

rap='0,2.4'
SiMuPtCut='3.5'
mBin='8,14'

low10sf='0.5'
low11sf='0.5'
low12sf='0.5'
high10sf='2.5'
high11sf='2.5'
high12sf='2.5'

c_='0.031'
c1_='-0.13142884'
c2_='-0.0114097'
c3_='0.021'

#void doSimultaneousV2MassFit_SysEff_AllState_hckim(int cLow = 20, int cHigh = 180,
#                      float ptLow = 0, float ptHigh = 50,
#                      float yLow = 0, float yHigh=2.4,
#                      float SiMuPtCut = 3.5, float massLow = 8, float massHigh =14, bool dimusign=true, int ibkg_vn_sel = fpol2,
#                      double low10sf=0.5, double low11sf=0.5, double low12sf=0.5,
#                      double high10sf=2.5, double high11sf=2.5, double high12sf=2.5,
#                      double c_ = 0.031, double c1_ = -0.13142884, double c2_ = -0.0114097, double c3_ = 0.021)

#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,180,0,3,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+3.0,'$high11sf'+3.0,'$high12sf'+3.0,'$c_'+0.02,'$c1_'+0.01,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,180,3,6,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+3.0,'$high11sf'+3.0,'$high12sf'+3.0,'$c_'+0.3,'$c1_','$c2_','$c3_')' ### fitstatus = 4 
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,180,6,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.11,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,180,0,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.10,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,100,0,3,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+0.3,'$high11sf'+0.3,'$high12sf'+0.3,'$c_'+0.0,'$c1_','$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,100,3,6,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.05,'$c2_'-0.05,'$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,100,6,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.05,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,100,0,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_','$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(100,180,0,3,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+2.0,'$high11sf'+2.0,'$high12sf'+3.0,'$c_','$c1_','$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(100,180,3,6,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf'+3.0,'$c_','$c1_','$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(100,180,6,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf'+4.5,'$c_'-0.05,'$c1_'+0.05,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(100,180,0,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.09,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(0,20,0,3,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+3.0,'$high11sf'+3.0,'$high12sf'+3.0,'$c_','$c1_'+0.13,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(0,20,3,6,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+4.0,'$high11sf'+4.0,'$high12sf'+4.0,'$c_'-0.05,'$c1_'+0.05,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(0,20,6,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_','$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(0,20,0,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_','$c2_','$c3_')'

#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,0,3,0,2.4,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.12,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,3,6,0,2.4,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.05,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,6,15,0,2.4,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.11,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,0,3,0,1.2,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+5.0,'$high11sf'+5.0,'$high12sf'+5.0,'$c_','$c1_'+0.10,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,3,6,0,1.2,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_'+0.09,'$c2_'-0.05,'$c3_')' ### failed
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,6,15,0,1.2,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_'+0.01,'$c1_'+0.105,'$c2_','$c3_')' ### failed
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,0,3,1.2,2.4,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+5.0,'$high11sf'+5.0,'$high12sf'+5.0,'$c_','$c1_'+0.05,'$c2_','$c3_')'
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,3,6,1.2,2.4,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf'+5.0,'$high11sf'+5.0,'$high12sf'+5.0,'$c_','$c1_'+0.00,'$c2_'+0.00,'$c3_')' ### failed
#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(10,120,6,15,1.2,2.4,'$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_','$c1_','$c2_','$c3_')'

#root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(20,60,0,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf','$c_'+0.01,'$c1_','$c2_','$c3_')'
root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C(60,100,0,50,'$rap','$SiMuPtCut','$mBin',true,fpol2,'$low10sf','$low11sf','$low12sf','$high10sf','$high11sf','$high12sf'+3.0,'$c_'+0.01,'$c1_','$c2_'+0.05,'$c3_')'

#for pt in '0,3' '0,50' '3,6' '6,50' 
#do
#  for cBin in '20,180' '20,100' '100,180' '0,20' 
#  do
#    for rap in '0,2.4' '0,1.2' '1.2,2.4' 
#    do
#      for SiMuPtCut in '3.5' 
#      do
#        for mBin in '8,14'
#        do
# 			  root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C('$cBin','$pt','$rap','$SiMuPtCut','$mBin',true,fpol2)'
#        done
#      done
#    done
#  done
#done

#for pt in '0,50' '0,3' '3,6' '6,15'
#do
#  for cBin in '10,120'
#  do
#    for rap in '0,2.4' '0,1.2' '1.2,2.4' 
#    do
#	   for SiMuPtCut in '3.5'
#      do
#        for mBin in '8,14'
#        do
# 			   root -l -q -b 'doSimultaneousV2MassFit_SysEff_AllState_hckim.C('$cBin','$pt','$rap','$SiMuPtCut','$mBin',true,fpol2)'
#        done
#      done
#    done
#  done
#done
