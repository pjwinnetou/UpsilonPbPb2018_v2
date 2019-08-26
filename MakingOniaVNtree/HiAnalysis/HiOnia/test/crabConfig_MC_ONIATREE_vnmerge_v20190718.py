from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
#config.General.requestName = "HIDoubleMu_Run2018A_PromptAOD_OniaTree_vnmerge_v1_notFinished_Cert_Run_326381_327560"
##config.General.requestName = "HIDoubleMu_Run2018A_Promptv1AOD_OniaTreeGF_vnmerge_v2_notFinished_CertMuonPhys_Run_326381_327560_v20190210"
###config.General.requestName = "Ups1SMM_5p02TeV_TuneCP5_HydjetDrumMB_MC_OniaTreeGF_vnmerge_v2_v20190710"

#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC-v1_OniaTreeGFvnmerge_p20190718"
#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC_ext1-v1_OniaTreeGFvnmerge_p20190718"
#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC-v1_OniaTreeGFvnmerge_p20190718_r2"
#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC_ext1-v1_OniaTreeGFvnmerge_p190718_r2"
#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC_ext1-v1_OniaTreeGFvnmerge_p190718_r3"
#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC-v1_OniaTreeGFvnmerge_p190718_r3"

#config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC_ext1-v1_OniaTreeGFvnmerge_p190718_r4"
config.General.requestName = "Ups1SMM_5p02TeV_pThat-2_TuneCP5_HydjetDrumMB_officialPythia8MC-v1_OniaTreeGFvnmerge_p190718_r4"





config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
#config.JobType.psetName = "hioniaanalyzer_PbPbPrompt_103X_DATA_from_Run327327_vnmerge_v4_cfg.py"
###config.JobType.psetName = "hioniaanalyzer_PbPbPrompt_103X_MC_addvn_v20190410_cfg.py"
config.JobType.psetName = "hioniaanalyzer_PbPbPrompt_103X_MC_addvn_v20190718_jbfixed_cfg.py"
config.JobType.maxMemoryMB = 2500         # request high memory machines.
config.JobType.maxJobRuntimeMin = 2750    # request longer runtime, ~48 hours.
config.JobType.inputFiles = ['HeavyIonRPRcd_PbPb2018_offline.db']
config.JobType.allowUndistributedCMSSW = True

## software : CMSSW_10_3_1

config.section_("Data")
#config.Data.userInputFiles = open('PR_DoubleMu_381_943.txt').readlines()
#config.Data.userInputFiles = open('Inputfiles_DoubleMu_PromptAOD_test.txt').readlines() ##CHECK IT!!## 
##config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v1/AOD'
#config.Data.inputDataset = '/Ups1SMM_5p02TeV_TuneCP5_Embd/anstahll-Ups1SMM_5p02TeV_TuneCP5_Embd_RECO_20190401-5db5dfa073297cb96791f14c622e83e2/USER'
#config.Data.inputDataset = '/Upsilon1S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11_ext1-v1/AODSIM'

config.Data.inputDataset = '/Upsilon1S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM'
###config.Data.inputDataset = '/Upsilon1S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11_ext1-v1/AODSIM'

#config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v2/AOD'
#config.Data.inputDataset = '/HIDoubleMuonPsiPeri/HIRun2018A-PromptReco-v2/AOD'
config.Data.inputDBS = 'global'
###config.Data.splitting = "LumiBased"
#config.Data.lumiMask = 'Cert_326381-327560_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
#config.Data.lumiMask = 'crab_projects/crab_HIDoubleMu_Run2018A_PromptAOD_OniaTree_vnmerge_v1_Cert_Run_326381_327560/results/notFinishedLumis.json'
###config.Data.lumiMask = 'json_DCSONLY_HI.txt'
#config.Data.lumiMask = 'Cert_326381-327560_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
#config.Data.runRange = '326381-326998' ## for 326XXX - Gyeonghwan
#config.Data.runRange = '327004-327560' ## for 327XXX - Hanseul
###config.Data.unitsPerJob = 3 ## After test, if you want please optimize
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1 

config.Data.publication = False
#config.Data.outputPrimaryDataset = "AOD"
#config.Data.outputDatasetTag = '326822'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD' ##CHECK IT!!##
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/dileptons/goni'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbData/326XXX' ## for 326XXX - Gyeonghwan
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbData/327XXX' ## for 327XXX - Hanseul

#config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_pThat_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8' ## for 327XXX - Hanseul
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/gbak/2018PbPbMC/Ups1S_pThat_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8' ## for 327XXX - Hanseul

###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC-v1_p20190718'
###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC_ext1-v1_p20190718'

###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC-v1_p20190718_r2'
###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC_ext1-v1_p20190718_r2'

###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC-v1_p20190718_r3'
###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC_ext1-v1_p20190718_r3'

config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC-v1_p20190718_r4'
###config.Data.outLFNDirBase = '/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC_ext1-v1_p20190718_r4'


config.section_("Site")
#config.Site.storageSite = "T2_CH_CERN"
config.Site.storageSite = "T2_KR_KNU"
##config.Site.whitelist = ["T1_FR_*","T2_FR_*","T3_FR_*","T2_CH_*","T2_CH_CERN_HLT"]
###config.Site.whitelist = ["T2_KR_*","T3_KR_*"]
#config.Site.whitelist = ['T3_KR_KNU"]
###config.Data.ignoreLocality = True

#config.section_("Debug")
#config.Debug.extraJDL = ["+CMS_ALLOW_OVERFLOW=False"]
