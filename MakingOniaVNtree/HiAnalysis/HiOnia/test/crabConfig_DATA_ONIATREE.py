from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = "HIDoubleMu_Run2018A_PromptAOD_OniaTree_Run_326381_326943"
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "hioniaanalyzer_PbPbPrompt_103X_DATA_cfg.py"
config.JobType.maxMemoryMB = 2500         # request high memory machines.
config.JobType.maxJobRuntimeMin = 2750    # request longer runtime, ~48 hours.

## software : CMSSW_10_3_1

config.section_("Data")
config.Data.userInputFiles = open('PR_DoubleMu_381_943.txt').readlines()
#config.Data.userInputFiles = open('Inputfiles_DoubleMu_PromptAOD_test.txt').readlines() ##CHECK IT!!## 
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1

config.Data.publication = False
#config.Data.outputPrimaryDataset = "AOD"
#config.Data.outputDatasetTag = '326822'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD' ##CHECK IT!!##
config.Data.outLFNDirBase = '/store/group/phys_heavyions/dileptons/goni'


config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
#config.Site.whitelist = ["T2_CH_CERN"]

#config.section_("Debug")
#config.Debug.extraJDL = ["+CMS_ALLOW_OVERFLOW=False"]
