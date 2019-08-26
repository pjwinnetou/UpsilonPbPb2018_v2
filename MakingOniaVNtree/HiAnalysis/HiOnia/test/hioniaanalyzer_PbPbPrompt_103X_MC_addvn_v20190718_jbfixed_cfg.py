import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
###//from VNA start
import os
import sys 
ivars = VarParsing.VarParsing('standard')

ivars.register ('lumifile',
#                'json_DCSONLY_HI.txt',
#                'Cert_326381-327560_HI_PromptReco_Collisions18_JSON_MuonPhys.txt',
                 'Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="lumi file")

ivars.register ('offset',
                'offset_PbPb2018_1_600000.root',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="offset file")

ivars.register ('dbfile',
                'HeavyIonRPRcd_PbPb2018_offline.db',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="dbfile file")

ivars.register ('eff',
                'NULL',
                mult=ivars.multiplicity.singleton,
                mytype=ivars.varType.string,
                info="efficiency file")

ivars.parseArguments()

###//from VNA end


#----------------------------------------------------------------------------

# Setup Settings for ONIA TREE: PbPb 2018

HLTProcess     = "HLT" # Name of HLT process
isMC           = True # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "GlbTrk" # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker), GlbOrTrk, TwoGlbAmongThree (which requires two isGlobal for a trimuon, and one isGlobal for a dimuon) are available
applyEventSel  =True # Only apply Event Selection if the required collections are present
OnlySoftMuons  = False # Keep only isSoftMuon's (without highPurity, and without isGlobal which should be put in 'muonSelection' parameter) from the beginning of HiSkim. If you want the full SoftMuon selection, set this flag false and add 'isSoftMuon' in lowerPuritySelection. In any case, if applyCuts=True, isSoftMuon is required at HiAnalysis level for muons of selected dimuons.
applyCuts      = False # At HiAnalysis level, apply kinematic acceptance cuts + identification cuts (isSoftMuon (without highPurity) or isTightMuon, depending on TightGlobalMuon flag) for muons from selected di(tri)muons + hard-coded cuts on the di(tri)muon that you would want to add (but recommended to add everything in LateDimuonSelection, applied at the end of HiSkim)
SumETvariables = True  # Whether to write out SumET-related variables
SofterSgMuAcceptance = False # Whether to accept muons with a softer acceptance cuts than the usual (pt>3.5GeV at central eta, pt>1.5 at high |eta|). Applies when applyCuts=True
doTrimuons     = False # Make collections of trimuon candidates in addition to dimuons, and keep only events with >0 trimuons
atLeastOneCand = False # Keep only events that have one selected dimuon (or at least one trimuon if doTrimuons = true). BEWARE this can cause trouble in .root output if no event is selected by onia2MuMuPatGlbGlbFilter!
OneMatchedHLTMu = -1   # Keep only di(tri)muons of which the one(two) muon(s) are matched to the HLT Filter of this number. You can get the desired number in the output of oniaTree. Set to -1 for no matching.
#############################################################################
keepExtraColl  = False # General Tracks + Stand Alone Muons + Converted Photon collections
saveHLT        = False # wheter to save the HLT trees or not 
#----------------------------------------------------------------------------

# Print Onia Tree settings:
print( " " )
print( "[INFO] Settings used for ONIA TREE: " )
print( "[INFO] isMC                 = " + ("True" if isMC else "False") )
print( "[INFO] applyEventSel        = " + ("True" if applyEventSel else "False") )
print( "[INFO] keepExtraColl        = " + ("True" if keepExtraColl else "False") )
print( "[INFO] SumETvariables       = " + ("True" if SumETvariables else "False") )
print( "[INFO] SofterSgMuAcceptance = " + ("True" if SofterSgMuAcceptance else "False") )
print( "[INFO] muonSelection        = " + muonSelection )
print( "[INFO] onlySoftMuons        = " + ("True" if OnlySoftMuons else "False") )
print( "[INFO] doTrimuons           = " + ("True" if doTrimuons else "False") )
print( "[INFO] saveHLT              = " + ("True" if saveHLT else "False") )
print( " " )

# set up process
process = cms.Process("HIOniaVN")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "Oniatree_addvn_GFversion_MC.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles =[
#  '/store/user/anstahll/Dilepton/MC/Embedded/JPsiMM_5p02TeV_TuneCP5_Embd_RECO_20190326/JPsiMM_5p02TeV_TuneCP5_Embd/JPsiMM_5p02TeV_TuneCP5_Embd_RECO_20190326/190327_060527/0008/HIN-HINPbPbAutumn18DRHIMix-00008_step2_8495.root',
#  '/store/user/anstahll/Dilepton/MC/Embedded/JPsiMM_5p02TeV_TuneCP5_Embd_RECO_20190326/JPsiMM_5p02TeV_TuneCP5_Embd/JPsiMM_5p02TeV_TuneCP5_Embd_RECO_20190326/190327_060527/0008/HIN-HINPbPbAutumn18DRHIMix-00008_step2_8496.root'
   '/store/himc/HINPbPbAutumn18DR/Upsilon1S_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11_ext1-v1/110000/0F1FDD00-E1FA-FD4B-BD73-1EE538C7D996.root'
]
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

triggerList    = {
		# Double Muon Trigger List
		'DoubleMuonTrigger' : cms.vstring(
			"HLT_HIL1DoubleMuOpen_v1",#0
			"HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v1", #1
			"HLT_HIL1DoubleMuOpen_Centrality_50_100_v1", #2
			"HLT_HIL1DoubleMu10_v1", #3
			"HLT_HIL2_L1DoubleMu10_v1",#4
			"HLT_HIL3_L1DoubleMu10_v1", #5
			"HLT_HIL2DoubleMuOpen_v1", #6
			"HLT_HIL3DoubleMuOpen_v1", #7
			"HLT_HIL3DoubleMuOpen_M60120_v1", #8
			"HLT_HIL3DoubleMuOpen_JpsiPsi_v1", #9
			"HLT_HIL3DoubleMuOpen_Upsi_v1", #10
			"HLT_HIL3Mu0_L2Mu0_v1", #11
			"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",#12
			"HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1",#13
			"HLT_HIL3Mu3_L1TripleMuOpen_v1",#14
			"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1_L1step",#15
			"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1_L2step",#16
			"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1_L3step",#17
                        ),
		# Double Muon Filter List
		'DoubleMuonFilter'  : cms.vstring(
			"hltL1fL1sL1DoubleMuOpenL1Filtered0",
			"hltL1fL1sL1DoubleMuOpenOSCentrality40100L1Filtered0",
			"hltL1fL1sL1DoubleMuOpenCentrality50100L1Filtered0",
			"hltL1fL1sL1DoubleMu10L1Filtered0",
			"hltL2fL1sL1DoubleMu10L1f0L2Filtered0",
			"hltDoubleMuOpenL1DoubleMu10Filtered",
			"hltL2fL1sL1DoubleMuOpenL1f0L2Filtered0",
			"hltL3fL1DoubleMuOpenL3Filtered0",
			"hltL3fL1DoubleMuOpenL3FilteredM60120",
			"hltL3fL1DoubleMuOpenL3FilteredPsi",
			"hltL3fL1DoubleMuOpenL3FilteredUpsi",
			"hltL3f0L3Mu0L2Mu0Filtered0",
                        "hltL3f0L3Mu0L2Mu0DR3p5FilteredNHitQ10M1to5",
			"hltL3f0L3Mu2p5NHitQ10L2Mu2FilteredM7toinf",
                        "hltL3fL1sL1DoubleMuOpenL1fN3L2f0L3Filtered3",
                        "hltL1fL1sL1DoubleMuOpenMAXdR3p5L1Filtered0",#L1 step for Jpsi trigger
                        "hltL2fDoubleMuOpenL2DR3p5PreFiltered0",#L2 step
                        "hltL3f0L3Mu0L2Mu0DR3p5FilteredNHitQ10M1to5"#"hltL3f0DR3p5L3FilteredNHitQ10"#L3 step
			),
                # Single Muon Trigger List
                'SingleMuonTrigger' : cms.vstring(
                        "HLT_HIL1MuOpen_Centrality_70_100_v1",
                        "HLT_HIL1MuOpen_Centrality_80_100_v1",
                        "HLT_HIL2Mu3_NHitQ15_v1",
                        "HLT_HIL2Mu5_NHitQ15_v1",
                        "HLT_HIL2Mu7_NHitQ15_v1",
                        "HLT_HIL3Mu3_NHitQ10_v1",
                        "HLT_HIL3Mu5_NHitQ10_v1",
                        "HLT_HIL3Mu7_NHitQ10_v1",
                        "HLT_HIL3Mu12_v1",
                        "HLT_HIL3Mu15_v1",
                        "HLT_HIL3Mu20_v1",
                        "HLT_HIL2Mu3_NHitQ15_v2",
                        "HLT_HIL2Mu5_NHitQ15_v2",
                        "HLT_HIL2Mu7_NHitQ15_v2",
                        "HLT_HIL3Mu3_NHitQ10_v2",
                        "HLT_HIL3Mu5_NHitQ10_v2",
                        "HLT_HIL3Mu7_NHitQ10_v2",
                        "HLT_HIL3Mu12_v2",
                        "HLT_HIL3Mu15_v2",
                        "HLT_HIL3Mu20_v2",
                        "HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_L2Filter_v1",
                        "HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_L3Filter_v1",
			),
	        # Single Muon Filter List
	        'SingleMuonFilter'  : cms.vstring(
                        "hltL1fL1sL1MuOpenCentrality70100L1Filtered0",
                        "hltL1fL1sL1MuOpenCentrality80100L1Filtered0",
                        "hltL2fL1sMu3OpenL1f0L2Filtered3NHitQ15",
                        "hltL2fL1sMu3OpenL1f0L2Filtered5NHitQ15",
                        "hltL2fL1sMu3OpenL1f0L2Filtered7NHitQ15",
                        "hltL3fL1sL1SingleMu3OpenL1f0L2f0L3Filtered3NHitQ10",
                        "hltL3fL1sL1SingleMu3OpenL1f0L2f0L3Filtered5NHitQ10",
                        "hltL3fL1sL1SingleMu3OpenL1f0L2f0L3Filtered7NHitQ10",
                        "hltL3fL1sL1SingleMu3OpenL1f7L2f0L3Filtered12",
                        "hltL3fL1sL1SingleMu3OpenL1f7L2f0L3Filtered15",
                        "hltL3fL1sL1SingleMu3OpenL1f7L2f0L3Filtered20",
                        "hltL2fL1sMuOpenL1f0L2Filtered3NHitQ15",
                        "hltL2fL1sMuOpenL1f0L2Filtered5NHitQ15",
                        "hltL2fL1sMuOpenL1f0L2Filtered7NHitQ15",
                        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered3NHitQ10",
                        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered5NHitQ10",
                        "hltL3fL1sL1SingleMuOpenL1f0L2f0L3Filtered7NHitQ10",
                        "hltL3fL1sL1SingleMuOpenL1f7L2f0L3Filtered12",
                        "hltL3fL1sL1SingleMuOpenL1f7L2f0L3Filtered15",
                        "hltL3fL1sL1SingleMuOpenL1f7L2f0L3Filtered20",
								"hltL2fDoubleMuOpenL2PreFiltered0",
                        "hltL3f0L3Filtered2p5NHitQ10",
			)
                }

## Global tag
if isMC:
  globalTag = '103X_upgrade2018_realistic_HI_v11'
else:
  globalTag = '103X_dataRun2_Prompt_v3'

#----------------------------------------------------------------------------

# load the Geometry and Magnetic Field
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

###///
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.load("HeavyIonsAnalysis.VNAnalysis/vnanalyzer_cfi")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("CondCore.CondDB.CondDB_cfi")
process.load("HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi")
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load("HeavyIonsAnalysis.Configuration.analysisFilters_cff")
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
###///

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')
### For Centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
#process.centralityBin.Centrality = cms.InputTag("hiCentrality")
#process.centralityBin.centralityVariable = cms.string("HFtowers")
print('\n\033[31m~*~ USING CENTRALITY TABLE FOR HYDJET DRUM5EV8 TUNE~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_v1032x02_mc"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
###///
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.hiCentrality.srcTracks = cms.InputTag("generalTracks")
process.hiCentrality.srcVertex = cms.InputTag("offlinePrimaryVertices")
###///
#----------------------------------------------------------------------------

# For OniaTree Analyzer
from HiAnalysis.HiOnia.oniaTreeAnalyzer_cff import oniaTreeAnalyzer
oniaTreeAnalyzer(process,
                 muonTriggerList=triggerList,# HLTProName=HLTProcess,
                 muonSelection=muonSelection, useL1Stage2=True, isMC=isMC, outputFileName=options.outputFile, doTrimu=doTrimuons)

#process.onia2MuMuPatGlbGlb.dimuonSelection       = cms.string("8 < mass && mass < 14 && charge==0 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
#process.onia2MuMuPatGlbGlb.lowerPuritySelection  = cms.string("")
#process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("") ## No need to repeat lowerPuritySelection in there, already included
if applyCuts:
  process.onia2MuMuPatGlbGlb.LateDimuonSel         = cms.string("userFloat(\"vProb\")>0.01")
process.onia2MuMuPatGlbGlb.onlySoftMuons         = cms.bool(OnlySoftMuons)
process.hionia.minimumFlag      = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.useGeTracks      = cms.untracked.bool(keepExtraColl) #for Reco_trk_*
process.hionia.fillRecoTracks   = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.CentralitySrc    = cms.InputTag("hiCentrality")
process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
process.hionia.srcTracks        = cms.InputTag("generalTracks")
#process.hionia.muonLessPV       = cms.bool(False)
process.hionia.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")
process.hionia.genParticles     = cms.InputTag("genParticles")
process.hionia.SofterSgMuAcceptance = cms.bool(SofterSgMuAcceptance)
process.hionia.SumETvariables   = cms.bool(SumETvariables)
process.hionia.applyCuts        = cms.bool(applyCuts)
process.hionia.AtLeastOneCand   = cms.bool(atLeastOneCand)
process.hionia.OneMatchedHLTMu  = cms.int32(OneMatchedHLTMu)
process.hionia.oniaPDG           = cms.int32(553)

#----------------------------------------------------------------------------

# For the primary vertex recovery
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")

#----------------------------------------------------------------------------

# For HLTBitAnalyzer
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltbitanalysis.HLTProcessName              = HLTProcess
process.hltbitanalysis.hltresults                  = cms.InputTag("TriggerResults","",HLTProcess)
process.hltbitanalysis.l1tAlgBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.l1tExtBlkInputTag           = cms.InputTag("hltGtStage2Digis","",HLTProcess)
process.hltbitanalysis.gObjectMapRecord            = cms.InputTag("hltGtStage2ObjectMap","",HLTProcess)
process.hltbitanalysis.gmtStage2Digis              = cms.string("hltGtStage2Digis")
process.hltbitanalysis.caloStage2Digis             = cms.string("hltGtStage2Digis")
process.hltbitanalysis.UseL1Stage2                 = cms.untracked.bool(True)
process.hltbitanalysis.getPrescales                = cms.untracked.bool(False)
process.hltbitanalysis.getL1InfoFromEventSetup     = cms.untracked.bool(False)
process.hltbitanalysis.UseTFileService             = cms.untracked.bool(True)
process.hltbitanalysis.RunParameters.HistogramFile = cms.untracked.string(options.outputFile)
process.hltbitanalysis.RunParameters.isData        = cms.untracked.bool(not isMC)
process.hltbitanalysis.RunParameters.Monte         = cms.bool(isMC)
process.hltbitanalysis.RunParameters.GenTracks     = cms.bool(False)
if (HLTProcess == "HLT") :
	process.hltbitanalysis.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO")
	process.hltbitanalysis.gmtStage2Digis    = cms.string("gtStage2Digis")
	process.hltbitanalysis.caloStage2Digis   = cms.string("gtStage2Digis")

##----------------------------------------------------------------------------

# For HLTObject Analyzer
process.load("HeavyIonsAnalysis.EventAnalysis.hltobject_cfi")
process.hltobject.processName = cms.string(HLTProcess)
process.hltobject.treeName = cms.string(options.outputFile)
process.hltobject.loadTriggersFromHLT = cms.untracked.bool(False)
process.hltobject.triggerNames = triggerList['DoubleMuonTrigger'] + triggerList['SingleMuonTrigger']
process.hltobject.triggerResults = cms.InputTag("TriggerResults","",HLTProcess)
process.hltobject.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","",HLTProcess)

#---------------------------------------------------------------------------

#For the main analysis list
process.oniaTreeAna.replace(process.hionia, process.centralityBin * process.hionia )

if saveHLT:
  process.oniaTreeAna = cms.Path(process.offlinePrimaryVerticesRecovery * process.hltbitanalysis * process.hltobject * process.oniaTreeAna )
else:
  process.oniaTreeAna = cms.Path(process.offlinePrimaryVerticesRecovery * process.oniaTreeAna )

if atLeastOneCand:
  process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter)
  if doTrimuons:
    process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlbFilter3mu * process.onia2MuMuPatGlbGlb)

if applyEventSel:
  process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
  process.load('HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi')
  process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
###///  process.oniaTreeAna.replace(process.hionia, process.centralityBin * process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter * process.hionia )
  process.oniaTreeAna.replace(process.hionia, process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter * process.hiEvtPlane * process.hiEvtPlaneFlat* process.hionia * process.vnanalyzer)


#---------------------------------------------------------------------------


#----------------------------------------------------------------------------
#Options:
process.source = cms.Source("PoolSource",
#process.source = cms.Source("NewEventStreamFileReader", # for streamer data
		fileNames = cms.untracked.vstring( options.inputFiles ),
###///
      inputCommands=cms.untracked.vstring(
         'keep *',
         'drop *_hiEvtPlane_*_*'
      )
###///
		)
process.TFileService = cms.Service("TFileService",
		fileName = cms.string( options.outputFile )
		)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#######################Offline Primary Vertices
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"


###///
process.MessageLogger.cerr.FwkReport.reportEvery=1000

process.CondDB.connect = "sqlite_file:"+ivars.dbfile
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
#                                                                  tag = cms.string('HeavyIonRPRcd')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')

import FWCore.PythonUtilities.LumiList as LumiList
goodLumiSecs = LumiList.LumiList(filename = ivars.lumifile ).getCMSSWString().split(',')

# MinBias trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltSelect = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltSelect.HLTPaths = [
        "HLT_HIMinimumBias_*",
    ]
process.hltSelect.andOr = cms.bool(True)
process.hltSelect.throw = cms.bool(False)

# Event Selection
process.primaryVertexFilter.src = cms.InputTag("offlinePrimaryVertices")
#process.towersAboveThreshold.minimumE = cms.double(4.0)
#process.eventSelection = cms.Sequence(
#    process.hfCoincFilter2
#    + process.clusterCompatibilityFilter
#    + process.primaryVertexFilter
#)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.hiEvtPlane.trackTag = cms.InputTag("generalTracks")
process.hiEvtPlane.vertexTag = cms.InputTag("offlinePrimaryVertices")
process.hiEvtPlane.loadDB = cms.bool(True)
process.hiEvtPlane.useNtrk = cms.untracked.bool(False)
process.hiEvtPlane.caloCentRef = cms.double(-1)
process.hiEvtPlane.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRef = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlinePrimaryVertices")
process.hiEvtPlaneFlat.useNtrk = cms.untracked.bool(False)
process.vnanalyzer.trackTag_ = cms.InputTag("generalTracks")
process.vnanalyzer.vertexTag_ = cms.InputTag("offlinePrimaryVertices")
process.vnanalyzer.useNtrk = cms.untracked.bool(False)
process.vnanalyzer.offsetFile = cms.untracked.string( ivars.offset )
process.vnanalyzer.effFile = cms.untracked.string( ivars.eff )
process.vnanalyzer.EPLevel = cms.untracked.int32(2)
process.vnanalyzer.Recenter = cms.untracked.bool(True)
process.vnanalyzer.chi2_ = cms.untracked.double(9.0)
process.vnanalyzer.dzdzerror_pix_ = cms.untracked.double(6.0)
process.vnanalyzer.dzdzerror_ = cms.untracked.double(3.0)
process.vnanalyzer.d0d0error_ = cms.untracked.double(3.0)
process.vnanalyzer.pterror_ = cms.untracked.double(0.1)
process.vnanalyzer.flatnvtxbins = cms.int32(10);
process.vnanalyzer.flatminvtx = cms.double(-25);
process.vnanalyzer.flatdelvtx = cms.double(5.0);
process.vnanalyzer.minrun_ = cms.untracked.int32(1);##326500
process.vnanalyzer.maxrun_ = cms.untracked.int32(1000000);##328500
process.vnanalyzer.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")

#process.schedule  = cms.Schedule( process.oniaTreeAna )
###///
