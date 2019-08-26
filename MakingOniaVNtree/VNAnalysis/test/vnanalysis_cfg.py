import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
ivars = VarParsing.VarParsing('standard')

ivars.register ('lumifile',
                'json_DCSONLY_HI.txt',
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

process = cms.Process("VNANAL")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
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

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_v3', '')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
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

process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
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

#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 
#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(),
#                             inputCommands=cms.untracked.vstring(
#        'keep *',
#        'drop *_hiEvtPlane_*_*'
#        )
#)


process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/ED38ED56-282F-1F44-8EE1-14C7A88E695C.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/E7874EB8-3B57-F74E-9FE2-B0F5DCF5A30A.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/E10FD51A-3868-7647-BBE6-F37048E63429.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/DEEBA4AA-4BBF-5843-A015-1240D38EA8DF.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/DBD1813E-C390-2C4C-8211-EA4FBB3917AA.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/D603AE46-05C7-894C-AC57-2F3387B3AAD8.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/D53479A7-0758-C746-8204-9735BEF4B46D.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/CA266A9E-A159-E142-8F3C-24768C91F2E0.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/C4C88B44-9BE9-EF45-AA98-CCD543503337.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/C3B1B535-A46C-444D-A63E-3626CBCB9AE9.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/C352710F-F9DD-3F4A-A7E8-477A578C0E72.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/BD5CCF34-D0AE-E841-8F3B-A771243BC3BC.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/B7900E54-B1F9-9A47-A80A-977104FEBAC9.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/A5BA7419-7BE6-0D44-8239-E69932197A15.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/99990A72-C35C-1143-AD1E-64F21A2BA84E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/96FAB6BA-50F7-724D-AB60-74E701D741D4.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/942A5B84-8378-D045-AD36-AA3F075B124C.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/91AAF6AF-0EE7-9345-B14B-BCFB0BCD50C5.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/90A53F37-2CE7-D84D-B774-C3979F845102.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/8C1FA525-9EE3-8E4C-876C-8BB086B706BD.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/889C126C-F0E7-534E-9C28-246AEC2D3193.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/74847CE7-D5F1-A143-9F65-4EC5114F8208.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/74427A76-83A6-BA43-8075-D66F199973CF.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/6DF921F6-CDDF-E049-BF5B-7734FE7D1DA6.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/68B6E95B-64D8-9F43-8086-E99D6EF847AC.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/6077DC41-0575-5F48-A57C-938DD5EB5AF4.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/60322F25-B61E-3146-A35C-540993BDC8EC.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/5DF10A7D-C102-B644-B98C-E75EEDF33B2E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/5CEE0A08-7B90-C54B-9AB5-3659BE8E0BD0.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/5AB1CB27-69D6-0A46-84AB-9835594FB4F2.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/58F879EE-F411-6F40-8576-CA32E17921DA.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/5521A72D-94CB-6242-A807-62CA33209B53.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/4F76FD34-D4F7-C543-B53A-2F28C236655F.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/4C2945A8-0E8D-5C42-9C07-9AD650A3980E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/46E29E17-8DB8-5747-BEC9-376A80F2C1EE.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/42E6DA5B-0033-6542-B2C7-A9A7D40B84E8.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/3C6FE385-1B3C-D24E-9EBC-73C8E407A6AE.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/374FDBBC-ADEF-C94A-A6C6-9A8DC1EAA07E.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/373ADDCA-5EEC-CB47-A019-D1C99D6839A5.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/35D97B5D-6EE2-544B-A605-91A3DD0F8C96.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/2E80D52C-08E8-894E-B882-19B8887B0F56.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/2D56F5CE-0D93-2642-9AA7-A811E22712FA.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/29FD2F95-6702-FA4D-BE14-3D22A7EE7215.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/125A07E2-7D81-9D47-A20B-9098C3402747.root',
        'root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIMinimumBias0/AOD/PromptReco-v1/000/326/534/00000/047760DC-191B-7B4B-AAF0-4BD760443BC1.root'),
                              inputCommands=cms.untracked.vstring(
                                 'keep *',
                                 'drop *_hiEvtPlane_*_*'
                             )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("vnanal.root")
)

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
process.vnanalyzer.minrun_ = cms.untracked.int32(326500);
process.vnanalyzer.maxrun_ = cms.untracked.int32(328500);
process.p = cms.Path(process.hltSelect*process.collisionEventSelectionAODv2*process.centralityBin* process.hiEvtPlane * process.hiEvtPlaneFlat*process.vnanalyzer)
