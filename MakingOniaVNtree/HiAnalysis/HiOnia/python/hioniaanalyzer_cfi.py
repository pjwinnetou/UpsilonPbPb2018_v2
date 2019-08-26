import FWCore.ParameterSet.Config as cms

hionia = cms.EDAnalyzer('HiOniaAnalyzer',
                        srcMuon          = cms.InputTag("patMuonsWithTrigger"),
                        srcMuonNoTrig    = cms.InputTag("patMuonsWithoutTrigger"),
                        srcDimuon        = cms.InputTag("onia2MuMuPatGlbGlb"),
                        srcTracks        = cms.InputTag("hiGeneralTracks"),
                        genParticles     = cms.InputTag("genParticles"),
                        EvtPlane         = cms.InputTag("hiEvtPlane","recoLevel"),
                        primaryVertexTag = cms.InputTag("hiSelectedVertex"),
                        
                        triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),

                        CentralitySrc    = cms.InputTag("hiCentrality"),
                        CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
                        
                        #-- Reco Details
                        useBeamSpot = cms.bool(False),
                        useRapidity = cms.bool(True),
                        
                        #--
                        maxAbsZ = cms.double(24.0),
                        
                        pTBinRanges = cms.vdouble(0.5, 3.0, 6.0, 8.0, 10.0, 15.0, 35.0),
                        etaBinRanges = cms.vdouble(0.0, 2.5),
                        centralityRanges = cms.vdouble(20,40,100),
                        
                        onlyTheBest = cms.bool(False),		
                        applyCuts = cms.bool(True),
			selTightGlobalMuon = cms.bool(False), 
                        storeEfficiency = cms.bool(False),
                        SofterSgMuAcceptance = cms.bool(False),
                        SumETvariables = cms.bool(True),
                        OneMatchedHLTMu = cms.int32(-1),
                        storeSameSign = cms.bool(False),
                        AtLeastOneCand = cms.bool(False),
                        doTrimuons = cms.bool(False),
                        removeSignalEvents = cms.untracked.bool(False),
                        removeTrueMuons = cms.untracked.bool(False),
                        
                        muonLessPV = cms.bool(False),
                        
                        #-- Gen Details
                        BcPDG = cms.int32(541),
                        oniaPDG = cms.int32(443),
                        muonSel = cms.string("GlbGlb"),
                        isHI = cms.untracked.bool(True),
                        isPA = cms.untracked.bool(False),
                        isMC = cms.untracked.bool(False),
                        isPromptMC = cms.untracked.bool(False),
                        useEvtPlane = cms.untracked.bool(True),
                        useGeTracks = cms.untracked.bool(False),

                        #-- Histogram configuration
                        combineCategories = cms.bool(False),
                        fillRooDataSet = cms.bool(False),
                        fillTree = cms.bool(True),
                        fillHistos = cms.bool(True),
                        minimumFlag = cms.bool(False),
                        fillSingleMuons = cms.bool(True),
                        fillRecoTracks = cms.bool(False),
                        histFileName = cms.string("Jpsi_Histos.root"),		
                        dataSetName = cms.string("Jpsi_DataSet.root"),
                        
                        #--
                        dblTriggerPathNames    = cms.vstring(),
                        dblTriggerFilterNames = cms.vstring(),
                        sglTriggerPathNames    = cms.vstring(),
                        sglTriggerFilterNames = cms.vstring()
                        )
