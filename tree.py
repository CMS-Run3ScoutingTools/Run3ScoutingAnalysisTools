import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2023C/ScoutingPFRun3/RAW/v1/000/367/696/00000/47d704c0-3d61-47c1-b4c4-abf9dd43f316.root',
        #'/store/data/Run2023C/ScoutingPFRun3/RAW/v1/000/367/232/00000/23176b63-d9f1-4fc2-b5f6-ab93f81c835e.root',
        #'/store/data/Run2023C/ScoutingPFRun3/RAW/v1/000/367/232/00000/4144e595-1af4-4f32-8b8e-9c104289634c.root'
 )
)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scout.root")
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '') # Run 3 2022
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v2', '') # Run 3 2023

L1Info = ["L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2", "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu8_SQ"]

process.scoutingTree = cms.EDAnalyzer('ScoutingTreeMakerRun3',
                                      muons             = cms.InputTag("hltScoutingMuonPacker"),
                                      electrons         = cms.InputTag("hltScoutingEgammaPacker"),
                                      photons           = cms.InputTag("hltScoutingEgammaPacker"),
                                      pfcands           = cms.InputTag("hltScoutingPFPacker"),
                                      pfjets            = cms.InputTag("hltScoutingPFPacker"),
                                      tracks            = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices   = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
                                      pfMet             = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
                                      pfMetPhi          = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
                                      rho               = cms.InputTag("hltScoutingPFPacker","rho"),
                                      triggerAlias = cms.vstring(["Run3_ScoutingDoubleMuNoVtx"]),
                                      triggerSelection = cms.vstring(["DST_Run3_DoubleMu3_PFScoutingPixelTracking_v*"]), # DST_Run3_PFScoutingPixelTracking_v* pre Run367622
                                      triggerConfiguration = cms.PSet(
                                          hltResults            = cms.InputTag('TriggerResults','','HLT'),
                                          l1tResults            = cms.InputTag('','',''),
                                          l1tIgnoreMaskAndPrescale = cms.bool(False),
                                          throw                 = cms.bool(True),
                                          usePathStatus = cms.bool(False),
                                      ),
                                      doL1 = cms.bool(True),
                                      doTriggerObjects = cms.bool(False),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1Seeds           = cms.vstring(L1Info)
                                      )

process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)
