import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2025C/ScoutingPFRun3/HLTSCOUT/v1/000/392/925/00000/b95d5cc9-62b2-4b3b-a0f9-d0d79b52a85d.root'
 )
)

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scout.root")
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_HLT_v1', '') # Run 3 2025
                                                                                                                                                                      
L1Info = ["L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu8_SQ","L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0_Upt6_SQ_er2p0","L1_DoubleMu0_Upt7_SQ_er2p0","L1_DoubleMu0_Upt8_SQ_er2p0", "L1_DoubleMu0_Upt15_Upt7", "L1_DoubleMu0_Upt6_IP_Min1_Upt4", "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6"]

process.scoutingTree = cms.EDAnalyzer('ScoutingTreeMakerRun3',
                                      muonsNoVtx        = cms.InputTag("hltScoutingMuonPackerNoVtx"), # hltScoutingMuonPackerNoVtx && hltScoutingMuonPackerVtx 
                                      muonsVtx          = cms.InputTag("hltScoutingMuonPackerVtx"), 
                                      electrons         = cms.InputTag("hltScoutingEgammaPacker"),
                                      photons           = cms.InputTag("hltScoutingEgammaPacker"),
                                      pfcands           = cms.InputTag("hltScoutingPFPacker"),
                                      pfjets            = cms.InputTag("hltScoutingPFPacker"),
                                      tracks            = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices   = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices = cms.InputTag("hltScoutingMuonPackerNoVtx","displacedVtx"),
                                      pfMet             = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
                                      pfMetPhi          = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
                                      rho               = cms.InputTag("hltScoutingPFPacker","rho"),
                                      triggerAlias = cms.vstring(["Run3_ScoutingDoubleMuNoVtx", "Run3_ScoutingDoubleMuVtx"]),
                                      triggerSelection = cms.vstring(["DST_PFScouting_DoubleMuonNoVtx_v*", "DST_PFScouting_DoubleMuonVtx_v*"]), 
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
