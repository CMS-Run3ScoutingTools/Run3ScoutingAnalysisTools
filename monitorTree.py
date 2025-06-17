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
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/a69cc2e9-6972-4b45-8fb2-b45f10c3f01e.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/524df273-5459-4fab-beac-f235a04c2091.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/672/00000/370a2ee2-d21a-40e0-a65c-7225d49dec7f.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/672/00000/42b03268-e2d8-4e2f-aa1f-dadd592f071c.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/98338bba-d4e9-40d3-be46-76464e3c4bfc.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/c94b5fe2-573c-4f41-9933-8c962ca7e48f.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/211accd1-c029-4da8-a30a-d87e0540acaa.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/e987ae2f-e44c-4a26-b258-0dcf68c33994.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/700/00000/f90ebb38-48e9-409a-ac81-0790182d4e7f.root",
                                "/store/data/Run2025C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/392/702/00000/04f282ae-dd0f-458b-be72-4bc195598b96.root"
                            ),                                                                 
                            lumisToProcess = cms.untracked.VLuminosityBlockRange('392700:1-392700:1000')
)


process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "rawDataCollector", "", "LHC" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scoutMonitor.root")
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_HLT_v1', '') # Run 3 2025

L1Seeds = ["L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu8_SQ","L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0_Upt6_SQ_er2p0","L1_DoubleMu0_Upt7_SQ_er2p0","L1_DoubleMu0_Upt8_SQ_er2p0", "L1_DoubleMu0_Upt15_Upt7", "L1_DoubleMu0_Upt6_IP_Min1_Upt4"]

process.scoutingTree = cms.EDAnalyzer('ScoutingTreeMakerRun3Monitor',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Seeds),
                                      muons             = cms.InputTag("hltScoutingMuonPackerVtx","","HLT"), #hltScoutingMuonPackerNoVtx
                                      offlineMuons      = cms.untracked.InputTag("slimmedMuons"),
                                      electrons         = cms.InputTag("hltScoutingEgammaPacker"),
                                      offlinePhotons    = cms.untracked.InputTag("slimmedPhotons"),
                                      photons           = cms.InputTag("hltScoutingEgammaPacker"),
                                      pfcands           = cms.InputTag("hltScoutingPFPacker"),
                                      pfjets            = cms.InputTag("hltScoutingPFPacker"),
                                      tracks            = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices   = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices = cms.InputTag("hltScoutingMuonPackerNoVtx","displacedVtx"),
                                      pfMet             = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
                                      pfMetPhi          = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
                                      rho               = cms.InputTag("hltScoutingPFPacker","rho"),
                                  )

process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)
