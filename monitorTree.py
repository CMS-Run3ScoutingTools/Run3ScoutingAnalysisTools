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
                                "/store/data/Run2024F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/382/262/00000/97ff508a-3a26-4cbe-a205-e500efb89a26.root", 
                                "/store/data/Run2024F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/382/262/00000/3d3eeb8b-98b5-411f-9a4b-e1791f7da3dd.root",  
                                "/store/data/Run2024F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/382/262/00000/f646d5cf-a8f6-4702-adb2-4c552ebaacbb.root",
                                "/store/data/Run2024F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/382/262/00000/c5f8b10d-de32-4405-958a-a2e2b6a1de29.root",  
                                ),                                                                 
                            lumisToProcess = cms.untracked.VLuminosityBlockRange('382262:1-382262:273')
)

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "rawDataCollector", "", "LHC" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scoutMonitor.root")
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v1', '') # Run 3 2024

L1Seeds = ["L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu8_SQ","L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0_Upt6_SQ_er2p0","L1_DoubleMu0_Upt7_SQ_er2p0","L1_DoubleMu0_Upt8_SQ_er2p0", "L1_DoubleMu0_Upt15_Upt7", "L1_DoubleMu0_Upt6_IP_Min1_Upt4"]

L1MonitorSeeds = ["L1_SingleMu22","L1_SingleMu25","L1_HTT200er","L1_HTT255er","L1_HTT280er","L1_HTT320er","L1_HTT360er","L1_HTT400er","L1_HTT450er","L1_ETT2000", "L1_HTT280er_QuadJet_70_55_40_35_er2p5","L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3","L1_SingleEG34er2p5", "L1_SingleEG36er2p5","L1_SingleEG38er2p5","L1_SingleEG40er2p5","L1_SingleJet160er2p5","L1_SingleJet180","L1_SingleJet200","L1_SingleTau120er2p1","L1_SingleTau130er2p1","L1_SingleEG42er2p5","L1_SingleEG45er2p5","L1_SingleEG60", "L1_DoubleEG_LooseIso18_LooseIso12_er1p5","L1_DoubleEG_LooseIso20_LooseIso12_er1p5","L1_DoubleEG_LooseIso22_LooseIso12_er1p5"]

process.scoutingTree = cms.EDAnalyzer('ScoutingTreeMakerRun3Monitor',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Seeds),
                                      l1MonitorSeeds    = cms.vstring(L1MonitorSeeds),
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
