import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",                            
                            fileNames = cms.untracked.vstring(
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/361/971/00000/6c45ebf6-052c-456f-b085-1def7f38830c.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/361/971/00000/2629b39a-a626-4774-be78-e1e6783b755c.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/361/971/00000/6c45ebf6-052c-456f-b085-1def7f38830c.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/361/971/00000/71a94f03-6acf-46fb-af77-5e29709199c3.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/362/167/00000/22266449-6097-4832-9eef-92c0c5e0f8dc.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/362/167/00000/27f02c1d-04d6-4af3-bf91-9fb5ba828997.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/362/167/00000/5e0ec603-8813-4227-9dc4-95b46e987588.root"
                                #
                                "/store/data/Run2023C/ScoutingPFMonitor/MINIAOD/PromptReco-v4/000/368/823/00000/cccd7032-a259-4e60-a21b-1b88c041c769.root",
                            ),
                            secondaryFileNames=cms.untracked.vstring(
                                #"/store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/361/971/00000/ba851896-96bc-4c21-b45a-dcb8ac57e272.root",
                                #"/store/data/Run2022F/ScoutingPFMonitor/RAW/v1/000/362/167/00000/bce2e687-410b-4a10-bbed-60c82d5ddf40.root"),
                                #
                                "/store/data/Run2023C/ScoutingPFMonitor/RAW/v1/000/368/823/00000/388bf0bd-900f-410a-864c-b77b243d11ff.root",
                                "/store/data/Run2023C/ScoutingPFMonitor/RAW/v1/000/368/823/00000/9da4e47e-8e12-4a7e-8c71-88458eef0563.root",
                                "/store/data/Run2023C/ScoutingPFMonitor/RAW/v1/000/368/823/00000/66c984ff-6f83-4c2e-9c77-34b6714fe479.root"
                            ),
                            lumisToProcess = cms.untracked.VLuminosityBlockRange('368823:3-368823:35')
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('361971:2303-361971:2327','362167:87-362167:96')
)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
#process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.gtStage2Digis.InputLabel = cms.InputTag( "rawDataCollector", "", "LHC" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scoutMonitor.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '') # Run 3 2022
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v2', '') # Run 3 2023

L1Seeds = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2", "L1_DoubleMu8_SQ"]
#L1Seeds = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2", "L1_SingleLooseIsoEG28er2p1", "L1_SingleLooseIsoEG28er1p5", "L1_SingleLooseIsoEG30er1p5", "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG32er2p1", "L1_DoubleEG_LooseIso16_LooseIso12_er1p5", "L1_DoubleEG_LooseIso18_LooseIso12_er1p5", "L1_DoubleEG_LooseIso20_LooseIso12_er1p5", "L1_DoubleEG_LooseIso22_LooseIso12_er1p5", "L1_SingleJet180", "L1_SingleJet200", "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", "L1_HTT280er", "L1_HTT320er", "L1_HTT360er", "L1_ETT2000"]

L1MonitorSeeds = ["L1_SingleMu22","L1_SingleMu25","L1_HTT200er","L1_HTT255er","L1_HTT280er","L1_HTT320er","L1_HTT360er","L1_HTT400er","L1_HTT450er","L1_ETT2000", "L1_HTT280er_QuadJet_70_55_40_35_er2p5","L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3","L1_SingleEG34er2p5", "L1_SingleEG36er2p5","L1_SingleEG38er2p5","L1_SingleEG40er2p5","L1_SingleJet160er2p5","L1_SingleJet180","L1_SingleJet200","L1_SingleTau120er2p1","L1_SingleTau130er2p1","L1_SingleEG42er2p5","L1_SingleEG45er2p5","L1_SingleEG60", "L1_DoubleEG_LooseIso18_LooseIso12_er1p5","L1_DoubleEG_LooseIso20_LooseIso12_er1p5","L1_DoubleEG_LooseIso22_LooseIso12_er1p5"]

#L1MonitorSeeds = ["L1_HTT200er","L1_HTT255er","L1_HTT280er","L1_HTT320er","L1_HTT360er","L1_HTT400er","L1_HTT450er","L1_ETT2000","L1_SingleJet180","L1_SingleJet200","L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5","L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5","L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5","L1_SingleLooseIsoEG28er2p1","L1_SingleLooseIsoEG28er1p5","L1_SingleLooseIsoEG30er1p5","L1_SingleIsoEG28er2p1","L1_SingleIsoEG30er2p1","L1_SingleIsoEG32er2p1","L1_DoubleEG_LooseIso16_LooseIso12_er1p5","L1_DoubleEG_LooseIso18_LooseIso12_er1p5","L1_DoubleEG_LooseIso20_LooseIso12_er1p5","L1_DoubleEG_LooseIso22_LooseIso12_er1p5"]

process.scoutingTree = cms.EDAnalyzer('ScoutingTreeMakerRun3Monitor',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Seeds),
                                      l1MonitorSeeds    = cms.vstring(L1MonitorSeeds),
                                      muons             = cms.InputTag("hltScoutingMuonPacker","","HLT"),
                                      offlineMuons      = cms.untracked.InputTag("slimmedMuons"),
                                      electrons         = cms.InputTag("hltScoutingEgammaPacker"),
                                      offlinePhotons    = cms.untracked.InputTag("slimmedPhotons"),
                                      photons           = cms.InputTag("hltScoutingEgammaPacker"),
                                      pfcands           = cms.InputTag("hltScoutingPFPacker"),
                                      pfjets            = cms.InputTag("hltScoutingPFPacker"),
                                      tracks            = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices   = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
                                      pfMet             = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
                                      pfMetPhi          = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
                                      rho               = cms.InputTag("hltScoutingPFPacker","rho"),
                                  )

process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)
