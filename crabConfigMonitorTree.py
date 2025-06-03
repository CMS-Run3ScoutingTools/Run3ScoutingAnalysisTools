from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'monitorSkim_muVtx_2024I'
config.JobType.pluginName = 'Analysis'

# Name of the CMSSW configuration file
# ------------------------------------
config.JobType.psetName = 'monitorTree.py'

config.Data.inputDataset = '/ScoutingPFMonitor/Run2024I-PromptReco-v1/MINIAOD'
#config.Data.secondaryInputDataset = '/ScoutingPFMonitor/Run2024E-v1/RAW'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
# This string is used to construct the output dataset name
# --------------------------------------------------------
config.Data.outputDatasetTag = 'monitorSkim_muVtx_2024I'

# These values only make sense for processing data: Select input data based on a lumi mask
# ----------------------------------------------------------------------------------------
# Golden JSON 2023: https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json
#config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json'
# Golden JSON 2024: https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json 
config.Data.lumiMask = 'Cert_Collisions2024_378981_386951_Golden.json'

#config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T3_CH_CERNBOX'
