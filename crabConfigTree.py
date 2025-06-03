from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'mmSkim_data_2024I_bothNovtxVtxMu_v1'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
# ------------------------------------
config.JobType.psetName = 'tree.py'

config.Data.inputDataset = '/ScoutingPFRun3/Run2024I-v1/HLTSCOUT'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
# This string is used to construct the output dataset name
# --------------------------------------------------------
config.Data.outputDatasetTag = 'mmSkim_data_2024I_bothNovtxVtxMu_v1'

# These values only make sense for processing data: select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2024_378981_386951_Golden.json' # Run 3 2024
config.Data.lumiMask = 'Cert_Collisions2024_F-I_Golden.json'

# Where the output files will be transmitted to
# ---------------------------------------------
config.Site.storageSite = 'T3_CH_CERNBOX'
