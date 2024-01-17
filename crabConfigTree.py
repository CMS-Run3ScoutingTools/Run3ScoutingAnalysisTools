from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'mmSkim_data_2023D'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'tree.py'

config.Data.inputDataset = '/ScoutingPFRun3/Run2023D-v1/RAW'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'mmSkim_data_2023D'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json' # Run 3 2022 
config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json' # Run 3 2023

# Where the output files will be transmitted to
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_CH_CERNBOX'
