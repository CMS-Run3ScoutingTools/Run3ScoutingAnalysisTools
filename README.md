# Run3ScoutingAnalysisTools
### Repository for run3 scouting analysis tools

#### Setup
Setup a CMSSW working area and clone the `Run3ScoutingAnalysisTools` repo in the specific branch `ScoutingPaper`:
```
cmsrel CMSSW_13_0_6
cd CMSSW_13_0_6/src
cmsenv
git cms-init
git clone git@github.com:elfontan/Run3ScoutingAnalysisTools.git -b ScoutingPaper
scram b -j 8
```

#### Run on the PFMonitor dataset
Run a basic example on a file of the PFMonitor dataset:
```
voms-proxy-init --voms cms --valid 168:00
cmsRun monitorTree.py
```

#### Run on the scouting dataset
Run a basic example on a file of the Scouting dataset:
```
voms-proxy-init --voms cms --valid 168:00
cmsRun tree.py
```
