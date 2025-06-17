# Run3ScoutingAnalysisTools
### Repository for Run 3 scouting analysis tools

#### Setup
Setup a CMSSW working area and clone the `Run3ScoutingAnalysisTools` repo in the specific branch `Run3`:
```
cmsrel CMSSW_15_0_6
cd CMSSW_15_0_6/src
cmsenv
git cms-init
git clone git@github.com:CMS-Run3ScoutingTools/Run3ScoutingAnalysisTools.git -b Run3
scram b -j 8
```

#### Run on the PFMonitor dataset
Run a basic example on a file of the PFMonitor dataset:
```
voms-proxy-init --voms cms --valid 168:00
cmsRun monitorTree.py
```

If no errors are observed, you can proceed to the crab submission taking care of updating properly the configuration file chenging the dataset name and output name and checking the certification json file:
``` 
crab submit crabConfigMonitorTree.py
``` 
#### Run on the scouting dataset
Run a basic example on a file of the Scouting dataset:
```
voms-proxy-init --voms cms --valid 168:00
cmsRun tree.py
```
If no errors are observed, you can proceed to the crab submission taking care of updating properly the configuration file chenging the dataset name and output name and checking the certification json file:
``` 
crab submit crabConfigTree.py
``` 

#### Check production status and luminosity
The status of the jobs can be monitored and failed jobs can be resubmitted:

``` 
crab status CRABDIR
crab resubmit CRABDIR
```
After having run successfully all jobs, the integrated luminosity corresponding to the processed data can be checked using [brilcalc](https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart) in the following way:
``` 
crab report CRABDIR
source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
brilcalc --version
brilcalc lumi -i CRABDIR/results/processedLumis.json
``` 
