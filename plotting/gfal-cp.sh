#!/bin/bash

# set the base path
BASE_PATH="/cms/store/user/elfontan/ScoutingPFMonitor/monitorSkim_13Feb2023_2022/230505_183721/0000"

# set the destination path
DESTINATION_PATH="/eos/user/e/elfontan/ScoutingParkingPaper/scoutMon_2022F_0000"

# set the range of values for i
for i in {469..999}; do
    # construct the path to the file using the value of i
    FILE_PATH="$BASE_PATH/scoutMonitor_$i.root"

    # construct the gfal-copy command
    GFAL_COMMAND="gfal-copy -t 720000 gsiftp://se01.cmsaf.mit.edu:2811////$FILE_PATH $DESTINATION_PATH"
    
    # print out the command
    echo "Executing command: $GFAL_COMMAND"
    
    # execute the command
    $GFAL_COMMAND    
done
