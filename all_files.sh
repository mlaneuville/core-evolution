#!/bin/bash
# The command line argument 'folder' is optional. If set, the script produces the output
# of all files in that subfolder.

FOLDER=$1
echo $FOLDER

for file in ./out/${FOLDER}/*output.txt; do
    essai=${file%-out*}
    essai=${essai##*/}

    if [ $FOLDER ]; then
        python plot.py -s $FOLDER $essai -n5
    else
        python plot.py $essai -n5
    fi
    
    echo "===="
    echo $essai
done
