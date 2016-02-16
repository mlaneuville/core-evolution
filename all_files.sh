#!/bin/bash

for file in ./out/*output.txt; do
    essai=${file%-out*}
    essai=${essai##*/}

    echo "===="
    echo $essai
    python plot.py $essai 10
done
