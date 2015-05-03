#!/usr/bin/env bash

for i in test_*
do
    rm $i/*.tipsy 2>/dev/null
    echo "Running test on $i"
    cd $i
    ../../IC paramfile.txt > IC_output.txt 2>&1
    cd ..
    echo "Testing output"
    python compare.py $i/*.tipsy $i/reference_output
    if [ $? -ne 0 ]
    then
        echo "--> Test failed."
        exit
    fi
done

echo "Tests seem OK"
