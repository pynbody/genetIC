#!/usr/bin/env bash

for i in test_*
do
    rm $i/*.tipsy 2>/dev/null
    echo "Running test on $i"
    cd $i
    time ../../IC paramfile.txt > IC_output.txt 2>&1
    if [ $? -ne 0 ]
    then
        cat $i/IC_output.txt
    fi
    cd ..
    echo "Testing output"
    python compare.py $i/*.tipsy $i/reference_output
    if [ $? -ne 0 ]
    then
        echo
        echo "--> TEST FAILED"
        echo
        echo "The IC generator output follows"
        echo
        cat $i/IC_output.txt
        exit
    fi
done

echo "Tests seem OK"
