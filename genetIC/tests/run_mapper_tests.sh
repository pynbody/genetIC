#!/usr/bin/env bash

function runtest {
    rm $i/map_out.txt 2>/dev/null
    echo "Running test on $i"

    cd $i
    IC_mapper=${IC_mapper:-../../genetIC_mapper}
    # if paramfile_a.txt exists, use it, otherwise use paramfile.txt
    PARAM_A=paramfile_a.txt
    if [ ! -f $PARAM_A ]; then
        PARAM_A=paramfile.txt
    fi

    head -1 $PARAM_A

    # if paramfile_b.txt exists, use it, otherwise use paramfile.txt
    PARAM_B=paramfile_b.txt
    if [ ! -f $PARAM_B ]; then
        PARAM_B=paramfile.txt
    fi
    #if ID_a.txt exists, use it; otherwise, assume we just want to output flags at end of processing paramfile_a
    ID_A=ID_a.txt
    if [ ! -f $ID_A ]; then
      $IC_mapper $PARAM_A output.txt > IC_output.txt 2>&1
    else
      $IC_mapper $PARAM_A $PARAM_B $ID_A output.txt > IC_output.txt 2>&1
    fi
    if [ $? -ne 0 ]
    then
        echo "TEST FAILED"
        cat IC_output.txt
        exit
    fi
    cd ..
    echo "Testing output"
    DIFF_RESULT=$(diff $i/reference.txt $i/output.txt)
    if [ "$DIFF_RESULT" != "" ]
    then
        echo $DIFF_RESULT
        echo
        echo "--> TEST FAILED"
        echo
        echo "The IC generator output follows"
        echo
        cat $i/IC_output.txt
        exit 1
    fi
}


if [ "$#" -eq 0 ]; then
  for i in mapper_test_*
  do
    runtest $i
  done
else
  for i in $@
  do
    runtest $i
  done
fi


echo "Tests seem OK"
