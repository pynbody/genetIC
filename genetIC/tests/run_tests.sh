#!/usr/bin/env bash

function runtest {
  command -v python >/dev/null && PYTHON=python || PYTHON=python3
  rm $1/*.tipsy 2>/dev/null
  echo -n "Running test on $1   "
  head -1  $1/paramfile.txt
  cd $1 || exit
  IC=${IC:-../../genetIC}
  command="$IC paramfile.txt > IC_output.txt 2>&1"
  eval "$command"
  if [[ $? -ne 0 && "$1" != *error* ]]
  then
      echo "--> TEST ERRORED"
      cat IC_output.txt
      exit 1
  fi
  cd ..
  echo "Testing output"
  $PYTHON ../../tools/compare.py $1/
  if [ $? -ne 0 ]
  then
      echo
      echo "--> TEST FAILED"
      echo
      echo "The IC generator output follows"
      echo
      cat $1/IC_output.txt
      exit 1
  else
      echo
  fi
}


if [ "$#" -eq 0 ]; then
  for i in test_*
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
