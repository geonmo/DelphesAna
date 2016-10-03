#!/usr/bin/env python
import os,sys

fileName = sys.argv[1]
datasetList = open(fileName).readlines()
print "File is opened."
for dataset in datasetList :

  dataset=  dataset.strip()
  if dataset=="" : continue
  #print dataset
  cmd="crab submit -c crabConfig.py Data.inputDataset=\'%s\'"%(dataset)
  print cmd
  os.system(cmd)

