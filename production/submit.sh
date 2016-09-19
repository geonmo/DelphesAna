
#!/usr/bin/env python
import os

datasetList = open("dataset_miniaod2.txt").readlines()

for dataset in datasetList :

  dataset=  dataset.strip()
  if dataset=="" : continue
  #print dataset
  cmd="crab submit -c crabConfig.py Data.inputDataset=\'%s\'"%(dataset)
  print cmd
  os.system(cmd)

