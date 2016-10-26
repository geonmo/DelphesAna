#!/usr/bin/env python
import os,sys,subprocess

dat="2016-10-26"
lists = [ "TTLL","DY50","DY10to50","ST_tW_top","ST_tW_atop","TT_mass_1665","TT_mass_1695","TT_mass_1755","TT_mass_1785"]
filelists =[]
idxJump=1

for x in lists :
  filelists.append( "filelist_%s.txt"%(x))

for idx, filelist in enumerate(filelists) :
  cmd = "crab out --xrootd -d crab_projects/crab_Production_%s_%d > %s"%(dat, idx+idxJump, filelist)
  print cmd
  #os.system(cmd)
  # check line
  cmd = "wc -l %s"%(filelist)
  print cmd
  fd = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  #print fd.stdout.readlines()[0].split()[0].strip()
  lines = fd.stdout.readlines()
  if ( len(lines) ==0 ) : continue 
  value=lines[0].split()[0].strip()
  if ( int(value) < 2) :
    print "Program is failed. Terminiated!" 
    sys.exit(-1) 
  cmd="sed -i 's#/cms-xrd-global.cern.ch#/cms-xrdr.sdfarm.kr//xrd#g' %s"%(filelist)
  print cmd
  os.system(cmd)

