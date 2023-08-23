#!/usr/bin/python

# Khalid Jawed
# 03/06/2013
# An attempt to replace MATLAB's parfor by python's parallel processing
# nProcesses is the number of processes that should be run simultaneously
# A list of system commands to be run are written in Commands.txt

# Don't use a command like 'cmd &> output' -> the process won't wait to finish in that case

# Updated: 07/24/2013

# Updated: 08/18/2013
# Usage: python RunBASim3.py [optional_filename]
# filename should have a list of commands one in each line

from multiprocessing import Pool
import subprocess
import time
import os
import sys

#############################################################
# RunBASim function will simply execute a single command
def RunBASim(cmdline):
    BASimProcess = subprocess.Popen(cmdline, shell=True)
    print (cmdline)
    BASimProcess.communicate()
#############################################################

#############################################################
# The main part

nProcesses = 30 #Number of processes to be run simultaneously

last_time = time.time()

if len(sys.argv) == 1:
  CmdFile = 'Commands.txt'
else:
  CmdFile = sys.argv[1]

f = open(CmdFile, 'r')
filecontent = f.readlines()
f.close()
	
p = Pool(processes = nProcesses)

p.map(RunBASim, filecontent)

now = time.time() - last_time
minute = now%60
print ('time taken to run simulations ', minute, ' minutes')
