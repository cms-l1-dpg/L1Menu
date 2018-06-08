#!/usr/bin/env python
# encoding: utf-8

# File        : CopyToLPC.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2016 Apr 11
#
# Description :


import subprocess
import subprocess
import os
from multiprocessing import Pool


#filename = "CERN_file_list/run_316380.list"
filename = "CERN_file_list/run_317087.list"

def fwork(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    p = Pool(8)
    cmdlist = []
    f = open(os.path.basename(filename), "w")
    #f = open(filename, "w")
    with open("%s" % filename) as file:
        for line_ in file.readlines():
            line =line_.strip()
            newloc = line.replace("eoscms.cern.ch", "cmseos.fnal.gov")
            newloc = newloc.replace("dpg_trigger", "lpctrig")
            newloc = newloc.replace("/cms/", "/uscms/")
            print newloc
            f.write(newloc+"\n")
            cmdlist.append("xrdcp %s %s" % ( line, newloc ))
    f.close()
    print p.map(fwork, cmdlist)

