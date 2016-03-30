#!/bin/env python
# Example PBS cluster job submission in Python

import os
import time
import glob
import re
import subprocess


###############################
RunProxy  = True
DelDir = None #Auto pick up by CMSSW_BASE
DryRun = False
DelExe    = 'testMenu2016'
OutDir  = 'Output/TSGv4Study'
Analysis  = 'METAct'
MenuFile = [
  #"menu/Menu_MuonStudy.txt",
  "menu/Menu_None.txt"
]
Ntuplelist = [
  #"ntuple/r259721_tsgv4Latest.list"
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSG-v4 ~~~~~
  #"ntuple/MC20PU_tsgv4.list",
  #"ntuple/MC30PU_tsgv4.list",
  #"ntuple/MC50PU_tsgv4.list",
  #"ntuple/r258425_tsgv4.list",
  #"ntuple/r258427_tsgv4.list",
  #"ntuple/r258428_tsgv4.list",
  #"ntuple/r258434_tsgv4.list",
  #"ntuple/r258440_tsgv4.list",
  #"ntuple/r258445_tsgv4.list",
  #"ntuple/r258448_tsgv4.list",
  #"ntuple/r259626_tsgv4.list",
  #"ntuple/r259721_tsgv4.list",
  #"ntuple/SingleMuZmu_tsgv4.list",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSGv4-METFix ~~~~~
  "ntuple/r259721_tsgv4METfix.list",
  "ntuple/r259721_tsgv4.list",
  "ntuple/r259721_gomber.list",
]
#GlobalOpt = "--doPlotEff"
#GlobalOpt = "--doPlotRate"
GlobalOpt = "--doPlotTest"
#GlobalOpt = "--doPlotRate --doPrintPU"
Options = {
  "":""
  #"test" : "-n 10"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Muon ER Study ~~~~~
  #"MuER0p8" : "--SetMuonER 0.8",
  #"MuER1p2" : "--SetMuonER 1.25",
  #"MuER2p1" : "--SetMuonER 2.1",
  #"MuER2p2" : "--SetMuonER 2.2",
  #"MuER2p3" : "--SetMuonER 2.3",
  #"MuER2p4" : "--SetMuonER 2.4",
  #"MuER2p5" : "--SetMuonER 2.5",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Muon Turn on Study ~~~~~
  #"r258425" : "--SelectRun 258425",
  #"r258427" : "--SelectRun 258427",
  #"r258428" : "--SelectRun 258428",
  #"r258434" : "--SelectRun 258434",
  #"r258440" : "--SelectRun 258440",
  #"r258445" : "--SelectRun 258445",
  #"r258448" : "--SelectRun 258448",
  #"r259626" : "--SelectRun 259626",
  #"r259721" : "--SelectRun 259721",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MET Cross Check ~~~~~
  #"Ntuple"    : "",
  #"Layer1"    : "--UseUpgradeLyr1",
  #"CaloTower" : "--UseL1CaloTower",

}
#LSFque = '8nm'
LSFque = '8nh'
nBunches = 3963.7 

def BSUB(Analysis, Menu_, Ntuple_, Option):
    global LSFque

    # Open a pipe to the qsub command.
    #output, input = popen2('echo')

    print DelDir, Menu_
    Menu = "%s/%s" % (DelDir, Menu_)
    Ntuple ="%s/%s" % (DelDir, Ntuple_)
    stripMenu = os.path.splitext(os.path.basename(Menu))[0]
    stripNtuple = os.path.splitext(os.path.basename(Ntuple))[0]
    job_name = "%s_%s_%s_%s" % (Analysis, stripMenu, stripNtuple, Option[0])
    out_name = "%s_%s_%s" % (stripMenu, stripNtuple, Option[0])

    cmd = "./%s %s -m %s -l %s -o %s -b %f %s" % (os.path.basename(DelExe), Option[1], Menu, Ntuple,out_name, nBunches, GlobalOpt)
    if DryRun:
       print cmd
       return

    p = subprocess.Popen('bsub', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    output, input = (p.stdout, p.stdin)
    ## -J Job name
    ## -o output -e error
    ## -q Que:
    ### 8nm (8 minutes)
    ### 1nh (1 hour)
    ### 8nh
    ### 1nd (1day)
    ### 2nd
    ### 1nw (1 week)
    ### 2nw
    job_string = """#!/bin/tcsh
    #BSUB -q %s
    #BSUB -J %s
    #BSUB -o %s/%s_stdout
    #BSUB -e %s/%s_stderr
    date
    cp %s ./
    ## Need this to fix the problem with batch node
    setenv LD_LIBRARY_PATH ''
    set TOP="$PWD"
    set CMSSW=%s
    cd $CMSSW/src
    eval `scramv1 runtime -csh`
    cd $TOP
    mkdir menu
    cp $CMSSW/src/L1TriggerDPG/L1Menu/macros/menu/run_lumi.csv menu/
    %s
    ls
    foreach i (`ls results`)
      rfcp results/$i %s
    end
    date""" % (LSFque, job_name, OutDir, job_name, OutDir, job_name, DelExe, os.environ['CMSSW_BASE'], cmd, OutDir)

    #print job_string
    # Send job_string to qsub
    input.write(job_string)
    input.close()

    # Print your job and the response to the screen
    print job_string
    print output.read()

    time.sleep(1)

def my_process():
    global OutDir
    ## Some checking first
    my_CheckFile()

    ## Create the output directory
    OutDir = "%s/%s/%s" % (DelDir , OutDir, Analysis)
    try:
        os.makedirs(OutDir)
    except OSError:
        pass

    for menu in MenuFile:
        for ntuple in Ntuplelist:
            for opt in Options.items():
                #print Analysis, menu, ntuple, opt
                BSUB(Analysis, menu, ntuple, opt)

def my_CheckFile():
    global DelDir
    global DelExe
    ## Check the Delphes Dir
    DelDir = "%s/src/L1TriggerDPG/L1Menu/macros/" % os.environ['CMSSW_BASE']

    ## Check DelFill to be execute
    DelExe = DelDir + "/" + DelExe
    if os.path.isfile(DelExe) and os.access(DelExe, os.X_OK):
        #print "Found DelFill"
        pass
    else:
        DelExe = DelDir + "/" + DelExe.split("/")[-1]
        if os.path.isfile(DelExe) and os.access(DelExe, os.X_OK):
            #print "Found DelFill"
            pass
        else:
            print "Please locate %s" % DelExe
            quit()


if __name__ == "__main__":
    my_process()
