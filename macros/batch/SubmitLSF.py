#!/bin/env python
# Example PBS cluster job submission in Python

import os
import time
import glob
import re
import subprocess


###############################
DelDir = None #Auto pick up by CMSSW_BASE
DryRun = False
DelExe    = 'testMenu2016'
# OutDir = 'Output/ITGv27/'
# OutDir  = 'Output/ITGv19Study/'
# OutDir  = 'Output/TSGv4Vsv3/'
OutDir  = 'Output/ITGv59'
# Analysis  = 'MuonPU2'
# Analysis  = 'MenuPU3'
# Analysis  = 'PFMET'
Analysis  = 'doubleMu2'
MenuFile = [
  # "menu/Menu_MuonStudy.txt",
  # "menu/Menu_MuonStudy.txt",
  # "menu/Menu_None.txt"
  # "menu/Menu_259721_TSGv4_Prescales.txt",
  # "menu/Menu_259721_TSGv3_FixPre_EG.txt",
  "menu/Lumi5E33_TSGv5_Prescales.txt",
  # "menu/Menu_259721_TSGv5_Prescales.txt",
  # "menu/Menu_ETMStudy.txt",
]
Ntuplelist = [
  #"ntuple/r259721_tsgv4Latest.list"
  # "ntuple/r259721_tsgv3.list",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSG-v4 ~~~~~
    # "ntuple/MC10PU_tsgv4.list",
    # "ntuple/MC20PU_tsgv4.list",
    # "ntuple/MC30PU_tsgv4.list",
    # "ntuple/MC40PU_tsgv4.list",
    # "ntuple/MC50PU_tsgv4.list",
    # "ntuple/r258425_tsgv4.list",
    # "ntuple/r258427_tsgv4.list",
    # "ntuple/r258428_tsgv4.list",
    # "ntuple/r258434_tsgv4.list",
    # "ntuple/r258440_tsgv4.list",
    # "ntuple/r258445_tsgv4.list",
    # "ntuple/r258448_tsgv4.list",
    # "ntuple/r259626_tsgv4.list",
    # "ntuple/r259721_tsgv4.list",
    # "ntuple/SingleMuZmu_tsgv4.list",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSGv4-METFix ~~~~~
    #"ntuple/r259721_tsgv4METfix.list",
    #"ntuple/r259721_tsgv4.list",
    #"ntuple/r259721_gomber.list",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Integration v19 ~~~~~
    # "ntuple/r258440_itgv19Layer1.list",
    # "ntuple/r258440_itgv19.list",
    # "ntuple/r259626_itgv19Layer1.list",
    # "ntuple/r259626_itgv19.list",
    # "ntuple/r259721_itgv19Layer1.list",
    # "ntuple/r259721_itgv19.list",
    # "ntuple/SingleMuZmu_itgv19Layer1.list",
    # "ntuple/SingleMuZmu_itgv19.list",
    # "ntuple/r259721_itgv19Layer1_Reco.list",
    # "ntuple/r259721_itgv19_Reco.list",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSGv4p1 ~~~~~
    # "ntuple/MC20PU_tsgv4p1.list",
    # "ntuple/MC30PU_tsgv4p1.list",
    # "ntuple/r258440_tsgv4p1.list",
    # "ntuple/r259626_tsgv4p1.list",
    # "ntuple/r259721_tsgv4p1.list",
    # "ntuple/SingleMuZmu_tsgv4p1.list",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Integration v27 ~~~~~
    # "ntuple/MC20PU_itgv27.list",
    # "ntuple/MC30PU_itgv27.list",
    # "ntuple/MC40PU_itgv27.list",
    # "ntuple/r258440_itgv27.list",
    # "ntuple/r259626_itgv27.list",
    # "ntuple/r259721_itgv27.list",
    # "ntuple/SingleMuZmu_itgv27.list",
    # "ntuple/r259721_itgv35.list"
    # "ntuple/r259721_itgv37p2_Opt5.list"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Com ~~~~~
    # "ntuple/r271306_unpack.list",
    # "ntuple/r271306_TSGv5p1.list"
    # "ntuple/r258440_itgv42p1.list",
    # "ntuple/r259626_itgv42p1.list",
    # "ntuple/r259721_itgv42p1.list",

    # "ntuple/r271071_itgv42p1.list",
    # "ntuple/r271074_itgv42p1.list",
    # "ntuple/r271075_itgv42p1.list",
    # "ntuple/r271084_itgv42p1.list",
    # "ntuple/r271306_itgv42p1.list",
    # "ntuple/r271336_itgv42p1.list",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Itg-v42.1 ~~~~~
    # "ntuple/r271336_itgv42p1.list",
    # "ntuple/r271336_unpack.list",

    # "ntuple/r272011_itgv42p1.list",
    # "ntuple/r272011_unpack.list",
    # "ntuple/r272022_itgv42p1.list",
    # "ntuple/r272022_unpack.list",

    # "ntuple/r258440_itgv46p0.list",
    # "ntuple/r259626_itgv46p0.list",
    # "ntuple/r259721_itgv46p0.list",
    # "ntuple/r272011_itgv46p0.list",
    # "ntuple/r272022_itgv46p0.list",


    # "ntuple/r272011_itgv48p2.list",
    # "ntuple/r272022_itgv48p2.list",
    # "ntuple/r272784_itgv48p2.list",
    # "ntuple/r272798_itgv48p2.list",
    # "ntuple/r272812_itgv48p2.list",
    # "ntuple/r272818_itgv48p2.list",
    # "ntuple/r272828_itgv48p2.list",
    # "ntuple/r272011_itgv48.list",
    # "ntuple/r272022_itgv48.list",
    # "ntuple/r272784_itgv48.list",
    # "ntuple/r272798_itgv48.list",
    # "ntuple/r272812_itgv48.list",
    # "ntuple/r272818_itgv48.list",
    # "ntuple/r272011_unpack.list",
    # "ntuple/r272022_unpack.list",
    # "ntuple/r272784_unpack.list",
    # "ntuple/r272798_unpack.list",
    # "ntuple/r272812_unpack.list",
    # "ntuple/r272818_unpack.list",

    "ntuple/r273725_itgv59.list",
    "ntuple/r273728_itgv59.list",
    "ntuple/r273730_itgv59.list",
    # "ntuple/r273725_itgv58.list",
    # "ntuple/r273728_itgv58.list",
    # "ntuple/r272828_unpack.list",
    # "ntuple/r272022_unpack.list",
    # "ntuple/r272011_unpack.list",
    # "ntuple/r272784_unpack.list",
    # "ntuple/r272812_unpack.list",
    # "ntuple/r272818_unpack.list",
    # "ntuple/r272798_unpack.list",
]
GlobalOpt = ""
# GlobalOpt += " --doPlotEff"
GlobalOpt += " --doPlotRate"
# GlobalOpt += " --doPlotTest"
GlobalOpt += " --SetNoPrescale"
GlobalOpt += " --doPrintPU"
# GlobalOpt += " --UsePFMETNoMuon"
# GlobalOpt += " --doPlotRate --doPrintPU"
Options = {
  # None:"",
  #"test" : "-n 10"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Muon ER Study ~~~~~
  #"MuER0p8" : "--SetMuonER 0.8",
  #"MuER1p2" : "--SetMuonER 1.25",
  #"MuER1p5" : "--SetMuonER 1.5",
  #"MuER1p8" : "--SetMuonER 1.8",
  #"MuER2p0" : "--SetMuonER 2.0",
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MET Cross Check ~~~~~
  #"r259721" : "--SelectRun 259721",
  "EmuC++"    : "",
  "EmuGT"    : "--UseuGTDecision",
  # #"Bitwise"    : "--UseUpgradeLyr1",
  # "Float" : " --SumJetET  30 --UseL1CaloTower",
  "UnpackuGT"    : " --UseuGTDecision --UseUnpackTree",
  "UnpackC++"    : " --UseUnpackTree",

}
# LSFque = '8nm'
LSFque = '8nh'
# LSFque = '1nd'
# nBunches = 6840 ## From 12PU scaled to 30PU
nBunches = 1165
# nBunches = 8
# nBunches = 3963.7

def BSUB(Analysis, Menu_, Ntuple_, Option):
    global LSFque

    # Open a pipe to the qsub command.
    #output, input = popen2('echo')

    Menu = "%s/%s" % (DelDir, Menu_)
    Ntuple ="%s/%s" % (DelDir, Ntuple_)
    stripMenu = os.path.splitext(os.path.basename(Menu))[0]
    stripNtuple = os.path.splitext(os.path.basename(Ntuple))[0]
    job_name = "%s_%s_%s_%s" % (Analysis, stripMenu, stripNtuple, Option[0])
    if len(MenuFile) == 1:
      out_name = "_".join(filter(None, (stripNtuple, Option[0])))
    else:
      out_name = "_".join(filter(None, (stripMenu, stripNtuple, Option[0])))

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
