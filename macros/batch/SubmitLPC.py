#!/bin/env python
# Example PBS cluster job submission in Python

import os
import time
import glob
import re
import subprocess
import tarfile
import time

###############################
DelDir = None #Auto pick up by CMSSW_BASE
tempdir =None
ProjectName = None
DryRun = False
DelExe    = 'testMenu2016'
OutDir = '/store/user/benwu/L1MenuStage2/TSGv4'
Analysis  = 'test'
MenuFile = [
  #"menu/Menu_MuonStudy.txt",
  #"menu/Menu_None.txt"
  #"menu/Menu_259721_TSGv4_Riccardo.txt"
  # "menu/Menu_259721_TSGv4_Prescales.txt",
  # "menu/Menu_259721_TSGv3_FixPre_EG.txt",
  "menu/Menu_259721_TSGv4_FixPre.txt",
  # "menu/Menu_ETMStudy.txt",
]
Ntuplelist = [
  #"ntuple/r259721_tsgv4Latest.list"
  # "ntuple/r259721_tsgv3.list",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSG-v4 ~~~~~
    #"ntuple/MC10PU_tsgv4.list",
    #"ntuple/MC20PU_tsgv4.list",
    #"ntuple/MC30PU_tsgv4.list",
    #"ntuple/MC40PU_tsgv4.list",
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
    "ntuple/r259721_tsgv4ZB1.list",
    # "ntuple/SingleMuZmu_tsgv4.list",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSGv4-METFix ~~~~~
  #"ntuple/r259721_tsgv4METfix.list",
  #"ntuple/r259721_tsgv4.list",
  #"ntuple/r259721_gomber.list",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Integration v19 ~~~~~
  #"ntuple/r258440_itgv19Layer1.list",
  #"ntuple/r258440_itgv19.list",
  #"ntuple/r259626_itgv19Layer1.list",
  #"ntuple/r259626_itgv19.list",
  #"ntuple/r259721_itgv19Layer1.list",
  # "ntuple/r259721_itgv19.list",
  # "ntuple/SingleMuZmu_itgv19Layer1.list",
  # "ntuple/SingleMuZmu_itgv19.list",
]
GlobalOpt =  " "
GlobalOpt += " --doPlotEff"
GlobalOpt += " --doPlotRate"
GlobalOpt += " --doPlotTest"
#GlobalOpt += " --doPlotRate --doPrintPU"
# GlobalOpt = " --doPrintPU"
Options = {
  #None:""
  "test" : "-n 10000"
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
  #"r259721" : "--SelectRun 259721",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MET Cross Check ~~~~~
  #"Default"    : "",
  #"Bitwise"    : "--UseUpgradeLyr1",
  #"CaloTower" : "--UseL1CaloTower",

}
nBunches = 3963.7

def CondorSub(Analysis, Menu, Ntuple, Option):
    npro =[ "%s/ntuple.tgz" % tempdir, "%s/menu.tgz" % tempdir]
    npro.append( DelExe )

    ## Update RunHT.csh with DelDir and pileups
    RunHTFile = tempdir + "/" + "RunExe.csh"
    with open(RunHTFile, "wt") as outfile:
        for line in open("%s/RunExe.csh" % os.path.dirname(os.path.realpath(__file__)), "r"):
            #line = line.replace("DELDIR", os.environ['PWD'])
            line = line.replace("DELDIR", os.environ['CMSSW_BASE'])
            line = line.replace("DELEXE", DelExe.split('/')[-1])
            line = line.replace("OUTDIR", OutDir)
            outfile.write(line)


    stripMenu = os.path.splitext(os.path.basename(Menu))[0]
    stripNtuple = os.path.splitext(os.path.basename(Ntuple))[0]
    job_name = "%s_%s_%s_%s" % (Analysis, stripMenu, stripNtuple, Option[0])
    if len(MenuFile) == 1:
      out_name = "_".join(filter(None, (stripNtuple, Option[0])))
    else:
      out_name = "_".join(filter(None, (stripMenu, stripNtuple, Option[0])))

    arg = " %s -m %s -l %s -o %s -b %f %s" % (Option[1], Menu, Ntuple,out_name, nBunches, GlobalOpt)
    cmd = "./%s %s" % (os.path.basename(DelExe),  arg )
    if DryRun:
       print cmd
       return

    tranferfiles = ", ".join(npro)
    arg = "\nArguments = %s \n Queue\n" % arg
    ## Prepare the condor file
    condorfile = tempdir + "/" + "condor_" + ProjectName
    with open(condorfile, "wt") as outfile:
        for line in open("%s/condor_template" % os.path.dirname(os.path.realpath(__file__)), "r"):
            line = line.replace("EXECUTABLE", os.path.abspath(RunHTFile))
            #line = line.replace("DELDIR", os.environ['CMSSW_BASE'])
            line = line.replace("TARFILES", tranferfiles)
            line = line.replace("TEMPDIR", tempdir)
            line = line.replace("PROJECTNAME", ProjectName)
            line = line.replace("ARGUMENTS", arg)
            outfile.write(line)

    Condor_Sub(condorfile)
    time.sleep(1)

def my_process():
    global OutDir
    global tempdir
    global ProjectName
    ## Some checking first
    my_CheckFile()

    ## Create the output directory
    try:
        os.makedirs(OutDir)
    except OSError:
        pass

    ProjectName = time.strftime('%b%d') + Analysis
    tempdir = '/tmp/' + os.getlogin() + "/" + ProjectName +  "/"
    try:
        os.makedirs(tempdir)
    except OSError:
        pass
    ## Create the output directory
    OutDir = OutDir +  "/" + ProjectName + "/"
    try:
        os.makedirs(OutDir)
    except OSError:
        pass

    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    os.chdir("../")
    subprocess.call("tar -czf %s/ntuple.tgz ntuple/" % tempdir, shell=True)
    subprocess.call("tar -czf %s/menu.tgz menu/" % tempdir, shell=True)
    os.chdir(curdir)

    for menu in MenuFile:
        for ntuple in Ntuplelist:
            for opt in Options.items():
                #print Analysis, menu, ntuple, opt
                CondorSub(Analysis, menu, ntuple, opt)

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


def Condor_Sub(condor_file):
    ## Since we run with xrootd, setup proxy
    #hasproxy = False
    #proxyfile = ''
    #if os.environ.has_key('X509_USER_PROXY'):
        #proxyfile = os.path.abspath(os.environ['X509_USER_PROXY'])
    #else:
        #hasproxy = False
    #if not hasproxy or not os.path.exists(proxyfile) or (time.time() - os.stat(proxyfile).st_ctime) / 60/24 > 1:
        #print "Proxy file is at least one day old. Requesting new proxy"
        #os.system("voms-proxy-init -valid 168:00 -voms cms")

    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.path.dirname(condor_file))
    print "To submit condor with " + condor_file
    os.system("condor_submit " + condor_file)
    os.chdir(curdir)



if __name__ == "__main__":
    my_process()
