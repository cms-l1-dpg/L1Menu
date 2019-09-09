#!/bin/env python
# Example PBS cluster job submission in Python

import os
import time
import glob
import copy
import re
import subprocess
import tarfile
import time
import shutil
from collections import defaultdict

###############################
DelDir = None #Auto pick up by CMSSW_BASE
#tempdir = '/uscmst1b_scratch/lpc1/lpctrig/benwu/CondorTemp'

#bx_id =11

tempdir = '/uscms_data/d3/huiwang/condor_temp'
ProjectName = "Menu2017"
DryRun = False
splitline = 5
DelExe    = 'testMenu2016'
OutDir = '/store/user/huiwang/L1Menu2017'
Analysis  = 'fill_7118_and_7131_nanoDST_Prescale_2018_v2_1_0_Col_2.0_HATS'
#Analysis  = 'fill_7358_nanoDST_Prescale_2018_v2_1_0_Col_2.0_BX%d' %bx_id
#Analysis  = 'Marco_Official_collision_V3_menu_prescale_Col_2.0_run_317648_and_317649'
MenuFile = [
  #"menu/Prescale_Sets_RUN_306091_col_1.6.txt"
  "menu/Prescale_2018_v2_1_0_Col_2.0_HATS.txt"
  #"menu/Marco_Official_collision_V3_menu_prescale_Col_2.0.txt"
]
Ntuplelist = [
]
Ntupledict = {
    # "ntuple/Trains_v95p12p2.list" : " --SelectBX \\\"[[714, 761], [1875, 1922]]\\\"  -u menu/TrainPLTZ.csv ",
    # "ntuple/fill_6356_6360.list" : " -u menu/runlumi_fill_6358_and_more.csv ",
    # "ntuple/run_316216_nanodst.list" : " -u menu/runlumi_fill_6358_and_more.csv ",
     "ntuple/fill_7118_and_7131_nanoDST_LPC.list" : " -u menu/runlumi_fill_6358_and_more.csv ",
    # "ntuple/fill_7358_nanoDST.list" : " -u menu/PU_12b/PU_bx_%d.csv" % bx_id,
}
GlobalOpt =  " "
#GlobalOpt += " --SelectRun 299380"
#GlobalOpt += " --SetNoPrescale"
#GlobalOpt += " --IgnorePrescale"
#GlobalOpt += " --doScanLS --SelectLS '[151,200]' "
#GlobalOpt += " --doPlotRate --doPrintPU --UseUnpackTree --Select_BX_in_12b %d" % bx_id
GlobalOpt += " --doPlotRate --doPrintPU --UseUnpackTree"
#GlobalOpt += " --SelectCol 1.8E34 "
Options = {
  "Default"    : "",
  # None:" "
    #"MenuPU2t"    : " --doPrintPU --IgnorePrescale --SelectCol 2.01 ",
    # "MenuPU2p2"  : " --doPrintPU --IgnorePrescale --SelectCol 2.20 ",
    # "MenuPU2"    : " --doPrintPU --IgnorePrescale --SelectCol 2.00 ",
    #"MenuPU1p8t"  : " --doPrintPU --IgnorePrescale --SelectCol 1.81 ",
    # "MenuPU1p8"  : " --doPrintPU --IgnorePrescale --SelectCol 1.80 ",
    # "MenuPU1p81" : " --doPrintPU --IgnorePrescale --SelectCol 1.81 ",
    # "EmuPU"      : " --doPrintPU --SetNoPrescale  --SelectCol 2.00 ",
    # "UnpackPU"   : " --doPrintPU --UseUnpackTree --SetNoPrescale  --SelectCol 2.00 ",
    #"Scan44PU" : " --doScanLS --SelectLS \"[151, 200]\" ",    #44PU
    # "Scan38PU" : " --doScanLS --SelectLS \\\"[365,395]\\\" ", #38PU
    # "Scan45PU" : " --doScanLS --SelectLS \\\"[209,233]\\\" ", #45PU
    # "Scan47PU" : " --doScanLS --SelectLS \\\"[175,193]\\\" ", #47PU
    # # "Scan50PU" : " --doScanLS --SelectLS \\\"[126,144]\\\" ", #50PU
    # "Scan33p9PU" : " --doScanLS --SelectLS \\\"[520,530]\\\" ", #33.9PU

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MET Cross Check ~~~~~
  #"Bitwise"    : "--UseUpgradeLyr1",
  #"CaloTower" : "--UseL1CaloTower",

}
nBunches = 2544
#nBunches = 24
#nBunches = 1866


def CondorSub(Analysis, Menu, Ntuple, Option, cnt):
    npro =[ "%s/ntuple.tgz" % tempdir, "%s/menu.tgz" % tempdir]
    npro.append( DelExe )

    ## Update RunHT.csh with DelDir and pileups
    RunHTFile = tempdir + "/" + "RunExe.csh"
    with open(RunHTFile, "wt") as outfile:
        for line in open("%s/RunExe.csh" % os.path.dirname(os.path.realpath(__file__)), "r"):
            #line = line.replace("DELDIR", os.environ['PWD'])
            #line = line.replace("DELDIR", os.environ['CMSSW_BASE'])
            line = line.replace("DELSCR", os.environ['SCRAM_ARCH'])
            line = line.replace("DELDIR", os.environ['CMSSW_VERSION'])
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
    arg = "\nArguments = %s \n Queue %d\n" % (arg, cnt)
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
    tempdir += '/' + os.getlogin() + "/" + ProjectName +  "/"
    try:
        os.makedirs(tempdir)
    except OSError:
        pass
    ## Create the output directory
    OutDir = OutDir +  "/" + ProjectName + "/"
    try:
        #os.makedirs(OutDir)
	subprocess.call("eos root://cmseos.fnal.gov mkdir -p %s" % OutDir, shell=True)
    except OSError:
        pass

    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    os.chdir("../")
    subprocess.call("tar -czf %s/menu.tgz menu/" % tempdir, shell=True)

    tupleDict=SplitNtuple()
    os.chdir(tempdir)
    subprocess.call("tar -czf %s/ntuple.tgz ntuple" % tempdir, shell=True)

    os.chdir(curdir)

    for menu in MenuFile:
        for ntuple, nopt in tupleDict.items():
            for opt in Options.items():
                lopt = list(copy.copy(opt))
                lopt[1] += nopt["opt"]
                # print Analysis, menu, ntuple, opt
                CondorSub(Analysis, menu, ntuple, tuple(lopt), nopt["cnt"])


def SplitNtuple():
    global Ntupledict

    filelistdir = tempdir + '/' + "ntuple"
    try:
        os.makedirs(filelistdir)
    except OSError:
        pass

    for n in Ntuplelist:
        Ntupledict[n]=""


    if splitline <= 0:
        for ntuple in Ntupledict.keys():
            shutil.copyfile(ntuple, "%s/%s" % (tempdir, ntuple))
        return Ntupledict

    splitedfiles = defaultdict(dict)
    for ntuple, opt in Ntupledict.items():
        f = open(ntuple, 'r')
        lines = f.readlines()
        lineperfile = splitline
        fraction = len(lines) / lineperfile
        key = os.path.splitext(ntuple)[0]

        for i in range(0, fraction):
            wlines = []
            if i == fraction - 1 :
                wlines = lines[lineperfile*i :]
            else:
                wlines = lines[lineperfile*i : lineperfile*(i+1)]
            if len(wlines) > 0:
                outf = open("%s/%s_%d.list" % (tempdir, key, i), 'w')
                outf.writelines(wlines)
            outf.close()

        splitedfiles["%s_$(Process).list" % key] = {
            "opt" : opt,
            "cnt" : fraction,
        }
    return splitedfiles

def my_CheckFile():
    global DelDir
    global DelExe
    ## Check the Delphes Dir
    DelDir = "%s/src/L1TriggerDPG/L1Menu/macros" % os.environ['CMSSW_BASE']

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
