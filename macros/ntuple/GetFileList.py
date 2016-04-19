#!/usr/bin/env python
# encoding: utf-8

# File        : GetFileList.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2016 Apr 02
#
# Description :

import subprocess
import re
from makeFileList import EOS
from collections import defaultdict

ntupleMap = {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSG-v3 ~~~~~
    # "r*_tsgv3"          : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v3/ZB/crab_l1-tsg-v3__*_ZB/",
    # "SingleMuZmu_tsgv3" : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v3/SingleMuon/crab_l1-tsg-v3__SingleMuon_ZMu/160225_175820/0000/",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSG-v4 ~~~~~
    # "r*_tsgv4"          : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/ZB/crab_l1-tsg-v4__*_ZB/",
    # "SingleMuZmu_tsgv4" : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/SingleMuon/crab_l1-tsg-v4__SingleMuon_ZMu/160314_132813/0000/",
    # "MC10PU_tsgv4"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/SingleNeutrino/crab_l1-tsg-v4__SingleNeutrino_25nsPU10/160402_155410/0000",
    # "MC20PU_tsgv4"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/SingleNeutrino/crab_SingleNeutrino_PU20_FixedRecipe/160315_112126/0000",
    # "MC30PU_tsgv4"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/SingleNeutrino/crab_SingleNeutrino_PU30_FixedRecipe/160315_112559/0000/",
    # "MC40PU_tsgv4"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/SingleNeutrino/crab_l1-tsg-v4__SingleNeutrino_25nsPU40/160401_213652/0000",
    # "MC50PU_tsgv4"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/SingleNeutrino/crab_l1-tsg-v4__SingleNeutrino_25nsPU50/160315_113429/0000",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TSG-v4p1 ~~~~~
    # "r*_tsgv4p1"          : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4p1/ZB/crab_l1-tsg-v4p1__*_ZB/",
    # "SingleMuZmu_tsgv4p1" : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4p1/SingleMuon/crab_l1-tsg-v4p1__SingleMuon_ZMu/160406_202706/0000/",
    # "MC20PU_tsgv4p1"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4p1/SingleNeutrino/crab_l1-tsg-v4p1__SingleNeutrino_25nsPU20/160406_202848/0000/",
    # "MC30PU_tsgv4p1"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4p1/SingleNeutrino/crab_l1-tsg-v4p1__SingleNeutrino_25nsPU30/160406_202819/0000/",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Itg-v19 ~~~~~
    # "r*_itgv19"          : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v19pt0/ZB/crab_l1t-integration-v19pt0__*_ZB/",
    # "r*_itgv19_Reco"     : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v19pt0/ZB/crab_l1t-integration-v19pt0__*_ZB_wRECO/",
    # "SingleMuZmu_itgv19" : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v19pt0/SingleMuon/crab_l1t-integration-v19pt0__SingleMuon_ZMu/160330_205037/0000/",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Itg-v19-layer1 ~~~~~
    # "r*_itgv19Layer1"          : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v19pt0-layer1/ZB/crab_l1t-integration-v19pt0-layer1__*_ZB/",
    # "r*_itgv19Layer1_Reco"     : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v19pt0-layer1/ZB/crab_l1t-integration-v19pt0-layer1__*_ZB_wRECO/",
    # "SingleMuZmu_itgv19Layer1" : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v19pt0-layer1/SingleMuon/crab_l1t-integration-v19pt0-layer1__SingleMuon_ZMu/160330_214907/0000/",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Itg-v27 ~~~~~
    # "r*_itgv27"          : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v27p0/ZB/crab_l1t-integration-v27p0__*_ZB/",
    # "SingleMuZmu_itgv27" : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v27p0/SingleMuon/crab_l1t-integration-v27p0__SingleMuon_ZMu/160408_194131/0000/",
    # "MC20PU_itgv27"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v27p0/SingleNeutrino/crab_l1t-integration-v27p0__SingleNeutrino_25nsPU20/160408_195438/0000/",
    # "MC30PU_itgv27"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v27p0/SingleNeutrino/crab_l1t-integration-v27p0__SingleNeutrino_25nsPU30/160408_195203/0000/",
    # "MC40PU_itgv27"      : "/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v27p0/SingleNeutrino/crab_l1t-integration-v27p0__SingleNeutrino_25nsPU40/160410_164601/0000/",

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Itg-v37 ~~~~~
    # "r*_itgv37p2":"/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v37p2/ZB/crab_l1t-integration-v37p2__*_ZB/",
    "r259721_itgv37p2_Opt5":"/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v37p2_Opt5Tau/ZeroBias1/crab_l1t-integration-v37p2_Opt5Tau__259721_ZeroBias1/160418_132646/0000/",
    # "r*_itgv36":"/afs/cern.ch/user/b/benwu/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v36p0_Opt5Tau/ZB/crab_l1t-integration-v36p0_Opt5Tau__*_ZB/",
}

def GetEOSdir(v):
    redir =""
    if v.find("/eos") != -1:
        findstore = v.find("/store", v.find("/eos"))
        redir = v[findstore:]
    else:
        redir = v
    return redir

def GetList(k, v, opt='w+'):
    f = open(k+".list", opt)
    eosdir = GetEOSdir(v)
    print eosdir
    p = subprocess.Popen('./makeFileList.py %s' % eosdir, shell=True, stdout=f, close_fds=True)
    rec_code = p.wait()
    f.flush()
    return rec_code

def GetNtuleList():
    for k,v  in ntupleMap.items():
        if '*' in k:
            GetZBLists(k , v)
        else:
            GetList(k, v)

def EOSls(eosdir):
    p1 = subprocess.Popen([EOS, "ls", eosdir], shell=False, stdout=subprocess.PIPE)
    (stdout, stderr)=p1.communicate()
    if stderr is not None:
        print "Trouble executing the srmls command"
        sys.exit(1)
    return stdout.split()

def GetZBlist(k, v):
    tempvlist = []
    for tempv in v:
        jdate = EOSls(tempv)
        if len(jdate) ==1:
            tempv +="/"+jdate[0]
        else:
            print "Multiple jobs: ",
            for i in range(len(jdate)):
                print "(%d) %s" % (i, jdate[i]),
            print ""
            job = input("Which one?")
            tempv += "/"+jdate[job]
        jbatch = EOSls(tempv)
        for j in jbatch:
            tempvlist.append(tempv+"/"+j)
    for i in range(len(tempvlist)):
        if i==0:
            GetList(k, tempvlist[i])
        else:
            GetList(k, tempvlist[i], 'a')

def GetZBLists(k, v):
    zbhomedir = GetEOSdir(v[:v.find("/ZB/")])
    lsdir = EOSls(zbhomedir)
    runmap = defaultdict(list)
    for zbdir in lsdir:
        if "ZeroBias" not in zbdir:
            continue
        dirpattern = GetEOSdir(v).replace("ZB", zbdir).replace("*", "(\d+)")
        if dirpattern[-1] == "/":
            dirpattern = dirpattern[:-1] + "/*"
        else:
            dirpattern = dirpattern + "/*"
        pat = re.compile(dirpattern)

        rundirs = EOSls(zbhomedir+"/"+zbdir)
        for rundir in rundirs:
            thisdir ="/".join([zbhomedir, zbdir, rundir])
            mat = pat.match(thisdir)
            if mat is not None:
                run = mat.group(1)
                runmap[k.replace("*", run)].append(thisdir)
    for k, v in runmap.items():
        GetZBlist(k, v)
        # if len(v) == 1:
            # GetList(k, v[0])
            # continue





if __name__ == "__main__":
    GetNtuleList()
