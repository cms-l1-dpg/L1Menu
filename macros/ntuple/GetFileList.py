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
    # "r271336_itgv42p1" : "~/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1t-integration-v42p1/ZeroBias1/crab_l1t-integration-v42p1__271336_ZeroBias1/160426_220107/0000/",
    "r*_v78Calo" : "/eos/uscms/store/group/lpctrig/apana/L1Menu_2016/Stage2/Collision2016-noRECO-l1t-integration-v78p0_CaloStage2Params_v2_2/ZeroBias/crab_Collision2016-noRECO-l1t-integration-v78p0_CaloStage2Params_v2_2__*_ZeroBias/",
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
    cmd = " ".join( [EOS, "ls", eosdir] )
    p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    stdout, stderr= p1.communicate()
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
            # job = input("Which one?")
            # tempv += "/"+jdate[job]
            print max(jdate)
            tempv += "/"+max(jdate)
        jbatch = EOSls(tempv)
        for j in jbatch:
            tempvlist.append(tempv+"/"+j)
    for i in range(len(tempvlist)):
        if i==0:
            GetList(k, tempvlist[i])
        else:
            GetList(k, tempvlist[i], 'a')

def GetZBLists(k, v):
    import re
    rsearch = re.search('ZeroBias|ParkingZeroBias', v)
    if rsearch is None:
        return None
    zbhomedir = GetEOSdir(v[:rsearch.start()])
    lsdir = EOSls(zbhomedir)
    runmap = defaultdict(list)
    for zbdir in lsdir:
        if "ZeroBias" not in zbdir:
            continue
        dirpattern = GetEOSdir(v).replace("ZeroBias", zbdir).replace("*", "(\d+)")
        if dirpattern[-1] == "/":
            dirpattern = dirpattern[:-1] + "/*"
        else:
            dirpattern = dirpattern + "/*"
        pat = re.compile(dirpattern)

        rundirs = EOSls(zbhomedir+"/"+zbdir)
        for rundir in rundirs:
            thisdir ="/".join([zbhomedir, zbdir, rundir])
            thisdir = thisdir.replace("//", "/")
            mat = pat.match(thisdir)
            if mat is not None:
                run = mat.group(1)
                runmap[k.replace("*", run)].append(thisdir)
    for k, v in runmap.items():
        GetZBlist(k, v)

if __name__ == "__main__":
    GetNtuleList()
