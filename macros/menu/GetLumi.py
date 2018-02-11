#!/usr/bin/env python
# encoding: utf-8

# File        : GetLumi.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2016 Aug 10
#
# Description :

import subprocess
import pandas as pd
import sys
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

runlist = [
    # 259721,
    # 259626,
    # 258425,
    # 258427,
    # 258428,
    # 258434,
    # 258440,
    # 258445,
    # 258448,
    # 274968,
    # 274969,
    # 274971,
    # 274998,
    # 274999,
    # 275001,
    # 275067,
    # 275066,
    # 275067,
    # 275068,
    # 275073,
    # 275074,
    # 275124,
    # 275125,
    # 276653,
    # 277069,
    # 277194,

    # 277216,
    # 277217,
    # 277218,
    # 277219,
    # 277220,
    # 277305,
    # 277420,
    # 278345,
    # 278346,
    # 278349,
    # 278366,
    # 278406,
    # 278509,
    # 278820,
    # 278821,
    # 278822,
    # 278923,
    # 278969,
    # 278975,
    # 278976,
    # 278986,

    # 277069,
    # 277070,
    # 277071,
    # 277072,
    # 277076,
    # 277087,
    # 277094,
    # 277096,
    # 277112,
    # 277180,
    # 277194,
    # 277202,
    # 278769,
    # 278770,
    # 278801,
    # 278803,
    # 278808,


   # 279862,
   # 279931,
   # 279966,
   # 279975,
   # 279993,
   # 279994,
   # 279995,
   # 280002,
   # 280006,
   # 280007,
   # 280013,
   # 280014,
   # 280015,

	306091,
	306092,
	306093
]

def Runcmd(run):
    # cmd = "brilcalc lumi --byls"
    cmd = "brilcalc lumi --byls -u '1e30/cm2s' "
    cmd += " --output-style csv -b 'STABLE BEAMS' "
    cmd += " -r %d " % run
    #testcmd = cmd + " --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json "
    testcmd = cmd + " -i /afs/cern.ch/user/d/deguio/public/Certification/Cert_13TeV_2017_HCAL_DCS_GOOD.txt "
    pipe = subprocess.Popen(testcmd, shell=True, stdout=subprocess.PIPE)
    out, err = pipe.communicate()
    if len(out) == 221: # Empty output
        pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        out, err = pipe.communicate()
    Pandas(out)

def Pandas(file):
    # df = pd.read_csv(file, header=1, skipfooter=3)
    # df = pd.read_csv(StringIO(file), header=1, skipfooter=3)
    df = pd.read_csv(StringIO(file), header=1, skipfooter=3, engine='python')
    df.loc[:,"run"] = df["#run:fill"].str.split(':').str[0]
    df.loc[:,"fill"] = df["#run:fill"].str.split(':').str[1]
    df.loc[:,"LS"] = df["ls"].str.split(':').str[0]
    df2=df[['fill','run','LS', 'avgpu']]
    df2.to_csv('lumi5.csv', mode='a', header=False, index=False)
    # df2.to_csv('run_lumi_old.csv', mode='a', header=False, index=False)
    # df.to_csv('lumi_new.csv', mode='a', header=False, index=False)

if __name__ == "__main__":
    for run in runlist:
        Runcmd(run)

