#!/usr/bin/env python
# encoding: utf-8

# File        : Config.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2016 Feb 16
#
# Description : 


S1S2Map = {
    "L1_SingleIsoEG30er" : "L1_SingleIsoEG34er" ,
    "L1_SingleIsoEG25er" : "L1_SingleIsoEG27er" ,
    "L1_SingleIsoEG25"   : "L1_SingleIsoEG27"   ,
    "L1_SingleIsoEG22er" : "L1_SingleIsoEG24er" ,
    "L1_SingleIsoEG20er" : "L1_SingleIsoEG23er" ,
    "L1_SingleIsoEG20"   : "L1_SingleIsoEG23"   ,
    "L1_SingleIsoEG18er" : "L1_SingleIsoEG20er" ,
    "L1_SingleEG40"      : "L1_SingleEG45"      ,
    "L1_SingleEG35"      : "L1_SingleEG40"      ,
    "L1_SingleEG30"      : "L1_SingleEG34"      ,
    "L1_SingleEG25"      : "L1_SingleEG27"      ,
    "L1_SingleEG20"      : "L1_SingleEG23"      ,
    "L1_SingleEG15"      : "L1_SingleEG18"      ,
    "L1_SingleEG10"      : "L1_SingleEG17"      ,
    "L1_DoubleEG_22_10"  : "L1_DoubleEG_24_10"  ,
    "L1_DoubleEG_15_10"  : "L1_DoubleEG_18_10"  ,
    "L1_HTT100"  : "L1_HTT200"  ,
}

S2S1Map = dict((v, k) for k, v in S1S2Map.iteritems())
DualMap = S2S1Map.copy()
DualMap.update(S1S2Map)
