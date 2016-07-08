#!/usr/bin/env python
# encoding: utf-8

# File        : PUDep.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2016 Jul 01
#
# Description :



def Menu7E33(x):
    return 12.63 + 0.54 * x + 0.110 * x*x


def Menu8p5E33(x):
    # return -5.26 + 2.36 * x + 0.028 * x*x
    return -8.20 + 2.53 * x + 0.009 * x*x

def Menu1E34(x):
    # return -5.31 + 2.20 * x + 0.023 * x*x
    return -3.85 + 1.92 * x + 0.010 * x*x

def Menu1p15E34(x):
    return -1.70 + 1.46 * x + 0.011 * x*x

def Menu1p3E34(x):
    return -0.73 + 1.25 * x + 0.009 * x*x

def Menu1p5E34(x):
    return 2.78 + 0.74 * x + 0.013 * x*x

# def Menu1E34v2(x):
    # # return -6.55 + 2.42 * x + 0.009 * x*x
    # return -2.86 + 1.91 * x + 0.023 * x*x

def Menu1p1E34(x):
    # return -1.76 + 1.66 * x + 0.016 * x*x
    return -0.35 + 1.41 * x + 0.023 * x*x

def Menu1p2E34(x):
    return 0.62 + 1.24 * x + 0.022 * x*x


# print "7E33   :  %f" %  (Menu7E33(25)   / Menu8E33(22.75)  * 2064)
# print "8.5E33 :  %f" %  (Menu8E33(30)   / Menu8E33(30.9)   * 2064)
print "1E34   :  %f" %  (Menu1E34(35) / Menu1E34(30.9) * 2064)
print "1.15E34:  %f" %  (Menu1p15E34(40) / Menu1p15E34(30.9) * 2064)
print "1.3E34 :  %f" %  (Menu1p3E34(45) / Menu1p3E34(30.9) * 2064)
print "1.5E34 :  %f" %  (Menu1p5E34(52) / Menu1p5E34(30.9) * 2064)

# print Menu1E34v2(35)
# print Menu1E34v2(30.9)
# print "1E34   :  %f" %  (Menu1E34(35) / Menu1E34(30.9) * 2064)
# print "1.15E34   :  %f" %  (Menu1E34(40) / Menu1E34(30.9) * 2064)
# print "1.3E34   :  %f" %  (Menu1E34(45) / Menu1E34(30.9) * 2064)


# print "7E33   :  %f" %  (Menu7E33(25)   / Menu8E33(22.75)  * 1)
