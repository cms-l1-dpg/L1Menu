#!/bin/env python
#
# Author:   Takashi MATSUSHITA
#
# Usage:
# python getUnprescaledJson.py --json <Cert_..._JSON.txt> --user <database account> --seed <L1 seed name>
# to produce <L1 seed name>_unprescaled.json file
#

import argparse
import itertools
import json
import cx_Oracle


def toRange(i):
  for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
    b = list(b)
    yield b[0][1], b[-1][1]


def getMenu(data):
  data.id2name = {}
  data.mask = {}
  data.name2id = {}

  query = """select algo_index, algo_name, algo_mask from cms_ugt_mon.view_ugt_run_algo_setting where run_number = %s order by algo_index""" % data.run
  data.cursor.execute(query)
  rc = data.cursor.fetchall()

  for x in rc:
    name = x[1].replace('"', '')
    data.id2name[x[0]] = name
    data.mask[x[0]] = x[2]
    data.name2id[name] = x[0]


def getPrescales(data):
  data.prescales = {}

  for idx in data.id2name.keys():
    query = """select prescale_index, prescale from cms_ugt_mon.view_ugt_run_prescale where algo_index = %s and run_number = %s order by prescale_index""" % (idx, data.run)
    data.cursor.execute(query)
    rc = data.cursor.fetchall()

    data.prescales[idx] = {}
    for x in rc:
      data.prescales[idx][x[0]] = x[1] * data.mask[idx]


def getPrescaleColumnSequence(data):
  data.ps_sequence = {}

  query = """select lumi_section, prescale_index from cms_ugt_mon.view_lumi_sections where run_number = %s order by lumi_section""" % (run,)
  data.cursor.execute(query)
  rc = data.cursor.fetchall()

  end = 0
  for x in rc:
    data.ps_sequence[x[0]] = x[1]
    end = max(end, x[0])

  # fix for missing lumi-section/ps column idx
  prev = data.ps_sequence[max(data.ps_sequence.keys())]
  for ii in range(end, 0, -1):
    if ii in data.ps_sequence:
      if data.ps_sequence[ii] != None:
        prev = data.ps_sequence[ii]
      else:
        data.ps_sequence[ii] = prev # fixing None
    else:
      data.ps_sequence[ii] = prev # setting missing ls


def getUnprescaled(data, ls_ranges, seed):
  data.unprescaled[run] = None

  prescaled = []
  unprescaled = []

  prev = data.ps_sequence[min(data.ps_sequence.keys())]
  for x in ls_ranges:
    begin, end = x
    for ls in range(begin, end+1):
      if ls not in data.ps_sequence:
        if data.prescales[data.name2id[seed]][prev] != 1:
          prescaled.append(ls)
        else:
          unprescaled.append(ls)
      else:
        prev = data.ps_sequence[ls]
        if data.prescales[data.name2id[seed]][data.ps_sequence[ls]] != 1:
          prescaled.append(ls)
        else:
          unprescaled.append(ls)

  data.unprescaled[run] = list(toRange(unprescaled))

  if len(prescaled):
    print '     prescaled in', list(toRange(prescaled))


class Object:
  def __init__(self):
    pass


if __name__ == '__main__':
  database = 'cms_omds_lb'
  user = None
  passwd = None

  parser = argparse.ArgumentParser()

  parser.add_argument("--json", dest="json", default=None, type=str, action="store", required=True, help="json file")
  parser.add_argument("--seed", dest="seed", default=None, type=str, action="store", required=True, help="seed name")
  parser.add_argument("--db", dest="database", default=database, type=str, action="store", help="database connection")
  parser.add_argument("--user", dest="user", default=user, type=str, action="store", required=True, help="database account user name")
  parser.add_argument("--passwd", dest="passwd", default=passwd, type=str, action="store", help="password")

  options = parser.parse_args()

  fp = file(options.json)
  run_ls = json.load(fp)

  if not options.passwd:
    import getpass
    options.passwd = getpass.getpass()


  data = Object()
  con = cx_Oracle.connect('%s/%s@%s' % (options.user, options.passwd, options.database))
  data.cursor = con.cursor()

  data.unprescaled = {}

  hasPrescaled = False
  print 'inf> checking %s in ...' % (options.seed,)
  for run in sorted(run_ls.keys()):
    print 'inf> run %s' % (run,)
    ls_ranges = run_ls[run]
    data.run = run

    getMenu(data)
    getPrescales(data)
    getPrescaleColumnSequence(data)
    getUnprescaled(data, ls_ranges, options.seed)

  fp = open(options.seed + '_unprescaled.json', 'w') 
  json.dump(data.unprescaled, fp)

# eof
