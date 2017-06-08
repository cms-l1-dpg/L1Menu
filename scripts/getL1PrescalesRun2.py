#!/bin/env python
#
# Author:   Takashi MATSUSHITA
#
# Usage:
# python getL1PrescalesRun2.py --run 285091 --user <database account>
#

import argparse
import cx_Oracle

MAX_INDEX = 512


def mask_file(run, menu, mask):
  rows = ""
  for ii in range(MAX_INDEX):
    if ii in mask:
      rows += "        <row>%d,%d</row> <!-- %s -->\n" % (ii, mask[ii], menu[ii])
    else:
      rows += "        <row>%d,%d</row>\n" % (ii, 0)
  rows = rows.strip('\n')

  template = """<run-settings id="uGT">
  <context id="uGtProcessor">
    <param id="finorMask" type="table">
      <columns>
 	      algo,mask
      </columns>
      <types>
 	      uint,uint
      </types>
      <rows>
%s
      </rows>
    </param>
  </context>
</run-settings>
"""

  output = file('mask-run%d.xml' % run, 'w')
  output.write(template % rows)
  output.close()

  return


def prescale_file(run, menu, prescales):
  header = "      <columns>algo/prescale-index"
  types = "      <types>uint"
  for ii in range(len(prescales[prescales.keys()[0]].keys())):
    header += ", %d" % ii
    types += ", uint"
  header += "</columns>"
  types += "</types>"

  rows = ""
  for ii in prescales:
    if ii in menu:
      rows += "        <row>%d" % ii
      for key in prescales[ii]:
        ps = prescales[ii][key]
        rows += ",%d" % ps
      rows += "</row> <!-- %s -->\n" % menu[ii]
  rows = rows.strip('\n')

  prescale = {}
  prescale['header'] = header
  prescale['rows'] = rows
  prescale['types'] = types

  template = """<?xml version="1.0" encoding="UTF-8"?>
<run-settings id="uGT">
  <context id="uGtProcessor">
    <param id="index" type="uint">0</param>
    <param id="prescales" type="table">
%(header)s
%(types)s
      <rows>
%(rows)s
      </rows>
    </param>
  </context>
</run-settings>
"""

  output = file('prescale-run%d.xml' % run, 'w')
  output.write(template % prescale)
  output.close()

  return


def csv_file(run, menu, prescales):
  header = "id,name,proposed name"
  for ii in range(len(prescales[prescales.keys()[0]].keys())):
    header += ",%d" % ii
  header += ",xml,POG,PAG,comment"

  rows = ""
  for ii in prescales:
    if ii in menu:
      rows += "%d,%s," % (ii, menu[ii])
      for key in prescales[ii]:
        ps = prescales[ii][key]
        rows += ",%d" % ps
      rows += ",,,,\n"
  rows = rows.strip('\n')

  csv = {}
  csv['header'] = header
  csv['rows'] = rows

  template = """%(header)s
%(rows)s
"""

  output = file('run%d.csv' % run, 'w')
  output.write(template % csv)
  output.close()

  return


if __name__ == '__main__':
  run = None
  database = 'cms_omds_lb'
  user = None
  passwd = None

  parser = argparse.ArgumentParser()

  parser.add_argument("--run", dest="run", default=run, type=int, action="store", required=True, help="run number to process")
  parser.add_argument("--db", dest="database", default=database, type=str, action="store", help="database connection")
  parser.add_argument("--user", dest="user", default=user, type=str, action="store", required=True, help="database account user name")
  parser.add_argument("--passwd", dest="passwd", default=passwd, type=str, action="store", help="password")
  parser.add_argument("--apply-mask", dest="apply_mask", action="store_true", help="apply mask")
  parser.set_defaults(apply_mask = False)

  options = parser.parse_args()

  if not options.passwd:
    import getpass
    options.passwd = getpass.getpass()

  con = cx_Oracle.connect('%s/%s@%s' % (options.user, options.passwd, options.database))
  cur = con.cursor()


  menu = {}
  mask = {}
  query = """select algo_index, algo_name, algo_mask from cms_ugt_mon.view_ugt_run_algo_setting where run_number = %s order by algo_index""" % options.run
  cur.execute(query)
  rc = cur.fetchall()
  for x in rc:
    menu[x[0]] = x[1]
    mask[x[0]] = x[2]

  prescales = {}
  for idx in menu.keys():
    print "%3d %40s" % (idx, menu[idx]),
    query = """select prescale_index, prescale from cms_ugt_mon.view_ugt_run_prescale where algo_index = %s and run_number = %s order by prescale_index""" % (idx, options.run)
    cur.execute(query)
    rc = cur.fetchall()
    row = {}
    for x in rc:
      row[x[0]] = x[1]
      if options.apply_mask:
        row[x[0]] = x[1] * mask[idx]

    for k, v in row.iteritems():
      print "%6d" % v,
    print

    prescales[idx] = row

  mask_file(options.run, menu, mask)
  prescale_file(options.run, menu, prescales)
  csv_file(options.run, menu, prescales)

# eof
