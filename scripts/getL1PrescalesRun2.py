#!/bin/env python
#
# Author:   Takashi MATSUSHITA
#
# Usage:
# python getL1PrescalesRun2.py --run 285091 --user <database account>
#

import argparse
import cx_Oracle


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
    prescales = {}
    for x in rc:
      prescales[x[0]] = x[1]
      if options.apply_mask:
        prescales[x[0]] = x[1] * mask[idx]

    for k, v in prescales.iteritems():
      print "%6d" % v,
    print

# eof
