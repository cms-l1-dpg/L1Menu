#!/bin/env python
#
# Author:   Takashi MATSUSHITA
#
# Usage:
# python getL1PrescalesRun1.py --run 262163 --user <database account>
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

  options = parser.parse_args()

  if not options.passwd:
    import getpass
    options.passwd = getpass.getpass()

  con = cx_Oracle.connect('%s/%s@%s' % (options.user, options.passwd, options.database))
  cur = con.cursor()


  menu = {}
  query = """select algo_index, alias from cms_gt.gt_run_algo_view where runnumber = %s order by algo_index""" % options.run
  cur.execute(query)
  rc = cur.fetchall()
  for x in rc:
    menu[x[0]] = x[1]

  prescales = {}
  for idx in menu.keys():
    print "%3d %40s" % (idx, menu[idx]),
    query = """select prescale_index, prescale_factor_algo_%03d from cms_gt.gt_run_presc_algo_view where runnr = %s order by prescale_index""" % (idx, options.run)
    cur.execute(query)
    rc = cur.fetchall()
    prescales = {}
    for x in rc:
      prescales[x[0]] = x[1]

    for k, v in prescales.iteritems():
      print "%6d" % v,
    print

# eof
