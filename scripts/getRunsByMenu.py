#!/bin/env python
#
# Author:   Takashi MATSUSHITA
#
# Usage:
# python getRunsByMenu.py --menu menu
#

import argparse
import cx_Oracle


if __name__ == '__main__':
  database = 'cms_omds_lb'
  user = None
  passwd = None
  menu = None
  
  parser = argparse.ArgumentParser()
  parser.add_argument("--menu", dest="menu", default=menu, type=str, action="store", required=True, help="menu name")
  parser.add_argument("--db", dest="database", default=database, type=str, action="store", help="database connection")
  parser.add_argument("--user", dest="user", default=user, type=str, action="store", required=True, help="database account user name")
  parser.add_argument("--passwd", dest="passwd", default=passwd, type=str, action="store", help="passwd")
  options = parser.parse_args()
  
  if not options.passwd:
    import getpass
    options.passwd = getpass.getpass()
  
  con = cx_Oracle.connect('%s/%s@%s' % (options.user, options.passwd, options.database))
  cur = con.cursor()
  
  query = """select wbm.lhcfill, wbm.runnumber, wbm.gtkey, conf.id, menu.l1_menu
             from cms_wbm.runsummary wbm
             join cms_trg_l1_conf.l1_trg_conf_keys conf on wbm.tsckey in conf.id
             join cms_trg_l1_conf.ugt_keys menu on (conf.ugt_key in menu.id)
             where conf.id like 'collisions' || '%%' and menu.l1_menu like '%s' || '%%' order by wbm.runnumber""" % options.menu
  cur.execute(query)
  rc = cur.fetchall()
  print 'Run # for %s' % options.menu
  for x in rc:
    if x[0]:
      print '  ', x[1], 'with menu:', x[4], " (", x[2], ",", x[3], ")"
  print

# eof
