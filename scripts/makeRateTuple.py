#!/bin/env python
#
# Author:   Takashi MATSUSHITA
#
# Usage:
# python makeRateTuple.py --runs 278986 278976 --user <database account>
# for x in convert*.cc; do root -l -b -q $x; done   # in bash
#

import argparse
import sys

import cx_Oracle


def getLumiSections(options, run):
  query = """select lhcfill, lumisection, instlumi, pileup, physics_flag from cms_runtime_logger.lumi_sections where runnumber = {0} order by lumisection""".format(run)
  
  options.cursor.execute(query)
  rc = options.cursor.fetchall()
  lumi_section = {}
  lhc_fill = 99999999
  for x in rc:
    fill, ls, lumi, pu, flag = x
    data = {'lumi': lumi, 'pu': pu, 'flag': flag}
    lumi_section.update({ls : data})
    lhc_fill = min(lhc_fill, fill)

  return lhc_fill, lumi_section



def getFillInfo(options, fill):
  fill_summary = {}
  query = """select nbunchesbeam1, nbunchesbeam2, ncollidingbunches, ntargetbunches from cms_runtime_logger.runtime_summary where lhcfill = {0}""".format(fill)
  options.cursor.execute(query)
  rc = options.cursor.fetchall()
  for x in rc:
    nb1, nb2, ncolliding, ntarget = x
    fill_summary.update({'fill': fill, 'nb1': nb1, 'nb2': nb2, 'ncolliding': ncolliding, 'ntarget': ntarget})

  return fill_summary



def setTcdsCounts(options, run, lumi_section):
  query = """select section_number, trg_cnt_total, trg_cnt_tt1, sup_trg_cnt_total, sup_trg_cnt_tt1 from cms_tcds_monitoring.tcds_cpm_counts_v where run_number = {0} order by section_number""".format(run)
  
  options.cursor.execute(query)
  rc = options.cursor.fetchall()
  for x in rc:
    try:
      ls, total, ugt, sup_total, sup_ugt = x
      lumi_section[ls]['total'] = total
      lumi_section[ls]['ugt'] = ugt
      lumi_section[ls]['sup_total'] = sup_total
      lumi_section[ls]['sup_ugt'] = sup_ugt
    except:
      e = sys.exc_info()[0]
      v = sys.exc_info()[1]
      print "warning> {0}: {1} at LS={2} of Run={3}".format(e.__name__, v, ls, run)
 


def setPrescaleIndex(options, run, lumi_section):
  query = """select lumi_section, prescale_index from cms_ugt_mon.view_lumi_sections where run_number = {0} order by lumi_section""".format(run)
  options.cursor.execute(query)
  rc = options.cursor.fetchall()
  for x in rc:
    ls, ps = x
    lumi_section[ls]['ps'] = ps



def dumpText(run, fill_info, lumi_section, algo_counts):
  fp = open('{0}.txt'.format(run), 'wb')
  for ls, data in  lumi_section.iteritems():
    try:
      text = "{fill} {nb1} {nb2} {ncolliding} {ntarget} ".format(**fill_info)
      text += "{0} {1} ".format(run, ls)
      text += "{flag} {ps} {lumi} {pu} {total} {ugt} {sup_total} {sup_ugt}".format(**data)
      counts = algo_counts[ls]
      for name in sorted(counts.keys()):
        text += " {0}".format(counts[name])
      text += "\n";
      fp.write(text)
    except:
      e = sys.exc_info()[0]
      v = sys.exc_info()[1]
      print "warning> {0}: {1} at LS={2} of Run={3}".format(e.__name__, v, ls, run)
  fp.close()



def getMenu(options, run):
  query = """select algo_index, algo_name from cms_ugt_mon.view_ugt_run_algo_setting where run_number = {0} order by algo_index""".format(run)
  options.cursor.execute(query)
  rc = options.cursor.fetchall()
  l1menu = {}
  for x in rc:
    index, name = x
    l1menu[index] = name.strip('"')
  return l1menu



def getAlgoRates(options, run, l1menu):
  algo_counts = {}
  for index in l1menu.keys():
    print "      fetching algo rate for %s..." % l1menu[index]
    query = """select lumi_sections_id, algo_count from cms_ugt_mon.view_algo_scalers where lumi_sections_id like '%07d_' || '%%' and algo_index = %s and scaler_type = 1 order by lumi_sections_id""" % (run, index)
    options.cursor.execute(query)
    rc = options.cursor.fetchall()
    for x in rc:
      lumi_sections_id, algo_count = x
      n, ls = map(int, lumi_sections_id.split('_'))
      if ls not in algo_counts.keys():
        algo_counts.update({ls: {}})
      algo_counts[ls].update({l1menu[index]: algo_count})

  return algo_counts



def makeConverter(run, algo_counts):
  fp = open('convert%s.cc' % run, 'wb')
  line = """
  void convert%s()
  {
    TFile *tfile = new TFile("%s.root", "RECREATE");
    ttree = new TTree("rate", "");
  
    int fill;
    int nb1;
    int nb2;
    int ncolliding;
    int ntarget;
    int run;
    int ls;
    int flag;
    int ps;
    float lumi;
    float pu;
    int total;
    int ugt;
    int sup_total;
    int sup_ugt;
  
    ttree->Branch("lhcFill", &fill);
    ttree->Branch("nBunchesBeam1", &nb1);
    ttree->Branch("nBunchesBeam2", &nb2);
    ttree->Branch("nCollidingBunches", &ncolliding);
    ttree->Branch("nTargetBunches", &ntarget);
    ttree->Branch("runNumber", &run);
    ttree->Branch("lumiSection", &ls);
    ttree->Branch("physicsFlag", &flag);
    ttree->Branch("prescaleColumnIndex", &ps);
    ttree->Branch("luminosity", &lumi);
    ttree->Branch("pileUp", &pu);
    ttree->Branch("totalL1a", &total);
    ttree->Branch("l1aByUgt", &ugt);
    ttree->Branch("suppressedTotal", &sup_total);
    ttree->Branch("suppressedUgt", &sup_ugt);
  """
  fp.write(line % (run, run))

  declare = "";
  for name in sorted(algo_counts[1].keys()):
    declare += 'int {0}; ttree->Branch("{0}", &{0});\n'.format(name)
  fp.write(declare);
  
  read = "";
  for name in sorted(algo_counts[1].keys()):
    read += ' >> {0}'.format(name)

  line = """
    std::ifstream input("%s.txt");
    std::string line;
    while (std::getline(input, line))
    {
      std::istringstream iss(line);
      iss >> fill >> nb1 >> nb2 >> ncolliding >> ntarget >> run >> ls >> flag >> ps >> lumi >> pu >> total >> ugt >> sup_total >> sup_ugt
          %s;
      ttree->Fill();
    }
  
    tfile->Write();
    tfile->Close();
  }
  """ % (run, read)
  fp.write(line)
  fp.close()



if __name__ == '__main__':
  runs = None
  database = 'cms_omds_lb'
  user = None
  passwd = None

  parser = argparse.ArgumentParser()

  parser.add_argument("--runs", dest="runs", default=runs, nargs="+", type=int, action="store", required=True, help="run numbers to process [separated with space]")
  parser.add_argument("--db", dest="database", default=database, type=str, action="store", help="database connection")
  parser.add_argument("--user", dest="user", default=user, type=str, action="store", required=True, help="database account user name")
  parser.add_argument("--passwd", dest="passwd", default=passwd, type=str, action="store", help="password")

  options = parser.parse_args()

  if not options.passwd:
    import getpass
    options.passwd = getpass.getpass()

  options.connection = cx_Oracle.connect('%s/%s@%s' % (options.user, options.passwd, options.database))
  options.cursor = options.connection.cursor()

  for run in options.runs:
    print "inf> processing run %s" % run
    fill, lumi_sections = getLumiSections(options, run)
    fill_info = getFillInfo(options, fill)
    setTcdsCounts(options, run, lumi_sections)
    setPrescaleIndex(options, run, lumi_sections)
    l1menu = getMenu(options, run)
    algo_counts = getAlgoRates(options, run, l1menu)
    dumpText(run, fill_info, lumi_sections, algo_counts)
    makeConverter(run, algo_counts)

# eof
