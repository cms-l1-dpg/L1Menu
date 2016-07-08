L1TMenu
=======

## Check out the code
Package to put together code and configuration file to prepare first version of the 2016 L1 menu

#### The Stage2 emulator and L1Ntuples package

<pre><code>
cmsrel CMSSW_8_0_11
cd CMSSW_8_0_11/src
cmsenv
git cms-init
</code></pre>

#### L1Menu DPG package for menu-making 
<pre><code>
git clone -b 2016-Tune git@github.com:cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
cd L1TriggerDPG/L1Menu/macros
make -j 9
</code></pre>


## Menu Tuning

Beside the old code inherited from 2015 L1Menu, a new executable *testMenu2016* is rewritten for easier tunning process.

<pre><code>
$ ./testMenu2016 --help
Allowed options:
  -h [ --help ]                         produce help message
  -m [ --menufile ] arg (=menu/Menu_259721_TSGv3_FixPre_Mu18.txt)
                                        set the input menu
  -l [ --filelist ] arg (=ntuple/r259721_tsgv3.list)
                                        set the input ntuple list
  -o [ --outfilename ] arg (=Auto)      set output file name
  -d [ --outputdir ] arg (=results)     set output directory
  -t [ --writetext ] arg (=1)           write rate to output
  -c [ --writecsv ] arg (=1)            write rate to output in CSV format
  -p [ --writeplot ] arg (=1)           write plot to output
  --doPlotRate                          save rate plot to output
  --doPlotEff                           save efficiency plot to output
  --doPlotTest                          save testing plot to output
  --doPlotuGt                           save uGT histograms to output
  --doPlotLS                            save count per LS plot to output
  --doTnPMuon                           use tag & probe for muon efficiency
  --doPrintPU                           print out rate per PU to file
  --doCompuGT                           Compare emulation with uGT tree
  -n [ --maxEvent ] arg (=-1)           run number of events; -1 for all
  -b [ --nBunches ] arg                 set number of bunches
  --SumJetET arg                        PT threshold of Jet for HT
  --SetMuonER arg                       Set the ER in eta for Muon
  --SetNoPrescale                       Set all prescales to 1
  --UseUpgradeLyr1                      Use Upgrade Layer1 Tree
  --UseL1CaloTower                      Use Layer1 CaloTower Tree
  --UsePFMETNoMuon                      Use PFMET no Muon in SingleMu sample
  --UseuGTDecision                      Trigger seeds fired by uGT
  --UseUnpackTree                       Use unpacked tree in Ntuple
  --SelectRun arg (=-1)                 Select specific run
  --SelectEvent arg (=-1)               Select specific event
  --SelectLS arg                        Select specific LS ranges
</code></pre>

#### Tips
* `--nBunches float`: Positive number will be treated as number of colliding bunches for pre-deadtime rate estimation;
                   negative number will be treated as L1 and HLT ZB prescales for post-deadtime rate estimation.
* `--SelectLS string`: SelectLS allows JSON-like format, like "[1, 30], [34, 40]"



### Tuning for Post-ICHEP
Use plot/PUDep.py for PU scaling 

```
./testMenu2016 -d NewLumi --UseUnpackTree --nBunches 4143.995262 --SelectLS '[1, 96]' -m menu/Lumi1p5E34.txt  -l ntuple/r275066_parkV62.list -o r275066_1p5E34_0 >&! r275066_1p5E34_0.log &
```
