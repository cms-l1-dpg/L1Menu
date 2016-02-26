L1TMenu
=======

## Check out the code
Package to put together code and configuration file to prepare first version of the 2016 L1 menu

#### The Stage2 emulator and L1Ntuples package

<pre><code>
cmsrel CMSSW_8_0_0_pre6
cd CMSSW_8_0_0_pre6/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic cms-l1t-offline:l1t-tsg-v3
git cms-addpkg L1Trigger/L1TCommon
scram b -j 8
</code></pre>

#### L1Menu DPG package for menu-making 
<pre><code>
git clone -b 2016-Tune https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
cd L1TriggerDPG/L1Menu
scramv1 b -j 9
cd macros
make -j 9
</code></pre>


## Menu Tuning

Beside the old code inherited from 2015 L1Menu, a new executable *testMenu2016* is rewritten for easier tunning process.

<pre><code>
$ ./testMenu2016 --help
Allowed options:
  -h [ --help ]                         produce help message
  -m [ --menufile ] arg (=menu/Menu_259721_Stage2_Tune.txt)
                                        set the input menu
  -l [ --filelist ] arg (=ntuple/Run259721_stage2_Len_New.list)
                                        set the input ntuple list
  -o [ --outfilename ] arg (=Auto)      set output file name
  -d [ --outputdir ] arg (=results)     set output directory
  -t [ --writetext ] arg (=1)           write rate to output
  -c [ --writecsv ] arg (=1)            write rate to output in CSV format
  -p [ --writeplot ] arg (=1)           write plot to output
  --doPlotRate arg                      save rate plot to output
  --doPlotEff arg                       save efficiency plot to output
  --doPrintLS arg                       print out rate per LS to file
  --doPrintPU arg                       print out rate per PU to file
  -n [ --maxEvent ] arg (=-1)           run number of events; -1 for all
  -b [ --nBunches ] arg                 set number of bunches
  --SumJetET arg                        PT threshold of Jet for HT
</code></pre>

To reproduce the current menu for 1.2E34 Hz/cm^2, you can run the following command:
```
./testMenu2016 -b 4136
```
> Run259721 has 517 bunches, corresponding to 1.5E33 Hz/cm^2. To scale to
> 1.2E34Hz/cm^2, you need to scale 517 * (1.2E34/1.5E33) = 4136.
