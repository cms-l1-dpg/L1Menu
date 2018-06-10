L1TMenu
=======

### HowToL1TriggerMenu twiki:
### https://twiki.cern.ch/twiki/bin/view/CMS/HowToL1TriggerMenu
You can find useful information, like how to customise seeds, how to make menulib and locations of ntuples and so on

### Checkout the code
<pre><code>
cmsrel CMSSW_9_2_15
cd CMSSW_9_2_15/src
cmsenv
git clone --branch hui https://github.com/cms-l1-dpg/L1Menu.git  L1TriggerDPG/L1Menu
make -j 8
</code></pre>

### Pick your preferred ntuple (ntuple location is in the twiki) 
nanoDST ntuple is recommended for rate study because it has much more statistics  
Need to use ZeroBias ntuple if you want to run emulation 
<pre><code>
 cd ntuple
./makeFileList.py your_ntuple_location
</code></pre>

### Make Lumi Section (LS) information table with BrilCalc
This is important, because 2018 ntuples are produced without json file. So we need to pick our preferred json file to avoid bad LS  
Also, LS vs PU information is stored in this table, which will be used later, in the rate vs PU plots  
Go to the menu folder, edit GetLumi.py with the run number of your ntuple 
<pre><code>
source GetLumi_setup.sh
./GetLumi.py
</code></pre>

### Make prescale (PS) table
2018 official PS table is in https://github.com/cms-l1-dpg/L1Menu2018/tree/master/official/PrescaleTables  
The PS table is in .xlsx format, for better presentation  
Convert the PS table to either a tab separated .txt file, or a .CSV file  
menu/Prescale_2018_v1_0_0_Col_2.0.txt is an example of the latest PS table, column 2.0e34  
Add your customised seeds into the PS table, with your designed prescale value

### Check the arguments of the code
 ./testMenu2016 --help for all arguments  
some useful argument:  
-u the LS information table you just made  
-m the PS table you just customised  
-l the ntuple list you just made  
-o name of output files  
-b number of bunches. Usually is 2544 for 2018 data  
--UseUnpackTree to use UnpackTree, the default is EmuTree  
--SelectRun to select the run number if your ntuple list has multiple runs, the default is the whole ntuple list  
--SelectLS to select the LS, the defualt is the whole LS  
You can loop up the LS information table for help. For nanoDST ntuple, the LS corresponding to PU 55-57 (for col 2.0e34) sould be enough. If your ntuple does not have that high PU, you can do a linear exreapolation

### Run the code with arguments
<pre><code>
./testMenu2016 -u menu/run_lumi.csv -m menu/Prescale_2018_v1_0_0_Col_2.0.txt -l ntuple/your_ntuple.list -o name_of_output_files -b 2544 --doPlotRate --UseUnpackTree --SelectRun your_run_number --SelectLS '[start_LS,end_LS]'
</code></pre>

### Make rate vs PU plots
add argument --doPrintPU  
remove argument --SelectRun and --SelectLS, because we need to cover the full PU  
You can use batch/SubmitLPC.py or SubmitLSF.py to split the job  
edit plot/CompPUDep.py and run it!
