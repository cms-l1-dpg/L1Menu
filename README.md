
L1TMenu
=======

Package to put together code and configuration file to prepare first version of the 2016 L1 menu

Instructions for running the Stage2 emulator and create L1Ntuples:

<pre><code>
export MY_CMSSW_VERSION="CMSSW_7_6_0"
cmsrel $MY_CMSSW_VERSION 
cd $MY_CMSSW_VERSION/src
cmsenv
git cms-merge-topic cms-l1t-offline:l1t-dev-recipe-CMSSW_7_6_0
scram b -j7
cmsRun L1Trigger/L1TCommon/test/reEmul.py  max=1000
</code></pre>

##### L1Menu DPG package for menu-makeing ####################
<pre><code>
git clone -b 2016-devel https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
cd L1TriggerDPG/L1Menu
export USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable"
scramv1 b -j 9
</code></pre>

####################
In order to run the BasicRatePlots macro

<pre><code>
cd L1TriggerDPG/L1Menu/macros
root
.L BasicRatePlots.C+
goRatePlots("RUN256843_Stage2",0,20000)
</code></pre>
("RUN256843_Stage2" is the sample identifier, 0=rates (1=cross sections), 20000=number of events to process)
####################

####################
In order to run on the L1Menu macro

<pre><code>
cd L1TriggerDPG/L1Menu/macros
root
.x  ../../L1Ntuples/macros/initL1Analysis.C+
.L L1Menu2016_minbias_cross_section.C+
RunL1(true,true,4)
</code></pre>
(4 being the scenario you want to test)
####################
In order to run a simple macro that runs over the events in the ntuple and produces simple plots:

<pre><code>
cd L1TriggerDPG/L1Menu/macros
root
.x analise_L1.C
RunL1(4)
</code></pre>
(4 being the scenario you want to test)
####################
