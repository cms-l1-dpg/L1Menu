
L1TMenu
=======

Package to put together code and configuration file to prepare Winter13 version of the 2015 L1 menu

Installation instructions:

<pre><code>
export MY_CMSSW_VERSION="CMSSW_7_2_0_pre6"
cmsrel $MY_CMSSW_VERSION 
cd $MY_CMSSW_VERSION/src

cmsenv

git clone https://github.com/cms-l1-dpg/L1Ntuples.git L1TriggerDPG/L1Ntuples
git clone https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu

export USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable"

scramv1 b -j 9
</code></pre>

In order to run on the L1Menu macro

cd L1TriggerDPG/L1Menu/macros
root
.x  ../../L1Ntuples/macros/initL1Analysis.C+
.L L1Menu2015_minbias_cross_section.C+
RunL1(true,true,4)

(4 being the scenario you want to test)
