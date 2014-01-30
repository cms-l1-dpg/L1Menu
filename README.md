
L1TMenu
=======

Package to put togheter code and configuration file to prepare Winter13 version of the 2015 L1 menu

Installation instructions:

<pre><code>
export MY_CMSSW_VERSION="CMSSW_5_3_14_patch2"
cmsrel $MY_CMSSW_VERSION 
cd $MY_CMSSW_VERSION/src

cmsenv

git cms-cvs-history import  UCT2015_v4 L1Trigger/RegionalCaloTrigger
git clone https://github.com/uwcms/UCT2015.git L1Trigger/UCT2015
cd L1Trigger/UCT2015
git checkout 2014-Menus-V1
cd ../..

git clone https://github.com/cms-l1-dpg/L1Ntuples.git L1TriggerDPG/L1Ntuples
git clone https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu

export USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable"

scramv1 b -j 9
</code></pre>
