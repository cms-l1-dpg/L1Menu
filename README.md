
L1TMenu
=======

Package to put togheter code and configuration file to prepare Winter13 version of the 2015 L1 menu

Installation instructions:

<pre><code>
export MY_CMSSW_VERSION="CMSSW_6_2_5"
cmsrel $MY_CMSSW_VERSION 
cd $MY_CMSSW_VERSION/src

cmsenv

# The PR for the RCT in 62X from Maria
git cms-merge-topic 2525

# The UCT2015 code
git clone https://github.com/uwcms/UCT2015.git L1Trigger/UCT2015
cd L1Trigger/UCT2015
git checkout 2014-Menus-V45
cd ../..

git clone https://github.com/cms-l1-dpg/L1Ntuples.git L1TriggerDPG/L1Ntuples
git clone https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu

export USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable"

scramv1 b -j 9
</code></pre>
