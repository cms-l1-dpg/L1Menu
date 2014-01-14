L1TMenu
=======

Package to put togheter code and configuration file to prepare Winter13 version of the 2015 L1 menu

Installation instructions:

<pre><code>
cmsrel CMSSW_5_3_14_patch2
cd CMSSW_5_3_14_patch2/src

cmsenv

git cms-addpkg L1Trigger/RegionalCaloTrigger       
git cms-addpkg DataFormats/L1CaloTrigger
git cms-addpkg L1TriggerConfig/L1ScalesProducers
git clone https://github.com/uwcms/UCT2015.git L1Trigger/UCT2015
cd L1Trigger/UCT2015
git checkout 2014-Menus-V0
cd ../..
git clone https://github.com/cms-l1-dpg/L1Ntuples.git L1TriggerDPG/L1Ntuples
git clone https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu

patch -N -p0 < L1Trigger/UCT2015/eic9bit.patch

scramv1 b -j 9
</code></pre>
