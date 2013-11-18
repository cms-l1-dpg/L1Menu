L1TMenu
=======

Package to put togheter code and configuration file to prepare Winter13 version of the 2015 L1 menu

Installation instructions:

<pre><code>
cmsrel CMSSW_5_3_13_patch1
cd CMSSW_5_3_13_patch1/src

cmsenv

git cms-addpkg L1Trigger/RegionalCaloTrigger       
git cms-addpkg DataFormats/L1CaloTrigger
git cms-addpkg L1TriggerConfig/L1ScalesProducers
git clone https://github.com/uwcms/UCT2015.git L1Trigger/UCT2015
git clone https://github.com/cms-l1-dpg/L1Ntuples.git L1TriggerDPG/L1Ntuples
git clone https://github.com/cms-l1-dpg/L1TMenu.git L1TriggerDPG/L1TMenu

patch -N -p0 < L1Trigger/UCT2015/eic9bit.patch

scramv1 b -j 9
</code></pre>
