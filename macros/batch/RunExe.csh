#!/bin/csh -v

set CMSSW = DELDIR
set EXE   = DELEXE
set OUTPUT = OUTDIR

#============================================================================#
#-----------------------------   Setup the env   ----------------------------#
#============================================================================#
echo "============  Running on" $HOST " ============"
cd ${CMSSW}/src
eval `scramv1 runtime -csh`
cd ${_CONDOR_SCRATCH_DIR}

#============================================================================#
#--------------------------   To Run the Process   --------------------------#
#============================================================================#
ls
tar -xzf ntuple.tgz
tar -xzf menu.tgz
echo ./$EXE $argv
./$EXE $argv

if ($? == 0) then
  foreach outfile (`ls results`)
    echo "Copying ${outfile} to ${OUTPUT}"
    xrdcp results/${outfile} "root://cmseos.fnal.gov/${OUTPUT}"
    if ($? == 0) then
      rm $outfile
    endif
  end
endif
