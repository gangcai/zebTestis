#!/bin/bash
#-V:pass all environment variables to the job
#-N: job name
#-l: specify the amount of maximum memory required
#-d: working directory
work_dir=$PWD #default current working directory
echo $work_dir
#qsub -V -N qTest -l h_vmem=5G -d $work_dir test.sh
#qsub -V -N mt -l mem=10000MB -v "arg1=1000,arg2=444" -d $work_dir run1.sh
#Rscript run_base.R --max.feature ${maxfeature} --pc.num ${pcnum} --min.feature ${minfeature} --variable.features.n ${vn} --regress.type ${rtype}
pc=15
minfeature=200
maxfeature=5000
vn=1000
rtype=1
mincell=5
mtper=5
#qsub -V -N "pc${pc}_maxfeature${maxfeature}_minfeature${minfeature}_vn${vn}_rtype${rtype}" -v "vn=${vn},rtype=${rtype},pcnum=${pc},minfeature=${minfeature},maxfeature=${maxfeature}" -l mem=10000MB -d $work_dir run.sh
#--max.feature ${maxfeature} --pc.num ${pcnum} --min.feature ${minfeature} --variable.features.n ${vn} --regress.type ${rtype} --mt.per ${mtper} --min.cell ${mincell}
qsub -V -N "qvv2"  -v "vn=${vn},rtype=${rtype},pcnum=${pc},minfeature=${minfeature},maxfeature=${maxfeature},mtper=${mtper},mincell=${mincell}" -l mem=10000MB -d $work_dir run.sh

