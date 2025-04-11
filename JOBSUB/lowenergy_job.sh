#!/bin/bash
export UPS_OVERRIDE="-H Linux64bit+3.10-2.17"
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw v10_04_05d00 -q e26:prof
setup gallery v1_22_06 -q e26:prof
#setup ifdhc

timenow=`date +%b_%d_%Y_%H_%M`
LOGFILE=LOG_${timenow}_${PROCESS}.log
#some basic info
echo Start  `date` >> ${LOGFILE} 2>&1
echo "the worker node is " `hostname` "OS: "  `uname -a` >> ${LOGFILE} 2>&1
echo "the user id is " `whoami` >> ${LOGFILE} 2>&1
echo "the output of id is " `id` >> ${LOGFILE} 2>&1
USER_SCRATCH_DIR=/pnfs/dune/scratch/users/sfogarty/JOBSUB/
OUTDIR=${USER_SCRATCH_DIR}/CLMatching/${RUN}/
ifdh cp ${USER_SCRATCH_DIR}/CL_lowenergy_analysis.C $_CONDOR_SCRATCH_DIR
cd $_CONDOR_SCRATCH_DIR

RUN=$1
SINGLEHIT_FILELIST=$2
RECO_FILELIST=$3
OUTDIR=$4
CHANNELMAP=$5
OPTAG=$6
OPTWINDOW=$7
OPZWINDOW=$8
OPYWINDOW=$9
MATCHTWINDOWLOW=${10}
MATCHTWINDOWHIGH=${11}
QZWINDOW=${12}
QYWINDOW=${13}
PROCESS_OFFSET=${14}
PROCESS=$(( ${PROCESS} + ${PROCESS_OFFSET} ))

ifdh cp ${USER_SCRATCH_DIR}/${CHANNELMAP} $_CONDOR_SCRATCH_DIR
ifdh cp ${USER_SCRATCH_DIR}/${SINGLEHIT_FILELIST} $_CONDOR_SCRATCH_DIR
ifdh cp ${USER_SCRATCH_DIR}/${RECO_FILELIST} $_CONDOR_SCRATCH_DIR
SINGLEHIT_FILE=$(sed -n $((PROCESS + 1))p ${SINGLEHIT_FILELIST})
ifdh cp ${SINGLEHIT_FILE} $_CONDOR_SCRATCH_DIR

#mkdir -p ${USER_SCRATCH_DIR}/CLMatching/
#mkdir -p ${USER_SCRATCH_DIR}/CLMatching/${RUN}
OUTFILE=PDHD_CLMATCHING_RESULTS_${RUN}_jobnum${PROCESS}.root

echo "Run = ${RUN}" >> ${LOGFILE} 2>&1
echo "SINGLEHIT_RUNLIST = ${SINGLEHIT_FILELIST}" >> ${LOGFILE} 2>&1
echo "RECO_FILELIST = ${RECO_FILELIST}" >> ${LOGFILE} 2>&1
echo "OUTDIR = ${OUTDIR}" >> ${LOGFILE} 2>&1
echo "OUTFILE = ${OUTFILE}" >> ${LOGFILE} 2>&1
echo "CHANNELMAP = ${CHANNELMAP}" >> ${LOGFILE} 2>&1
echo "OPTAG = ${OPTAG}" >> ${LOGFILE} 2>&1
echo "OPTWINDOW = ${OPTWINDOW}" >> ${LOGFILE} 2>&1
echo "OPZWINDOW = ${OPZWINDOW}" >> ${LOGFILE} 2>&1
echo "OPYWINDOW = ${OPYWINDOW}" >> ${LOGFILE} 2>&1
echo "MATCHTWINDOWLOW = ${MATCHTWINDOWLOW}" >> ${LOGFILE} 2>&1
echo "MATCHTWINDOWHIGH = ${MATCHTWINDOWHIGH}" >> ${LOGFILE} 2>&1
echo "QZWINDOW = ${QZWINDOW}" >> ${LOGFILE} 2>&1
echo "QYWINDOW = ${QYWINDOW}" >> ${LOGFILE} 2>&1
echo "Process offset = ${PROCESS_OFFSET}" >> ${LOGFILE} 2>&1

echo "XROOTD_LIB: ${XROOTD_LIB}" >> ${LOGFILE} 2>&1
## run main script
LD_PRELOAD=$XROOTD_LIB/libXrdPosixPreload.so root -l -q -b "CL_lowenergy_analysis.C(${RUN}, ${PROCESS}, \"${SINGLEHIT_FILELIST}\", \"${RECO_FILELIST}\", \"${OUTFILE}\", \"${CHANNELMAP}\", \"${OPTAG}\", \"${OPTWINDOW}\", \"${OPZWINDOW}\", \"${OPYWINDOW}\", \"${MATCHTWINDOWLOW}\", \"${MATCHTWINDOWHIGH}\", \"${QZWINDOW}\", \"${QYWINDOW}\")"
##

ifdh cp $OUTFILE ${OUTDIR}

echo "End `date`" >> ${LOGFILE} 2>&1

ifdh cp ${LOGFILE} ${OUTDIR}
exit 0
