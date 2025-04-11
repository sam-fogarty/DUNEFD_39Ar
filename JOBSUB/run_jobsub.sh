USERDIR=/exp/dune/data/users/sfogarty/DUNEFD_39Ar/JOBSUB/
USER_SCRATCH_DIR=/pnfs/dune/scratch/users/sfogarty/JOBSUB/
RUN=29218
SINGLEHIT_FILELIST=${USERDIR}/singlehit_lists/singlehit_${RUN}_filepaths.txt
RECO_FILELIST=${USERDIR}/runlists/${RUN}_filepaths.txt
#OUTDIR=${USERDIR}/CLMatching/${RUN}
folder_name=${RUN}_test2
mkdir -p ${USER_SCRATCH_DIR}/CLMatching/${folder_name}
OUTDIR=${USER_SCRATCH_DIR}/CLMatching/${folder_name}

cp -f ${SINGLEHIT_FILELIST} ${USER_SCRATCH_DIR}
cp -f ${RECO_FILELIST} ${USER_SCRATCH_DIR}
cp -f ${USERDIR}/CL_lowenergy_analysis.C ${USER_SCRATCH_DIR}

SINGLEHIT_FILELIST=singlehit_${RUN}_filepaths.txt
RECO_FILELIST=${RUN}_filepaths.txt
#mkdir -p $OUTDIR
#mkdir -p ${USER_SCRATCH_DIR}/CLMatching/

CHANNELMAP=PDHD_PDS_ChannelMap.csv
OPTAG=pdhddaphne:daq
OPTWINDOW=0.1 #us
OPZWINDOW=600 #cm
OPYWINDOW=600 #cm
MATCHTWINDOWLOW=50 #us
MATCHTWINDOWHIGH=2200 #us
QZWINDOW=55 #cm
QYWINDOW=12 #cm
PROCESS_OFFSET=4372

NJOBS=10
#if [ ! -e "$SINGLEHIT_FILELIST" ]; then
#   echo "Singlehit filelist does not exist: ${SINGLEHIT_FILELIST}"
#   exit 1
#fi

#if [ ! -e "$RECO_FILELIST" ]; then
#   echo "Input filelist does not exist: ${RECO_FILELIST}"
#   exit 1
#fi
echo "submitting job..."
jobsub_submit -N ${NJOBS} \
--email-to fogar314@colostate.edu \
--memory=2000MB \
--expected-lifetime=1h \
--singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest \
-G dune -e IFDH_CP_MAXRETRIES=2 \
--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE -e UPS_OVERRIDE="-H Linux64bit+3.10-2.17" \
file://${USERDIR}/lowenergy_job.sh $RUN $SINGLEHIT_FILELIST $RECO_FILELIST $OUTDIR $CHANNELMAP $OPTAG $OPTWINDOW $OPZWINDOW $OPYWINDOW $MATCHTWINDOWLOW $MATCHTWINDOWHIGH $QZWINDOW $QYWINDOW $PROCESS_OFFSET