USERDIR=/exp/dune/data/users/sfogarty/DUNEFD_39Ar/JOBSUB/

RUN=28850
SINGLEHIT_FILELIST=${USERDIR}/singlehit_lists/singlehit_${RUN}_filepaths.txt
RECO_FILELIST=${USERDIR}/runlists/${RUN}_filepaths.txt
#OUTDIR=/pnfs/dune/scratch/users/sfogarty/CLMatching/${RUN}
OUTDIR=${USERDIR}/CLMatching/${RUN}
mkdir -p $OUTDIR
CHANNELMAP=/pnfs/dune/scratch/users/sfogarty/PDHD_PDS_ChannelMap.csv
OPTAG=pdhddaphne:daq
OPTWINDOW=0.1 #us
OPZWINDOW=600 #cm
OPYWINDOW=600 #cm
MATCHTWINDOWLOW=0 #us
MATCHTWINDOWHIGH=2200 #us
QZWINDOW=55 #cm
QYWINDOW=12 #cm
PROCESS_OFFSET=0

if [ ! -e "$SINGLEHIT_FILELIST" ]; then
   echo "Singlehit filelist does not exist: ${SINGLEHIT_FILELIST}"
   exit 1
fi

if [ ! -e "$RECO_FILELIST" ]; then
   echo "Input filelist does not exist: ${RECO_FILELIST}"
   exit 1
fi
PROCESS=0
OUTFILE=${OUTDIR}/PDHD_CLMATCHING_RESULTS_${RUN}_jobnum${PROCESS}.root
#kx509;voms-proxy-init -noregen -rfc -voms dune:/dune/Role=Analysis
#LD_PRELOAD=$XROOTD_LIB/libXrdPosixPreload.so root -l 'CL_lowenergy_analysis.C(${RUN}, ${PROCESS}, ${SINGLEHIT_FILELIST}, ${RECO_FILELIST}, ${OUTFILE}, ${CHANNELMAP}, ${OPTAG}, ${OPTWINDOW}, ${OPZWINDOW}, ${OPYWINDOW}, ${MATCHTWINDOWLOW}, ${MATCHTWINDOWHIGH}, ${QZWINDOW}, ${QYWINDOW})'
#LD_PRELOAD=$XROOTD_LIB/libXrdPosixPreload.so root -l "CL_lowenergy_analysis.C(${RUN}, ${PROCESS}, \"${SINGLEHIT_FILELIST}\", \"${RECO_FILELIST}\", \"${OUTFILE}\", \"${CHANNELMAP}\", \"${OPTAG}\", \"${OPTWINDOW}\", \"${OPZWINDOW}\", \"${OPYWINDOW}\", \"${MATCHTWINDOWLOW}\", \"${MATCHTWINDOWHIGH}\", \"${QZWINDOW}\", \"${QYWINDOW}\")"
LD_PRELOAD=$XROOTD_LIB/libXrdPosixPreload.so root -l -q -b "CL_lowenergy_analysis.C(${RUN}, ${PROCESS}, \"${SINGLEHIT_FILELIST}\", \"${RECO_FILELIST}\", \"${OUTFILE}\", \"${CHANNELMAP}\", \"${OPTAG}\", \"${OPTWINDOW}\", \"${OPZWINDOW}\", \"${OPYWINDOW}\", \"${MATCHTWINDOWLOW}\", \"${MATCHTWINDOWHIGH}\", \"${QZWINDOW}\", \"${QYWINDOW}\")"
