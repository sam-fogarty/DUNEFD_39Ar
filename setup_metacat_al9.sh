source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load metacat@4.0.0
export METACAT_AUTH_SERVER_URL=https://metacat.fnal.gov:8143/auth/dune
export METACAT_SERVER_URL=https://metacat.fnal.gov:9443/dune_meta_prod/app 
metacat auth login -m password sfogarty
