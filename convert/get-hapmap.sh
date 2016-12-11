mkdir -p hapmap
pushd hapmap
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz -O hapmap.tar.gz
tar -xvf hapmap.tar.gz
popd
