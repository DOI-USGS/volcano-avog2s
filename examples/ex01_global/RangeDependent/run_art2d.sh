#!/bin/bash

ART2DHOME="${HOME}/bin"

tmpfile="tmp.ctr"

envfile="InfraAtmos01.env"
infile="ClevDLL.in"
outfile="art2d_xsec.bin"
topofile="topo.in"

echo "Preparing files from Art2d"
# First build control file for g2s_env_convert
echo "$envfile" > $tmpfile
echo "$infile" >> $tmpfile
echo "$topofile" >> $tmpfile
${ART2DHOME}/g2s_env_convert < $tmpfile

echo "Running Art2d"
# First build control file art2d
echo "$infile" > $tmpfile          # sound speed file
echo "0 5000 0 140" >> $tmpfile        # x/z min/max [km]
echo "6371" >> $tmpfile                # Earth rad., flat earth transform
echo "$topofile" >> $tmpfile                  # topo file
echo "2500 1 0" >> $tmpfile             # x0,z0,t0 [km,km,s]
echo "0.1" >> $tmpfile                 # step size [km^2/s]
echo "3" >> $tmpfile                   # type of output (1,2,3)
echo "1 179 179" >> $tmpfile             # elev1,elev2,nray
echo "1" >> $tmpfile                   # irflct (1 or 0)
echo "4" >> $tmpfile                   # iord (order of RK; 2 or 4)
${ART2DHOME}/art2d < $tmpfile 
mv -f art2d.bin $outfile
# Now split off the first line of the infile
head -1 $infile > dims.txt
tail -n +2 < $infile > art2d_data.in
##./plot_art2d



