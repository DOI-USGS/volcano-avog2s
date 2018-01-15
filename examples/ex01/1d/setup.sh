ln -s ../../../bin/g2s_Extract_Grid .
ln -s ../G2S*raw .
ln -s ~/work/USGS/Ground2Space/GeoAc-v1.1.1/GeoAc2D .

echo "Now you can generate the profiles you need by running:"
echo "./g2s_Extract_Grid example1_ext3d.ctr "
echo " "
echo "And then run the foreward 1-d stratified model using:"
echo "./GeoAc2D -prop Clev0.met theta_min=0.0 theta_max=45.0 theta_step=0.5 azimuth=41 bounces=10 z_src=1.73 CalcAmp=False"
echo " "
echo "Finally, plot the results using plot_1d in octave or matlab"

