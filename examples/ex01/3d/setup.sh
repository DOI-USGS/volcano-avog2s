ln -s ../../../bin/g2s_Extract_Grid .
ln -s ../G2S*raw .
ln -s ~/work/USGS/Ground2Space/GeoAc-v1.1.1/GeoAcGlobal.RngDep .

echo "Now you can generate the profiles you need by running:"
echo "./g2s_Extract_Grid example1_ext3d.ctr "
echo " "
echo "And then run the foreward 3-d range-dependent model using:"
echo "./GeoAcGlobal.RngDep -prop Clev Clev.loclat Clev.loclon theta_min=0.0 theta_max=45.0 theta_step=0.5 \\"
echo "   phi_min=0.0 phi_max=30.0 phi_step=10.0 bounces=10 lat_src=52.8222 lon_src=-169.945 z_src=1.73 CalcAmp=False"
echo " "
echo "Finally, plot the results using plot_3d in octave or matlab"
