ln -s ../G2S*raw .
ln -s ../../PostProc_scripts/load_GeoAc_atmo.m .
ln -s ../../PostProc_scripts/plot_GeoAc_profile.m .
ln -s ../../PostProc_scripts/load_StationInfo.m .
ln -s ../../PostProc_scripts/load_GeoAc_raypaths.m .



echo "Now you can generate the profile you need (e.g. for Cleveland Volcano) by running:"
echo "../g2s_Extract_Sonde example2_step3_ext1d_FC.ctr"
echo " "
echo "-------------------------------------------------"
echo "   GeoAc "
echo "      https://github.com/LANL-Seismoacoustics/GeoAc"
echo "-------------------------------------------------"
echo "For the forward 1-d stratified GeoAc model, use:"
echo "GeoAc2D -prop Clev0.met theta_min=-30.0 theta_max=55.0 theta_step=1.0 azimuth=41 bounces=10 z_src=1.73 CalcAmp=False"
echo " "
echo "Finally, plot the results using plot_GeoAc_profile.m in octave or matlab"
echo " "

