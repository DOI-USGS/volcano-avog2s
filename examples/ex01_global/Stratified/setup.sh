ln -s ../G2S*raw .
ln -s ../../PostProc_scripts/load_GeoAc_atmo.m .
ln -s ../../PostProc_scripts/plot_GeoAc_profile.m .
ln -s ../../PostProc_scripts/load_StationInfo.m .
ln -s ../../PostProc_scripts/load_GeoAc_raypaths.m .


echo "Now you can generate the profile you need (e.g. for Cleveland Volcano) by running:"
echo "../g2s_Extract_Sonde example1_step3_ext1d.ctr"
echo " "
echo "-------------------------------------------------"
echo "   GeoAc "
echo "      https://github.com/LANL-Seismoacoustics/GeoAc"
echo "-------------------------------------------------"
echo "For the forward 1-d stratified GeoAc model, use:"
echo "GeoAcGlobal -prop Clev0.met theta_min=-30.0 theta_max=55.0 theta_step=1.0 azimuth=41 bounces=10 lat_src=52.8222 lon_src=-169.945 z_src=1.73 CalcAmp=False WriteAtmo=True"


echo "-------------------------------------------------"
echo "   NCPA codes for a stratified atmosphere "
echo "      https://github.com/chetzer-ncpa/ncpaprop"
echo "-------------------------------------------------"
echo "For the stratified Modess (Model Effective Sound Speed), use"
echo "  One profile:"
echo "    Modess --atmosfile Clev0.met --skiplines 0 --atmosfileorder ztuvdp --azimuth 41 --freq 0.1 --write_2D_TLoss --sourceheight_km 1.73"
echo "  Sweep of profiles:"
echo "    Modess --atmosfile Clev0.met --skiplines 0 --atmosfileorder ztuvdp --freq 0.1 --Nby2Dprop --azimuth_start 0 --azimuth_end 360 --azimuth_step 1 --sourceheight_km 1.73"
echo ""
echo "For the Wide-Angle High-Mach Modal Code, use"
echo "  One profile:"
echo "    WMod --atmosfile Clev0.met --skiplines 0 --atmosfileorder ztuvdp --azimuth 41 --freq 0.1 --write_2D_TLoss --sourceheight_km 1.73"
echo "  Sweep of profiles:"
echo "    WMod --atmosfile Clev0.met --skiplines 0 --atmosfileorder ztuvdp --freq 0.1 --Nby2Dprop --azimuth_start 0 --azimuth_end 360 --azimuth_step 1 --sourceheight_km 1.73"
echo "For the Complex Effective Sound Speed (single tone normal mode expansion with attenuation via complex eigenvalues), use"
echo "  One profile:"
echo "    CModess --atmosfile Clev0.met --skiplines 0 --atmosfileorder ztuvdp --azimuth 41 --freq 0.1 --write_2D_TLoss --sourceheight_km 1.73"

