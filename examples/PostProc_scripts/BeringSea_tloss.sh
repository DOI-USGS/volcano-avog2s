#!/bin/bash

rc=0

echo "Cleveland:-169.945:52.822:1.73" > Volcs.txt

nvolcs=`wc -l Volcs.txt | cut -c1-2`
nv=1
Volc=`head -n ${nv} Volcs.txt | tail -1 | cut -f1 -d':'`

# We need to know if we must prefix all gmt commands with 'gmt', as required by version 5
GMTv=5
type gmt >/dev/null 2>&1 || { echo >&2 "Command 'gmt' not found.  Assuming GMTv4."; GMTv=4;}
GMTpre=("-" "-" "-" "-" " "   "gmt ")
GMTelp=("-" "-" "-" "-" "ELLIPSOID" "PROJ_ELLIPSOID")
GMTnan=("-" "-" "-" "-" "-Ts" "-Q")
GMTrgr=("-" "-" "-" "-" "grdreformat" "grdconvert")

mapscale="6i"

lonw=183.0
lats=51.0
lone=202.0
latn=59.5

PROJp=-JM195.5/55.0/${mapscale}

AREAp="-R${lonw}/${lone}/${lats}/${latn}"
DETAILp=-Dh
COASTp="-G220/220/220 -W"
COASTp="-Ggrey90 -W -Slightsteelblue1"

zmin=-180.0
zmax=-80.0
dz=2.0
cpt="GMT_hotinv.cpt"
echo "0.0   255     255     255 0.375 255     255     0"   > ${cpt}
echo "0.375 255     255     0   0.75  255       0     0"  >> ${cpt}
echo "0.75  255     0       0   1.0     0       0     0"  >> ${cpt}
makecpt -C${cpt} -T${zmin}/${zmax}/${dz} > tl.cpt

#${GMTpre[GMTv]} pscoast $AREAp $PROJp $DETAILp -Sc -K > temp.ps
# Filter results to strip out 'inf', blank lines and the header
#grep '[0-9]' Bering_results.dat | grep -v '#' | grep -v inf > results.dat
# Plot bounce points colored by travel-time
#awk '{print $2, $1, $3}' Nby2D_tloss_1d2.llt | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O -Sc2.5p -Ctl.cpt  >> temp.ps
#awk '{print $5, $4, $6}' results.dat | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O -Sc1.5p -Ctt.cpt  >> temp.ps
${GMTpre[GMTv]} ${GMTrgr[GMTv]} "tloss.nc?tloss" tloss.grd
${GMTpre[GMTv]} grdimage ${GMTnan[GMTv]} tloss.grd -Ctl.cpt $AREAp $PROJp -K  > temp.ps
${GMTpre[GMTv]} pscoast $AREAp $PROJp $DETAILp -Ggrey90 -K -O >> temp.ps

# Add colorbar

# Plot all the infrasound arrays
Clevx=-169.94
Clevy=52.822
SndPtx=-160.49
SndPty=55.337
Aktx=-165.99
Akty=54.133
Okx=-168.1750
Oky=53.468
Dilx=-158.51
Dily=59.047
Adkx=-176.6581
Adky=51.88
echo "${Clevx} ${Clevy} 1.0"   | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Ss10.0p -Gblack -W0.5,0/0/0  >> temp.ps
echo "${SndPtx} ${SndPty} 1.0" | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Ss10.0p -Gblack -W0.5,0/0/0  >> temp.ps
echo "${Aktx} ${Akty} 1.0"     | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Ss10.0p -Gblack -W0.5,0/0/0  >> temp.ps
echo "${Okx} ${Oky} 1.0"       | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Ss10.0p -Gblack -W0.5,0/0/0  >> temp.ps
echo "${Dilx} ${Dily} 1.0"     | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Ss10.0p -Gblack -W0.5,0/0/0  >> temp.ps
echo "${Adkx} ${Adky} 1.0"     | ${GMTpre[GMTv]} psxy $AREAp $PROJp -K -O  -Ss10.0p -Gblack -W0.5,0/0/0  >> temp.ps

psscale -D2.75i/-0.4i/4i/0.15ih -Ctl.cpt -B20f20/:"TLoss (Db)": -O -K >> temp.ps

#${GMTpre[GMTv]} psbasemap -B5g5:."Strat-Modess-${Volc} ": $AREAp $PROJp -O >> temp.ps
${GMTpre[GMTv]} psbasemap -B5g5:."Strat-WMod-${Volc} ": $AREAp $PROJp -O >> temp.ps
ps2epsi temp.ps
epstopdf temp.epsi
mv temp.pdf ${Volc}_tloss.pdf
convert -density 400 temp.epsi -rotate 90 -resize x750 -background white -flatten ${Volc}_tloss.png

