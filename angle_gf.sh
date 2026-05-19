#!/bin/bash


x_side="170.0"
y_side="70.0"
z_side="40.0"
n_count="100000000"

for num in {0..359} ; do
    if [ $num -lt 10 ]
        then
        ext="_deg00$num"
    elif [ $num -lt 100 ]
        then
        ext="_deg0$num"
    else
        ext="_deg$num"
    fi

printf "2\n" > angle_in.txt
printf "$x_side $y_side 0.0\n" >> angle_in.txt
printf "$x_side $y_side $z_side\n" >> angle_in.txt
printf "$n_count\n" >> angle_in.txt
printf "$num.0\n" >> angle_in.txt
printf "$ext\n" >> angle_in.txt

head angle_in.txt

./geofac_2d_WAL_2026_norotation

./geom_fac_plot_2d_WAL_2_orientation $ext

mv plots/dgf_2D/*$ext* plots/dgf_2D/angle/.

mv *$ext.txt angle/.

done
