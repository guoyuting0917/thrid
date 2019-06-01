#!/bin/sh

echo -e "remove output file ? (y/[n]): \c"
read ISREMOVE

if [ "$ISREMOVE" = "y" ] || [ "$ISREMOVE" = "yes" ]; then
    rm  out_sum.dat
    rm  out_ene.dat
    rm  out_pos.dat
    rm  out_vel.dat
    if [ -e "out_for.dat" ]; then
        rm  out_for.dat
    fi
    rm  out_the.dat
    rm  out_bar.dat
    rm  out_pre.dat
    rm  out_pdb.pdb
    if [ -e "out_thc.dat" ]; then
        rm  out_thc.dat
    fi
    if [ -e "out_htf.dat" ]; then
        rm  out_htf.dat
    fi
    if [ -e "out_mtf.dat" ]; then
        rm  out_mtf.dat
    fi
    if [ -e "out_umb.dat" ]; then
        rm  out_umb.dat
    fi
fi

exit 0
