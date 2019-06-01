#!/bin/sh

### set src directory of peachgk_md
SRCDIR=$HOME/peachgk_md/peachgk-md-2.142

echo 'install executable and input file'

cp -i $SRCDIR/param/C7H15OH_cor.dat .
cp -i $SRCDIR/param/C7H15OH_top.dat .
cp -i $SRCDIR/param/C5H11OH_cor.dat .
cp -i $SRCDIR/param/C5H11OH_top.dat .
#cp -i $SRCDIR/param/H2O_top.dat .
#cp -i $SRCDIR/param/H2O_topOE.dat .
cp -i $SRCDIR/param/H2O_topOB.dat .
cp -i $SRCDIR/param/para_bond.dat .
cp -i $SRCDIR/param/para_vdw_s.dat .
cp -i $SRCDIR/param/para_const.dat .
#cp -i $SRCDIR/param/add_top.dat .
#cp -i $SRCDIR/param/para_cstmnb.dat .

cp -i $SRCDIR/peachgk_md.out .

cp -i $SRCDIR/peachgk.ini .
if [ "$1" = "-tf" ]; then
  cp -i $SRCDIR/transflux.ini .
fi
if [ "$1" = "-tc" ]; then
  cp -i $SRCDIR/tempcont.ini .
fi
if [ "$1" = "-hfc" ]; then
  cp -i $SRCDIR/hfcont.ini .
fi
if [ "$1" = "-lt" ]; then
  cp -i $SRCDIR/langecont.ini .
fi
if [ "$1" = "-umb" ]; then
  cp -i $SRCDIR/potbias.ini .
fi

cp -i $SRCDIR/rmout.sh .

exit 0
