PDIR=$PWD

WDIR=/Users/Hans/SCRATCH/SPOTDIAGRAMS ; [[ ! -d $WDIR ]] && mkdir -p $WDIR
DDIR=/Users/Hans/SPOTDIAGRAMS/DATA

nbeam=50000
nbeam=100000
colours='r g b d v'
design=$DDIR/sphere.dat
design=$DDIR/newton.dat
design=$DDIR/apoklaas.dat
chmod 444 $design

FC="gfortran -fdefault-integer-8  -fdefault-real-8"

FLIST="spotmod.F90 spotio.F90 spot.F90"
prog=spot.x

cd $PDIR ; cp $FLIST $WDIR

cd $WDIR
   [[ -f $prog ]] && rm -f $prog
   $FC -o $prog $FLIST

for c in $colours ; do
   fout=$WDIR/$(basename $design | awk -F.dat '{print $1}').$c
   ./$prog -i $design -o $fout -n $nbeam -c $c

   col='black'
   [[ $c == 'r' ]] && col='red'
   [[ $c == 'g' ]] && col='green'
   [[ $c == 'b' ]] && col='blue'
   [[ $c == 'd' ]] && col='yellow'
   [[ $c == 'v' ]] && col='violet'
   gnuplot -e "plot '$fout' with dots lc rgb '${col}'"
done
