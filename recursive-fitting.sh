#!/bin/bash

INTRJ=input/out.trj
INTOP=input/topology.gro
OUTDIR=output

./fitting.py -t $INTRJ -p $INTOP -o $OUTDIR/fitted0.trj -i
./fitting.py -t $OUTDIR/fitted0.trj -p $INTOP -o $OUTDIR/fitted1.trj
./fitting.py -t $OUTDIR/fitted1.trj -p $INTOP -o $OUTDIR/fitted2.trj
./fitting.py -t $OUTDIR/fitted2.trj -p $INTOP -o $OUTDIR/fitted3.trj
./fitting.py -t $OUTDIR/fitted3.trj -p $INTOP -o $OUTDIR/fitted4.trj
./fitting.py -t $OUTDIR/fitted4.trj -p $INTOP -o $OUTDIR/fitted5.trj