#!/bin/bash

INTRJ=workspace/01-centering/centering.npz
OUTDIR=workspace/02-fitting

./fitting.py -n $INTRJ -o $OUTDIR/fitted0 -i
./fitting.py -n $OUTDIR/fitted0.npz -o $OUTDIR/fitted1
./fitting.py -n $OUTDIR/fitted1.npz  -o $OUTDIR/fitted2
./fitting.py -n $OUTDIR/fitted2.npz -o $OUTDIR/fitted3
./fitting.py -n $OUTDIR/fitted3.npz -o $OUTDIR/fitted4
./fitting.py -n $OUTDIR/fitted4.npz -o $OUTDIR/fitted5