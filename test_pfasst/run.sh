bindir="../../release/src/examples/FE_heat2d/"

nSteps="2 4 8 16 32 64 128 256"
nFines="1 2 4"
nCoarses="1 2 4"
nIters="1 2 3 4 8"
nIters="3"
nCores="2 4"
#nCores="4"

binFix="${bindir}/FE_heat2d_pfasst color=false abs_res_tol=10e-17 --num_elements=10 tend=20.0 --nCoarse=1 --nFine=1"
for nStep in $nSteps; do
  for nCore in $nCores; do
    for nIter in $nIters; do
      clog="log_${nIter}_${nStep}_${nCore}"
      echo "$clog"
      nStep2=$(($nStep*$nCore))
      mpirun -np $nCore ${binFix} --num_iters=$nIter --num_steps=$nStep2 > $clog
    done
  done
done
