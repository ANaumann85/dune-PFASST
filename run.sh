bindir="./debug/src/examples/FE_heat2d"
mpirun -np 2 ${bindir}/FE_heat2d_pfasst color=false abs_res_tol=10e-17 --nFine=1 --nCoarse=1 --num_elements=10 tend=20.0 --num_iters=10 --num_steps=8 > log
#${bindir}/FE_heat2d_mlsdc color=false abs_res_tol=10e-17 --nFine=2 --nCoarse=2 --num_elements=10 --tend=20.0 --num_iters=4 --num_steps=8 > log
#gdb --args ${bindir}/FE_heat2d_mlsdc color=false abs_res_tol=10e-17 --nFine=1 --nCoarse=1 --num_elements=10 --tend=20.0 --num_iters=4 --num_steps=8 
#${bindir}/FE_heat2d_sdc color=false abs_res_tol=10e-17 --num_elements=20 --tend=20.0 --num_iters=10 --num_steps=8 > log
