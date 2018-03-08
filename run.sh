bindir="./debug/src/examples/FE_heat2d"
#mpirun -np 2 ${bindir}/FE_heat2d_pfasst abs_res_tol=10e-17 --num_elements=10 tend=20.0 --num_iters=10 num_steps=8 > log
${bindir}/FE_heat2d_mlsdc color=false abs_res_tol=10e-17 --num_elements=10 --tend=20.0 --num_iters=10 --num_steps=8 > log
#${bindir}/FE_heat2d_sdc color=false abs_res_tol=10e-17 --num_elements=10 --tend=20.0 --num_iters=10 --num_steps=8 > log
