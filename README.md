solving poisson problem
	$-\Delta u = f$,  in $[0, 1]^3$
	$u = u_\partial$, on $\partial[0, 1]^3$

require package AFEPack and Deal.II

./run file_name_of_mesh M file_name_of_output_solution:
	solving poisson problem with polynomial order M, 
	exact solution and right-hand-side are defined insider poisson.cpp

./mesh_generation/main a nx ny nz: 
	divive $[0, a]^3$ with equal partition number nx, ny and nz for each direction

in folder ./visualize
    ./run file_name_of_mesh file_name_of_solution N_local_refine tol:
	    generate opendx file for visualization
