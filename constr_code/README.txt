There's a Makefile, a parameterfile and two particle input files "h40.txt" and "h40_b8.txt".

---------------------------
Compilation: you can disable hdf5 by commenting out the HDF5_FLAGS line (it's only required if you want an hdf5 output file, right now it's set to produce only the gadget file -> see parameterfile)

---------------------------
Basic stuff:
Once compiled, you can run the code with

./IC_deriv paramfile.txt 0 0.3

where the first number is the number of iterations (-1 skips that step entirely), and the second number is the constraint number "d". This should take about 1 minute to run. (5 iterations take , 10 take 30 mins at 256^3 resolution).

You may want to change the "outdir" in the paramfile.

I have also included a file "sample_output.txt" with the stuff that the code writes to the shell (for the settings used above).

---------------------------
Code structure:
This is the full IC generator, and the main function is in desperate need to be broken up into smaller pieces ;-)

Lines 1370 - 1567 generate the initial density field (white noise + PS + normalisations)

1568 - 1625: alpha_mu calculation

1672 - 1718: Iteration

rest: Zeldovich, output

---------------------------
Important functions:

cen_deriv4_alpha: finite difference to fourth order

calcA : recursion/iteration over alpha and unconstrained field

Grid class: sets up a 3D grid and allows to find new indices when moving along the directions for the finite difference (this could probably be used at other places in the code, but didn't bother changing bits that work)

---------------------------
Advanced stuff:
If you don't change anything, the particle input file will be read and the constraints will apply to these particles. If you want to run at lower resolution you will have to make some changes around line 1589: adjust n_in_bin, comment out the GetBuffer... and uncomment the part_arr[i]=i in the loop. Then it will constrain the first n_in_bin particles, starting from the grid corner (unless you center, see below).

If you want to use the particles in bin 8 you have to change the IDfile in the parameterfile and change n_in_bin to 2215 in line 1590.

You can also play around with the "vector" which is the direction of AM to be constrained. You can choose any values there, it will be normalised within the code.

Currently there is no "centering" done in the constraint calculation. You can enable this by uncommenting line 1511.

You can output the angular momentum for each particle by enabling the ofs1/2/3 things around line 1600.
