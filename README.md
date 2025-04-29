Params:
- initialization functions: 0, 1, 2, 3, 4, 5 --> loading from existing file, sinusoid, rain drop, multiple rain drops, square, velocity
- keep adding new data during wave propaation, e.g., random rain drops: true | false
- iterative compression: true | false --> lossy compress the data after each update and propagate the erroneous data to future steps
- cdf condition: 1, 2, 3 --> Dirichlet, Mur absorption, Nueman (only 2D)
- Dx: x-axis domain size
- Dy: y-axis domain size
- Dz: z-axis domain size
- dh: simulation spatial resolution
- dt: simulation temporal resolution
- T: simulation duration
- C: wave speed
- gamma: wave energy dissipation rate
- tol: tolerance used for compression the timestep data (tol==0 means no compression)
- eb_type: ABS | REL --> absolute or relative error bound
- fname_wt: output file containing the timestep update data
- fname: pre-existing file for loading the checkpoint data
- init_ts: checkpoint restart timestep 

Scripts to reproduce the wave equation simulation and the C/R with lossy compression: 
1) Compression with mgard energy error control
    mgard_energy_eval.sh
2) Compression with regular mgard l2 error control
    mgard_l2_eval.sh
3) Compression with sz3 l-inf error control
    sz_eval.sh
5) Compression with zfp l-inf error control
    zfp_eval.sh
