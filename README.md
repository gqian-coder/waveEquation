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

Example: 
- 3D wave equation, multiple rain drops, Mur condition, dissipation rate of 5e-3, writing out the uncompressed data:
    ./waveProp_3d 3 false false 2 20 20 20 0.1 0.2 80 0.2 5e-3 0.0 ABS fiveRainDrop_20_20_20_h01_t02_Mur.bp nan 0 
- previous case, with compression:
    ./waveProp_3d 3 false false 2 20 20 20 0.1 0.2 80 0.2 5e-3 1e-3 ABS fiveRainDrop_20_20_20_h01_t02_Mur_mgr_abs_1e-3.bp fiveRainDrop_20_20_20_h01_t02_Mur.bp 0 
- checkpoint restart @ the 1000th timestep, propagating the erronous data to 3000 timesteps and writing out the timestep results:
    ./waveProp_3d 3 false false 2 20 20 20 0.1 0.2 60 0.2 1e-3 0.0 ABS fiveRainDrop_20_20_20_h01_t02_Mur_mgr_abs_1e-3_start_1000ts.bp fiveRainDrop_20_20_20_h01_t02_Mur_mgr_abs_1e-3.bp 1000 
