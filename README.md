Params:
- initialization functions: 0, 1, 2, 3, 4, 5 --> loading from existing file, sinusoid, rain drop, multiple rain drops, square, velocity
- keep adding new data during wave propaation, e.g., random rain drops: true / false
- iterative compression: true/false --> lossy compress the data after each update and propagate the erroneous data to future steps
- cdf condition: 1, 2, 3 --> Dirichlet, Mur absorption, Nueman (only 2D)
- Dx: x-axis domain size
- Dy: y-axis domain size
- Dz: z-axis domain size

3D wave equation:
writing out the uncompressed data: ./waveProp_3d 2 false false 2 20 20 20 0.1 0.2 0.2 80 0.0 ABS fiveRainDrop_20_20_20_h01_t02_Mur.bp fiveRainDrop_20_20_20_h01_t02.bp 0 
