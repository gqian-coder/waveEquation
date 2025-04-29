#!/bin/sh

# step 1: simulating wave equation and outputting the uncompressed data starting from step 2000 
# h=1 and delta_t = 1e-3, velocity condition was initialized from a mask --- CurveFault_A_t20_512_512_512_vmax_68.bin
./waveProp_twoStep_3d 6 0 3 512 512 512 1 0.001 4 0 CurveFault_A_t20_512_512_512_vmax_68.bin mask 0 0 0 ABS gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68.bp ./ 0 1 1 0 2000

# step 2: reading the wave field data at step 2000, compressing and saving it as checkpoint data using zfp 
# eb for l-inf error = 4e-3, and the achieved l2 error is around 1e-3 
./waveProp_twoStep_3d_zfp 0 0 3 512 512 512 1 1e-3 0.01 0 CurveFault_A_t20_512_512_512_vmax_68.bin mask 0 4e-3 gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68_zfp_eb_4e-3_2000ts.bp gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68.bp 0 1 0 2000 0

# step 3: loading the compressed checkpoint, restarting and continuing the simulation for another 2000 steps 
./waveProp_twoStep_3d_zfp 0 0 3 512 512 512 1 1e-3 2 0 CurveFault_A_t20_512_512_512_vmax_68.bin mask 0 0 gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68_eb_4e-3_zfp_cr2000ts.bp gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68_eb_4e-3_zfp_2000ts.bp 0 1 0 2000 0

# step 4: error evaluation 
./waveEnergy_3d gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68.bp gaussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68_eb_4e-3_zfp_cr2000ts.bp CurveFault_A_t20_512_512_512_vmax_68.bin 1 1e-3 1 gussian_pulse_src_wave_512_512_512_h1_t001_CurveFault-A_noSourceAfter025s_vmax_68_eb_4e-3_zfp_cr2000ts_nonOpt 0 1
