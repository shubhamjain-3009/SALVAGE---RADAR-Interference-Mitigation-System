# SALVAGE---RADAR-Interference-Mitigation-System
Basic Steps to Run code:

To run the code, you first need to run the scenario simulator called: FMCW_Simulate.m. This will generate a file called temp_file1.mat. This file has the results from simulation
You now need to run a file called  FMCW_Analyze.m. This file reads the temp_file1.mat and implements the receiver DSP. This includes the proposed SALVAGE Algorithm. 

Setting the Scenario:
To simulate a scenario, you will need to set values of distance, velocity, angle and slopes for the ego, targets and interferers. The variables have to be set in the FMCW_Simulate.m file. These are as follows: 

1. d_int: array containing distances of interferers. Preferably values between 5 and 100.
Note that every interferer is also a target here. If you want to turn off the effect of interference, you can either leave this array blank, or set int_present to 0. 

2. theta_int, v_int, K_int: As above, you can set the values in these variables for velocity, angle and slope of the interferers. Keep the same number of elements in these as you had int he distance vector. keep theta between -90 ad 90, and preferably between -40 and 40. keep slopes around 10-40 MHz/us. 

3. d_tar, v_tar, theta_tar: Just as above, these are the distances angles and velocities of the targets in the scene. Similar constraints apply

4. N_f: Number of frames to run the simulation for. TYpically 5 frames gives a good picture overall. 

5. Chirps_per_frame = Number of Chirps per frame. Typically 50. And its a value I dont change much. 


Reproducing Results:
Case 1: No Mitigation 
Set variable Inter_Frame_Gap to case 1 value in Simulation code (FMCW_Simulate.m). Set Ignore_NDMC_flag to 1 in Analysis Code (FMCW_Rev10_Analysis.m)

Case 2: Time Dithering 
Set variable Inter_Frame_Gap to case 3 value in Simulation code (FMCW_Simulate.m). Set Ignore_NDMC_flag to 1 in Analysis Code (FMCW_Rev10_Analysis.m)

Case 3: SALVAGE 
Set variable Inter_Frame_Gap to case 3 value in Simulation code (FMCW_Simulate.m). Set Ignore_NDMC_flag to 0 in Analysis Code (FMCW_Rev10_Analysis.m)

Note: You may want to change the NDMC_threshold variable in the analysis. It triggers the mechanism and should be set at a value (in dBm) which is just slightly above a no interference case. Typical value has been set already. 

Addtional Info: Following functions were used for:
1. Plot_Basic_FFT: It operates on the temp_file1.mat and can be used to see the effect of different k interference on individual chirps by setting the chp number on top between 1 and 50. It also shows the comparision between excision and no excision in time domain. To modify that, you can change the start and stop index in Excision section of code. These are the start and end indices of sample length you want to excise. Change chirp, look at time domain, figure out what samples you want to excise, and enter those values in Excision section. 

2. LMS: It operates on the temp_file1.mat and attempts to use LMS filtering to clean up a chirp. You can choose the chirp you want to examine by setting the frm and chp parameters(frame and chirp). You can play around with lenth of filter l, step size for MATLAB filter mu, and my own implementation step size with the variable "step".
