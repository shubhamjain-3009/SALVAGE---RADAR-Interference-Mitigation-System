# SALVAGE---RADAR-Interference-Mitigation-System

To run the code, you first need to run the scenario simulator called: FMCW_Rev_11.m. This will generate a file called temp_file1.mat. This file has the results from simulation
You now need to run a file called  FMCW_Rev10_Analysis.m. This file leads the temp_file1.mat and implements the receiver DSP. 
This includes the proposed SALVAGE Algorithm. To bypass the algorithm, you need to set ignore_NDMC_flag = 1.
You may want to change the NDMC_threshold variable on top. It triggers the mechanism and should be set at a value (in dBm) which is just slightly above a no interference case. Typical value has been set already. 

To simulate a scenario, you will need to set values of distance, velocity, angle and slopes for the ego, targets and interferers. The variables have to be set in the FMCW_Rev_11.m file. These are as follows: 
1. d_int: array containing distances of interferers. Preferably values between 5 and 100.
Note that every interferer is also a target here. If you want to turn off the effect of interference, you can either leave this array blank, or set int_present to 0. 

2. theta_int, v_int, K_int: As above, you can set the values in these variables for velocity, angle and slope of the interferers. Keep the same number of elements in these as you had int he distance vector. keep theta between -90 ad 90, and preferably between -40 and 40. keep slopes around 10-40 MHz/us. 

3. d_tar, v_tar, theta_tar: Just as above, these are the distances angles and velocities of the targets in the scene. Similar constraints apply

4. N_f: Number of frames to run the simulation for. TYpically 5 frames gives a good picture overall. 

5. Chirps_per_frame = Number of Chirps per frame. Typically 50. And its a value I dont change much. 


Addtional Info: Following functions were used for:
1. Plot_Basic_FFT: It operates on two .met files which were obtained with and without 1 interferer. It can be used to see the effect of different k interference on inndividual chirps by setting the chp number on top between 1 and 50. It also shows the comparision between excision and no excision in time domain. To modify that, you can change the indices 96 and 490 in line 29. These are the start and end indices of sample length you want to excise. Change chirp, look at time domain, figure out what youw ant to excise, and enter those values in line 29. 

2. LMS: It operates on the temp_file1.mat. You can choose the chirp you  want to  examine by setting the frm and chp parameters(frame and chirp). You can play aroudn with lenth of filter l, step size for MATLAB filter mu, and my own implementation step (line 22)

Following functions are critical: 
1. Create_Time_Matrix and Create_Time_Matrix_Int: Work to give you the exact times at which simulations are carried out
2. Update_Scene_Mod:  Generates d, theta for every sim time
3. Offset Chirp and Int offset chirp: Imlpement samlpe delays which emulate Time of flight. 



NOTE: Apologies for the efforts it might take to go through my code. Would have been more thourough with documentation but ran out of time. 
