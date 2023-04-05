clear all;
close all;

%%%% This code generates a scenario and creates RADAR's receive Waveform (Ego) based on
%%%% that. This waveform consists of reflection from targets, interferering
%%%% RADAR waveforms. It is captured per chirp, for a certain number of
%%%% chirps per frame (Chirps_Per_Frame) and for a certain number of frames (N_f)


%%%%%%%%% Defining Simulation Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_f = 5; %%% Number of frame for which the simulator needs to be run


%%%%% Defining Scenario %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_int = [18 41 50 91 98]; %%% distances of interferers from Ego Vehicle
N_int = length(d_int);
v_int = [20 40 50 90 100]; %%% velocities of interferers
rcs_int = ones(1,N_int); %%% rcs of interferers
theta_int = [-32 -11 3 17 31]; %%% AoA of interferers

int_present = 1; %1 = Present, 0 = Not 
% This control allows you to simulate every interfering vehicle as
% a target only. Kind of like having the interfering vehicles at
% the positions you specify above, but not having transmitting
% RADARs on any of them.

K_inter = [10e12 17e12 23e12 12e12 30e12]; %%% FMCW Sweep slopes of interferers
f_start_int = 8e09*ones(1,N_int); %%% Start frequencies of interferers
Inter_Frame_Gap_Int_val = 20e-06*ones(1,N_int); %%% Inter Frame Gap of interferers
Inter_Frame_Gap_Int = zeros(N_int,N_f); 
Tx_Start_Time_Int = 0*ones(1,N_int);

d_tar = [12 33 68 75 110];  %%% distances of targets from ego vehicle
v_tar = [-50 -20 -35 -65 -70]; %%% velocities of targets wrt ego vehicle
theta_tar = [-39 -25 -4 10 24]; %%% AoA of targets wrt ego vehicle
rcs_tar = ones(1,length(d_tar)); %%% rcs of targets
K_ego = 20e12; %%% FMCW  sweep slope of ego vehicle
f_start = 8e09; %%% Sweep start frequency of ego vehicle

%The three conditions below (uncomment any one at a time) are
%to specify the Inter Frame Gap (IFG)  of the Ego Vehicle. 
% Case 1: No Time dithering (IFG ego same as IFG Interferer) 
% Case 2: Fixed Time dithering (IFG different from Interferer)
% Case 3: Randomized Time dithering (IFG different and varying)

%Inter_Frame_Gap = 20e-06.*ones(1,N_f-1);
%Inter_Frame_Gap = 1200e-06.*ones(1,N_f-1);
Inter_Frame_Gap = [800e-06 1200e-06 500e-06 600e-06]; %%%% Random IFG Case

Tx_Start_Time = 0;
Noise_amp = 0.000000005;

%%%%%%%%%%%%%%%%%%% Sim Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This section sets up the simulation environment based on parameters
%%% that you have entered %%%%%

N_int = length(d_int);
att_factor_Int = 2; 
% exponent of distance for free space path loss for interferers
% set to 2 because its a direct ppath

N_tar = length(d_tar);
att_factor_tar = 4; 
% exponent of distance for free space path loss for targets
% set to 4 because it travels to and fro

dist = [d_tar d_int]; %%Just add elements for multiple targets
v = [v_tar v_int];
theta_initial = [theta_tar theta_int];%% between +90 and -90
rcs = [rcs_tar rcs_int];
N_ref = length(dist);
tar_present = 1;

Noise_present = 1;

%%%%%%%%% Defining System Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
global c;
c = 3e08; %% Speed of light
dsamp_rate = 40e06; %% ADC Sampling Rate

Dec_factor = 512; %% Digital to Analog Rate multiplication factor. As you will see in documentation,
% This code tries to emulate analog domain by generating samples with a
% very high sampling rate. The digital domain sampling rate is chosen to be
% similar to real world RADARs. This Dec_factor decides how high the "analog
% rate" will be. 

asamp_rate = dsamp_rate*Dec_factor; %% Equivalent Analog Sampling Rate
t_samp = 1/asamp_rate; %%Sample time in the equivalent Analog Domain

%%%%%%%%% Defining Ego RADAR system Variables (Rx) %%%%%%%%%%%%%%%%%%

Nrx = 16; %% Number of Antennae
N_sample = 512; %% Number of samples per chirp
adc_len = N_sample/dsamp_rate; %% Length of ADC capture
Angle_FFT_Len = 1024; 

%%%%%%%%% Defining Ego RADAR Waveform Variables %%%%%%%%%%%%%%%%%%%%%%%%
% detailed figure is present in the documentation%%%%

Ts = 28e-06; %% Ramp time. Total time = t_start+Ts 
BW = Ts*K_ego; %% Bandwidth
lambda = c/f_start; %% Wavelength
Srx = lambda/2; %% Spacing between antenna elements
adc_start = 10e-06; %% must be greater than t_start
k = K_ego; 
t_start = 8e-06; %% Silent time between chirps
Chirps_Per_Frame=50;
Chirp_Time = (Ts+t_start);
Frame_Time = Chirp_Time*Chirps_Per_Frame;

Sim_Times = Create_Time_Matrix(N_f,Chirps_Per_Frame,Tx_Start_Time,Frame_Time,Inter_Frame_Gap,Chirp_Time);
%This function will define the simulation times

EoSim = Sim_Times(end,end)+Chirp_Time;

%%%%%%%%% Defining Interferer Radar Waveform Variables %%%%%%%%%%%%%%%%%%%%%

%Ts_int = 28e-06*ones(1,N_int); %%% Ramp Time
%BW_int = Ts_int.*K_inter; %% Bandwidth
BW_int = 300e06*ones(1,N_int); %% Badnwidth of interferer fixed to 300 MHz
Ts_int = (BW_int./K_inter); %%% Ramp Time Calculated
k_int = K_inter;
t_start_int = 8e-06*ones(1,N_int); %% silence time between chirps
Chirps_Per_Frame_int = 50*ones(1,N_int); %% Chirps per frame for interferer
for i=1:N_int
    Inter_Frame_Gap_Int(i,:) = Inter_Frame_Gap_Int_val(i).*ones(1,N_f); 
end
Chirp_Time_int = (Ts_int+t_start_int);
Frame_Time_int = Chirp_Time_int.*Chirps_Per_Frame_int;

%%%%%%%%%%%%%%%%%%%%%%% Scenario Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = zeros(N_f,Chirps_Per_Frame,N_ref);
theta = zeros(N_f,Chirps_Per_Frame,N_ref);
vf = zeros(N_f,Chirps_Per_Frame,N_ref);
for i=1:size(Sim_Times,1)
    for j=1:size(Sim_Times,2)
        ToS = Sim_Times(i,j);
        [d(i,j,:),theta(i,j,:),vf(i,j,:)] = Update_Scene_mod(dist,theta_initial,v,ToS);
    end
end

%%%%%%%%%%%%%%% Interferer Transmit Waveform Creation %%%%%%%%%%%%%%%%%%%%%
% This section generates the FMCW waveform for every interferer along with
% a matrix of its chirp times called Int_Time_Matrix

inter = struct();
for i=1:N_int
    t = (0:t_samp:Ts_int(i)-t_samp);
    x_int = exp(1j*2*pi*(f_start_int(i) + 0.5*k_int(i).*t).*t); %% FMCW waveform expression
    x_int = [zeros(1,floor(t_start_int(i)/t_samp)) x_int];
    inter(i).x_int = x_int;
    inter(i).Int_Time_Matrix(:,:) = Create_Time_Matrix_Int(Chirps_Per_Frame_int(i),Tx_Start_Time_Int(i),Frame_Time_int(i),Inter_Frame_Gap_Int(i,:),EoSim,Chirp_Time_int(i));
end

%%%%%%%%%%%%%%%%%% Ego Vehicle and Receive Waveform Creation %%%%%%%%%%%%%%%%%

t = (0:t_samp:Ts-t_samp);
x = cos(2*pi*(f_start + 0.5*k.*t).*t); %% Real FMCW Sweep Waveform generated for one chirp 
x = [zeros(1,floor(t_start/t_samp)) x]; %% silence time between chirps incorporated
len = length(x);

xr = exp(1j*2*pi*(f_start + 0.5*k.*t).*t); %% Complex Receiver waveform (Quadrature LOs equivalent to multiplication by exponent)
xr = [zeros(1,floor(t_start/t_samp)) xr];

x_adc = xr(floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp)); 
%%% Typically ADCs start sampling a little after the chirp start time. This incorporates that

%%%%%%%%%%%% Antenna Array Phase Offsets Calculation %%%%%%%%%%%%%%%%%%%%%%
% This is very similar to theoretical expression used for modeeling channel
% at antenna elements

A = zeros(N_f,Chirps_Per_Frame,Nrx,N_ref+N_int);
mu_tar = ((-2*pi*Srx*(f_start))/c);

%%% Since interferers can have different frequencies, they can have
%%% different phae delays. This part models those
for i = 1:N_int
    mu_int(i) = ((-2*pi*Srx*(f_start_int(i)/2))/c);
end

mu  = zeros(N_f,Chirps_Per_Frame,N_ref+N_int);

mu(:,:,1:N_ref) = mu_tar*sind(theta); %%% For the reflections (target vehicles + interferering vehicles), the mu is calculated


for j = 1:N_int
    mu(:,:,N_ref+j) = mu_int(j)*sind(theta(:,:,N_tar+j)); %% For the interfering waveforms, mu is calculated.
end

for i = 1:Nrx
    A(:,:,i,:) = exp(1j*mu*(i-1)); %% THe antenna offset array is calculated by raising mu to exponent.
end

%%%%%%%%%% Doppler Offset Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Very similar to theoreticla channel modelling of doppler offsets

dopp_offset = zeros(N_f,Chirps_Per_Frame,N_ref+N_int);
f_dopp_tar = f_start*v/c;
f_dopp_int = f_start_int.*v_int/c;

for i = 1:N_f
    for j = 1:Chirps_Per_Frame
        for m = 1:N_ref
            dopp_offset(i,j,m) = exp(1j*2*pi*f_dopp_tar(m)*(j-1)*Chirp_Time);
        end
        for m = 1:N_int
            dopp_offset(i,j,m+N_ref) = exp(1j*2*pi*f_dopp_int(m)*(j-1)*Chirp_Time);
        end
    end
end
toc

%%%%%%%%%%%%%%%%%%% Main Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adcn_mat_ds_frame = zeros(Chirps_Per_Frame,Nrx,N_sample);
adcn_mat_ds_totalsim = zeros(N_f,Chirps_Per_Frame,Nrx,N_sample);
time_exec = zeros(1,N_f);
Interference = zeros(N_int,len);

for j = 1:N_f
    tic
    for m = 1:Chirps_Per_Frame

        %%%%%%%%%%%%%%%%% Receive Waveform Creation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y_mat = zeros(N_ref+N_int,len);
        
        % Generate and Attenuate Target Reflections
        for i=1:N_ref
            Reflected_Chirp = Offset_Chirp(x,(2*d(j,m,i)/c),(j-1)*Chirps_Per_Frame+m,t_samp); %% Generate Target Reflections
            xn = vf(j,m,i).*(Reflected_Chirp)/(d(j,m,i)^att_factor_tar); %% Attenuate reflections
            Y_mat(i,:) = rcs(i)*tar_present*xn;%%%  control to turn targets on or off     
        end
        
        % Generate and Attenuate Interferers
        for i = 1:N_int
            Interference_Chirp = Int_Offset_Chirp(inter(i).x_int,Sim_Times(j,m),Chirp_Time_int(i),Chirp_Time,t_samp,d(j,m,i+N_tar),c,inter(i).Int_Time_Matrix);
            Interference(i,:) = vf(j,m,i+N_tar)*Interference_Chirp/(d(j,m,i+N_tar)^att_factor_Int);
        end
        

        Y_mat(N_ref+1:end,:) = int_present*Interference(:,:); %%%  cotrol for interferer
        Noise_wvfm = Noise_amp*(2*rand(N_ref+N_int,len)-1+2j*rand(N_ref+N_int,len)-1j);
        
        %%% Add dopple
        Y_mat = Y_mat + Noise_present*Noise_wvfm;%%%%%%% control for noise       
        Y_mat_dopp = squeeze(dopp_offset(j,m,:)).*Y_mat;
        
        %%% Add AoA
        Y_Rx = squeeze(A(j,m,:,:))*Y_mat_dopp;
      
        %%%%%%%%%%%%% Downconversion (Mixer) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Y_Rx = Y_Rx(:,floor(adc_start/t_samp):floor((adc_len+adc_start)/t_samp));
        adcn_mat = conj(Y_Rx).*(x_adc);
        
        %%%%%%%%%%%%%%%%%%% ADC Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        adcn_mat_ds = zeros(Nrx,N_sample+1);
        for i=1:Nrx
        adcn_mat_ds(i,:) = decimate(adcn_mat(i,:),Dec_factor); 
        end
    
    adcn_mat_ds_frame(m,:,:)= adcn_mat_ds(:,1:N_sample); %% Chirp Waveform
    end

adcn_mat_ds_totalsim(j,:,:,:) = adcn_mat_ds_frame; %% Frame Waveform
time_exec(j) = toc
end

Avg_Sim_Time_Per_Frame = sum(time_exec)/N_f
save("temp_file1.mat");
save("temp_file2.mat","adcn_mat_ds_totalsim","N_sample");
save("temp_file3.mat","adcn_mat_ds_totalsim","dsamp_rate","BW","Ts","c","N_sample")

 