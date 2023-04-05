clear all;
close all;
load("temp_file1.mat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NDMC_threshold = -110; %% Use this flag to trigger the SALVAGE Algorithm
Ignore_NDMC_flag = 0; %% Use this flag to Bypass the SALVAGE Algorithm. 1  = Bypass, 0 = Consider
Diplay_results_in_polar = 1; %% Turn this off (set to 0) if you want to see results in x-y form (x : distance, y: Angle)

range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
Max_Range = (dsamp_rate*Ts*c)/(2*BW);
range_xaxis = (1:N_sample)*range_res; %%% Range axis for one sided FFT
RA_xaxis = (-N_sample/2:(N_sample/2-1))*range_res; %% Range axis for 2 sided FFT

u = -0.5+1/(2*N_sample)/2:1/N_sample:0.5-1/(2*N_sample); 
angle_vals = asind((lambda/Srx)*u).'; %%% Angle axis for Angle FFT

wi = (kaiser(N_sample,19))'; %% Window for FFTs
wi = wi/sum(wi); %% Normalizing Window

%%%%%%%%%%%$$$%%% Display Overview of chirps andd frames %%%%%%%%%%%%%%%%%%
% THis helps in getting a bird's eye view of how the chirps are looking

Num_Chirp_Disp = 5; %% Choose the number of chirps per frame you want to see for an overview of performance. 

tic
for frm = 1:N_f
    for chp = 1:Num_Chirp_Disp
        Return_Chirp_Number = chp*(Chirps_Per_Frame/Num_Chirp_Disp);
        adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(frm, Return_Chirp_Number,:,:));
        RangeFFT_mat = fftshift(fft(adcn_mat_ds.*wi,N_sample,2));
        
        %%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
        subplot(5,5,(frm-1)*5+chp); 
        plot(RA_xaxis,db(RangeFFT_mat(1,:)));
        grid on;
        title(sprintf('Frame %d Chirp %d',frm,Return_Chirp_Number));
    end
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%% Power Profile Overview %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section plots out the power in negative Bins. comparision with NDMC threshold allows user to see if the SALVAGE algorithm will be triggered 

for cnt = 1:N_f
    for chp = 1:Chirps_Per_Frame
        Chirp_Data = squeeze(adcn_mat_ds_totalsim(cnt, chp,:,:));
        [power_disp(cnt,chp),~] = Check_Int(NDMC_threshold,Chirp_Data);
    end
    figure(2)
    hold on
    plot(db(power_disp(cnt,:)),DisplayName="Frame "+cnt);
end
figure(2)
hold off
ylabel("dB")
xlabel("Chirp Number")
title("Magnitude of power in  negative frequency bins")
yline(NDMC_threshold,Color="red",LineStyle="--",DisplayName="NDMC Threshold");
ylim([-250 0]);
legend;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% System %%%%%%%%%%%%%%%%%%%

for cnt = 1:N_f
    Frame_Data = squeeze(adcn_mat_ds_totalsim(cnt,:,:,:));
    Frame_Data_singant = squeeze(adcn_mat_ds_totalsim(cnt,:,1,:));
    for chp = 1:Chirps_Per_Frame
        Chirp_Data = squeeze(adcn_mat_ds_totalsim(cnt, chp,:,:));
        [~,NDMC_flag] = Check_Int(NDMC_threshold,Chirp_Data); %%% NDMC Flag turns to 1 if the negative bins power indicates interference
        if(NDMC_flag ==1 && Ignore_NDMC_flag==0)
            break;
        else
            DataSet = Chirp_Data;
            RangeFFT_CUT = (fft(DataSet.*wi,N_sample,2));
            [Locations,AoA] = Find_Objects(N_sample,range_res,RangeFFT_CUT,angle_vals,f_start,Srx,Nrx,c); %% Finds object peaks from RangeFFT Data
            [Performance_chp(chp,:),Targets_mat,Eval_mat] = Evaluate_Results(Locations,AoA,squeeze(d(cnt,chp,:))',squeeze(theta(cnt,chp,:))',Nrx,range_res); %% Compares object peaks to actul targets
            if(chp == 1)
                Targets_obs = Eval_mat;
                Targets_act = Targets_mat;
            else
                Targets_obs = [Targets_obs Eval_mat];
                Targets_act = [Targets_act Targets_mat]; %% Accumulate target data from every chirp to overlay later
            end
        end
    end
    if(NDMC_flag ==1 && Ignore_NDMC_flag==0)
        %%% This section is the SALVAGE Algorithm %%%%
       [ind_min,CUT,chirp_final] = NDMC(Frame_Data_singant,wi,N_sample); %% Construct best case RangeFFT chirp based on NDMC. Also return the indices of min samples to plot histogram
       DataSet = squeeze(Frame_Data(CUT,:,:));
       RangeFFT_CUT = (fft(DataSet.*wi,N_sample,2));

        figure(3)
        subplot(N_f,3,3*(cnt-1)+1)
        plot(RA_xaxis,db(abs(fftshift(chirp_final))));
        xlabel("Distance (m)")
        ylabel("dB")
        title("Reconstructed Distance FFT by NDMC")
        subplot(N_f,3,3*(cnt-1)+2)
        histogram(ind_min);
        xlabel("Chirp Number")
        title("Histogram of Samples Contributed to Best Chirp")
        ylim([0 512])
        subplot(N_f,3,3*(cnt-1)+3)
        plot(RA_xaxis,db(abs(fftshift(RangeFFT_CUT(1,:)))));
        xlabel("Distance (m)")
        ylabel("dB")
        title("Frame Number "+cnt+", Chirp Number "+CUT)

       [Locations,AoA] = Find_Objects(N_sample,range_res,RangeFFT_CUT,angle_vals,f_start,Srx,Nrx,c);
       [Performance_frm(cnt,:),Targets_act,Eval_mat] = Evaluate_Results(Locations,AoA,squeeze(d(cnt,CUT,:))',squeeze(theta(cnt,CUT,:))',Nrx,range_res);
       Targets_obs = Eval_mat;
    else
       Performance_frm(cnt,:) = sum(Performance_chp,1)/Chirps_Per_Frame; %%% Find Average values in frame for recall and precision
    end
    
    if(Diplay_results_in_polar==1)
        figure(40+cnt)
        s = polarscatter(deg2rad(Targets_act(2,:)),Targets_act(1,:),18,'filled',MarkerFaceColor="red",MarkerFaceAlpha="flat");
        s.AlphaData = 0.005*ones(1,length(Targets_act(1,:)));
        hold on
        s2 = polarscatter(deg2rad(Targets_obs(2,:)),Targets_obs(1,:),9,'filled',MarkerEdgeColor="none",MarkerFaceAlpha="flat",MarkerFaceColor="blue");
        s2.AlphaData = 0.005*ones(1,length(Targets_obs(1,:)));
        title("Frame" + cnt)
        rlim([0 150])
        thetalim(([-90 90]))
        legend("Targets","Detections");
    else
      figure(4)
      subplot(1,N_f,cnt)
      s = scatter(Targets_act(1,:),Targets_act(2,:),18,'filled',MarkerFaceColor="red",MarkerFaceAlpha="flat");
      s.AlphaData = 0.005*ones(1,length(Targets_act(1,:)));
      hold on
      s2 = scatter(Targets_obs(1,:),Targets_obs(2,:),9,'filled',MarkerEdgeColor="none",MarkerFaceAlpha="flat",MarkerFaceColor="blue");
      s2.AlphaData = 0.005*ones(1,length(Targets_obs(1,:)));
      legend("Targets","Detections");
      xlabel("Distance (metres)")
      ylabel("AoA (Degrees)")
      title("Frame" + cnt)
      xlim([0 150])
      ylim([-90 90])
    end
end

figure(5)
plot(1:1:N_f,Performance_frm(:,1),Marker="*");
hold on
plot(1:1:N_f,Performance_frm(:,2),Marker="*");
hold off
title("Performance Metrics")
legend("Recall","Precision");
xlabel("Frame")
ylabel("Percentage(%)")
ylim([-5 105])

avg_recall = sum(Performance_frm(:,1),"all")/N_f
avg_precision = sum(Performance_frm(:,2),"all")/N_f