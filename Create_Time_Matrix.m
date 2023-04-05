% This function creates a matrix of time values at which the simulation
% needs to run. It is at these time values that the scene (distance, theta) is
% updated and a snapshot of the received wavefor is captured. This
% corresponds to the start of every chirp and that is what this function
% populates in the Sim_Times matrix

function [Sim_Times] = Create_Time_Matrix(N_f,Chirps_Per_Frame,Tx_Start_Time,Frame_Time,Inter_Frame_Gap,Chirp_Time)
Sim_Times = zeros(N_f,Chirps_Per_Frame);
k = length(Inter_Frame_Gap);
for i = 1:N_f
    if(i==1)
        Sim_Times(i,1)=Tx_Start_Time;
    else
        Sim_Times(i,1) = Sim_Times(i-1,1) + Frame_Time + Inter_Frame_Gap(double(mod(i-1,k)==0)*k+mod(i-1,k));
    end
    for j = 2:Chirps_Per_Frame
        Sim_Times(i,j) = Sim_Times(i,j-1) + Chirp_Time;
    end
end
