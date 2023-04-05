%% This function calculates power in negative FFT bins from time data. And compares them with NDMC threshold to see if they have been interfered with

function [power,NDMC_flag] = Check_Int(NDMC_threshold,Chirp_Data)

N_sample = length(Chirp_Data);
NDMC_flag = 0;

ipfft = fftshift(fft(Chirp_Data,N_sample,2));
d = ipfft(1,1:0.5*length(ipfft(1,:)));
power = sum(abs(d).^2);

if(db(power)>NDMC_threshold)
NDMC_flag = 1;
end
