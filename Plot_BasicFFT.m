clear all;
close all;
load("temp_file3.mat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
RA_xaxis = (-N_sample/2:(N_sample/2-1))*range_res;
wi = (kaiser(N_sample,19))';
wi = wi/sum(wi);        
frm = 1;
chp = 20;
adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(frm, chp,:,:));
RangeFFT_mat = fftshift(fft(adcn_mat_ds.*wi,N_sample,2));

figure(1)
subplot(1,2,1)
plot(RA_xaxis,db(RangeFFT_mat(1,:)));
grid on;
title(sprintf('RangeFFT: Frame %d Chirp %d',frm,chp));
ylabel("dB")
xlabel("bins")
subplot(1,2,2)
plot(real(adcn_mat_ds(1,:)));
grid on;
title(sprintf('Time Domain: Frame %d Chirp %d',frm,chp));
xlabel("samples")
 
figure(1)
subplot(1,2,1)
plot(RA_xaxis,db(RangeFFT_mat(1,:)));
grid on;
title(sprintf('RangeFFT: Frame %d Chirp %d',frm,chp));
ylabel("dB")
xlabel("bins")
subplot(1,2,2)
plot(real(adcn_mat_ds(1,:)));
grid on;
title(sprintf('Time Domain: Frame %d Chirp %d',frm,chp));
xlabel("samples")

%%%%%%%%%%%% Excision %%%%%%%%%%%%%%%%%%%%%%%%%%
start_index = 96;
stop_index = 490;
adcn_mat_ds2 = [adcn_mat_ds(1,1:start_index) adcn_mat_ds(1,stop_index:end)];
wi2 = (kaiser(length(adcn_mat_ds2),19))';
wi2 = wi2/sum(wi2); 
RangeFFT_mat_cut = fftshift(fft(squeeze(adcn_mat_ds2).*wi2,N_sample,2));

figure(2)
plot(RA_xaxis,db(RangeFFT_mat_cut));
hold on
plot(RA_xaxis,db(RangeFFT_mat(1,:)));
grid on;
title(sprintf('RangeFFT: Frame %d Chirp %d',frm,chp));
legend("With Interference+Excision","With interference")
ylabel("dB")
xlabel("bins")

