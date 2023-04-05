%%
clear all;
close all;
load("temp_file2.mat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LMS %%%%%%%%%%%%%%%%%%%
l = 10; %%%filter length
mu = 0.49; %%% Control Variable
frm = 1;
chp = 30;
wi = (kaiser(N_sample,19))';
wi = wi/sum(wi);    
Return_Frame_Number = frm;
Return_Chirp_Number = chp;

adcn_mat_ds_filt = squeeze(adcn_mat_ds_totalsim(Return_Frame_Number, Return_Chirp_Number,:,:));
RangeFFT_mat_filt = fftshift(fft(adcn_mat_ds_filt.*wi,N_sample,2));
RangeFFTn_filt = RangeFFT_mat_filt(1,:);  %% Since this is only Range FFT, we pick out the first antenna waveform directly
x = RangeFFTn_filt(0.5*length(RangeFFTn_filt)+1:end)/max(abs(RangeFFTn_filt(0.5*length(RangeFFTn_filt)+1:end)));
d = conj(fliplr(RangeFFTn_filt(1:0.5*length(RangeFFTn_filt))))/max(abs(RangeFFTn_filt(0.5*length(RangeFFTn_filt)+1:end)));
P = sum(abs(d).^2);
step = 2/(1500*P); %%% Step size formula picked up from reference paper

tot = RangeFFT_mat_filt(1,:);
wo = [1 zeros(1,l-1)]';
fi = zeros(1,l-1);
err = zeros(1,N_sample/2)';

for j = 1:N_sample/2
    fi = [d(j) fi(1:l-1)];
    fo = wo'*transpose(fi);
    err(j) = x(j)-fo;
    wo = wo + step.*fi'.*conj(err(j));
end
y_out_p = err';

lms = dsp.LMSFilter('Length',l,'StepSize',mu);
[y,e] = lms(x',d');
y_out = x'-y;

figure(12)
subplot(2,3,1)
plot(db(abs(x)));
title("RangeFFT Positive Bins")
ylabel("dB")
xlabel("bins")
subplot(2,3,2)
plot(db(abs(d)));
title("RangeFFT Negative Bins")
ylabel("dB")
xlabel("bins")
subplot(2,3,3)
plot(db(abs(tot)));
title("RangeFFT")
ylabel("dB")
xlabel("bins")

subplot(2,3,4)
plot(db(abs(y_out)));
title("MATLAB LMS Filter: length = "+l+", mu = "+mu)
ylabel("dB")
xlabel("bins")
subplot(2,3,5)
plot(db(abs(y_out_p)));
title("LMS Filter from [13] length = "+l+", mu = "+step)
ylabel("dB")
xlabel("bins")
figure(2)
hold on
plot(db(y_out_p));
plot(db(y_out));
plot(db(x));
hold off
legend("My LMS","MATLAB LMS","Unfiltered")
title("Comparision")
ylabel("dB")
xlabel("bins")