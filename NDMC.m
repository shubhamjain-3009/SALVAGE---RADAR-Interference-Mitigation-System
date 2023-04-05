% THis function constructs the Best case RangeFFT chirp based on NDMC

function [ind_min,CUT,chirp_final] = NDMC (adcn_mat_frame,wi,N_sample)
ind_min = zeros(1,N_sample);
chirp_final = zeros(1,N_sample);
RangeFFT_mat2 = (fft(adcn_mat_frame.*wi,N_sample,2));
for j = 1:N_sample
    temp = abs(RangeFFT_mat2(:,j));
    [~,ind_min(j)] = min(temp);
    chirp_final(j) = RangeFFT_mat2(ind_min(j),j);
end
CUT = mode(ind_min);