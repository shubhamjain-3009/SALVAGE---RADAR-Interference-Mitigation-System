%% This function rorates the original_chirp based ont he time offset provided by user

function [reflected_chirp] = Offset_Chirp(original_chirp,offset,chirp_num,t_samp)
n = length(original_chirp);
offset_n = floor(offset/t_samp);
if(chirp_num==1)
    reflected_chirp = [zeros(1,offset_n) original_chirp(1:n-offset_n)];
else
    reflected_chirp = [original_chirp(n-offset_n-1:end) original_chirp(1:n-offset_n)];
end

reflected_chirp = reflected_chirp(1:n);