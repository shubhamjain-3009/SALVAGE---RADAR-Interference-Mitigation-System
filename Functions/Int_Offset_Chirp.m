
% This function first find the closest match of an interferer chirp start with
% the simulation time under consideration. Based on that, it finds the
% offset by which the interferer chirp needs to be rotated to model it as an
% interfering waveform

function[int_chirp] = Int_Offset_Chirp(original_int_chirp,te,Chirp_Time_Int,Chirp_Time_Ego,t_samp,d,c,Time_Matrix_Int)
ti = Find_Time_Match(te,Time_Matrix_Int);
ti = ti + d/c;
offset = abs(ti-te);
offset_n = floor(offset/t_samp);
n = length(original_int_chirp);
len = Chirp_Time_Ego/t_samp; %%Analog sampling time since we are dealing with RF waveforms

if ((offset>Chirp_Time_Int)&&(offset>Chirp_Time_Ego))
    int_chirp = zeros(1:len);
    disp("here")
else
    core = [original_int_chirp(offset_n:n) original_int_chirp(1:offset_n-1)];
    int_chirp = core;
    while(length(int_chirp)<len)
        int_chirp = [int_chirp core];
    end
    int_chirp = int_chirp(1:len);
end
