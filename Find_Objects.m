% This function find object peaks from RangeFFT Data. This is done in to
% stages and explained in detail in documentation

function [Distances,Angles] = Find_Objects(N_sample,range_res,RangeFFT,angle_vals,f_start,Srx,Nrx,c) 

Range_axis = (1:N_sample)*range_res;
[pks,locns,~,p] = findpeaks(db(RangeFFT(1,:)),Range_axis);
Threshold = 0.15*max(p);
Distances = Peak_Finder_Rev3(pks,locns,p,Threshold); %% This function chooses relevant peaks out of all those returned by matlab in FFT.

if(Distances(1)==0)
    Angles = [0];
else
    [~,Distance_Indices] = intersect(Range_axis,Distances);
    AngleFFT_mat = zeros(length(Distance_Indices),N_sample);
    Angles = zeros(1,length(Distances));
    for i = 1:length(Distances)
        for co = 1:length(angle_vals)
            AngleFFT_mat(i,co) =  exp((+1j*2*pi*f_start*Srx*sind(angle_vals(co))*[0:Nrx-1])/c)*RangeFFT(:,Distance_Indices(i));
        end
        AngleFFT_mat(i,:) = abs(AngleFFT_mat(i,:))/max(abs(AngleFFT_mat(i,:)));
        [~,Temp] = findpeaks(AngleFFT_mat(i,:),angle_vals,"SortStr","descend");
        Angles(i) = Temp(1);
    end
end
    