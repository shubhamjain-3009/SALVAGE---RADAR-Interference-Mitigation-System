% THis function compares detected r and theta values to actual r and theta
% values. Tolerance is defined by range, angle resolution and modelling
% inaccuracy for ToF. Based on results, True positive rate, and Precision
% is calculated and returned

function [Performance_p,Targets_obs,Eval_mat] = Evaluate_Results(Locations,AoA,dist,theta_initial,Nrx,range_res)

Targets_det = [Locations;-AoA];
Targets_act = [dist;theta_initial];
tag = zeros(1,size(Targets_act,2));
tag_det = zeros(1,size(Targets_det,2));

[x_det,y_det] = pol2cart(Targets_det(2,:),Targets_det(1,:));
[x_act,y_act] = pol2cart(Targets_act(2,:),Targets_act(1,:));
theta_Res = asind((2/Nrx));
tol = 1;

dev_max = sqrt(asind((2/Nrx))^2+range_res^2);
for i = 1:length(Targets_act)
    [Diff,Ind] = min(sqrt((x_act(i)-x_det).^2+(y_act(i)-y_det).^2));    
    if(abs(Targets_det(1,Ind)-Targets_act(1,i))<(range_res+tol)&&abs(Targets_det(2,Ind)-Targets_act(2,i))<(theta_Res))
        tag(i) = Ind;
        tag_det(Ind) = i;
    end
end
Targets_obs = [Targets_act;tag];
Eval_mat = [Targets_det;tag_det];

%%%TP
temp = unique(tag);
if(temp(1)==0)
    TP = length(temp)-1;
else
    TP = length(temp);
end
TP_p = (TP/length(tag))*100;

%%%%FP

FP = length(tag_det)-nnz(unique(tag_det));
PP_p = (TP/(TP+FP))*100;

Performance_p(1:2) = [TP_p PP_p];
