%   radar_data = load("04_29_with_interference2_adcdata\adcData2.mat");
%   radar_data = radar_data.adcData;
%   radar_data = reshape(radar_data,512,size(radar_data,2),192);
%   radar_data = radar_data(:,:,antenna_azimuthonly);
%   real_data = radar_data;
%   radar_chirp = squeeze(radar_data(:,chirp_idx,:)).';
    load("C:\Shubham\RADAR TestBed\From Server\Data_Sets_Sim\Output_Para_17.mat")
    load("C:\Shubham\RADAR TestBed\From Server\Data_Sets_Sim\Input_Para17.mat")
    dsamp_rate = 20e06;
    BW = 300e06;
    Ts = 28e-06;
    c = 3e08;
    N_sample = 512;
    Angle_FFT_Len = 512;
    lambda = c/77e09;
    Srx = 16;
    range_res = (dsamp_rate*Ts*c)/(2*BW*N_sample);
    Max_Range = (dsamp_rate*Ts*c)/(2*BW);
    range_xaxis = (1:N_sample)*range_res;
    u = -0.5:1/Angle_FFT_Len:0.5-1/Angle_FFT_Len; %% x-axis in world of sin(theta)*(sep/lambda)
    angle_vals = asind((lambda/Srx)*u).';

    adcn_mat_ds = squeeze(adcn_mat_ds_totalsim(1,:,1,:));
    rangefft = fft(adcn_mat_ds,512,2);
    doppler_fft = fftshift(fft(rangefft,50,1),1);
    %doppler_fft = [doppler_fft; doppler_fft(1,:)];

    n_angle_fft_size = 512;
    indices_1D = 0:n_angle_fft_size-1;
    
    %range_resolution = 3e8/2/(0.512e9);
    range_resolution = range_res;
    range_vals = indices_1D*range_resolution;
    
%     sine_theta = -2*((-n_angle_fft_size/2:n_angle_fft_size/2)/n_angle_fft_size);
%     cos_theta = sqrt(1-sine_theta.^2);
%     theta = -asin(sine_theta);
%     
%     [R_mat, sine_theta_mat] = meshgrid(range_vals,sine_theta);
%     [~, cos_theta_mat] = meshgrid(indices_1D,cos_theta);
%     x_axis = R_mat.*cos_theta_mat;
%     y_axis = R_mat.*sine_theta_mat;
%     figure(1)
%     subplot(1,3,1)
%     surf(-y_axis, x_axis, db(abs(anglefft)), 'EdgeColor', 'None');
%     
%     ix = find(imregionalmax(db(abs(anglefft))));
%     [a,b] = ind2sub(size(db(abs(anglefft))),ix);
    
    
    subplot(2,3,1)
    plot(range_vals,db(abs(rangefft(1,:))));
    subplot(2,3,2)
    imagesc((db(abs(doppler_fft))));    %axis([0 300 40 120]);
    subplot(2,3,3)
    plot(db(abs(rangefft))','b');
    subplot(2,3,4)
    doppler_fft_2 = fftshift(fft(rangefft,50,1),1);
    imagesc(db(abs(doppler_fft_2))');