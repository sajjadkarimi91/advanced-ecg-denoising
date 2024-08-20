% electro-mechanicl modeling of heart rates

% A sample MATLAB script for reading the ECG-PCG data files and load the
% R-peaks to extract the S1/S2 indexes of the PCG channels
%
% Dependencies: This script uses some functions from the Open-Source
% Electrophysiological Toolbox (OSET): https://github.com/alphanumericslab/OSET.git
% PhysioNet-Cardiovascular-Signal-Toolbox: https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
% Ecg-Kit :https://github.com/marianux/ecg-kit

% NOTE: This code has not been optimized for any specific dataset and is
%
% CITE:
% 1- The EPHNOGRAM on PhysioNet
% 2- ?????
% 3- R. Sameni, The Open-Source Electrophysiological Toolbox (OSET), v 3.14, URL: https://github.com/alphanumericslab/OSET
%
% By: Sajjad Karimi and Reza Sameni
% Email: sajjadkarimi91@gmail.com
% Email: reza.sameni@gmail.com
% JAN 2024
%

clear;
close all;
clc;

root_path = './under-development-codes\03 ecg-denoising-rpeaks';
addpath(genpath(root_path))

oset_path = './OSET';% enter the path of OSET toolbox
addpath(genpath(oset_path))

lsim_path = 'D:\PHD codes\chmm-lsim-karimi-toolbox';  % download from https://github.com/sajjadkarimi91/chmm-lsim-matlab-toolbox
addpath(genpath(lsim_path))



db_folder = './under-development-codes/ephnogram/MAT'; % .mat files folder
rpeak_folder = './under-development-codes/01 paper ephnogram database/RPEAKS'; % .mat files folder
db_csvfile = './under-development-codes/ephnogram/ECGPCGSpreadsheet.csv'; % information for recordings
path_save = './under-development-codes/ephnogram/MATDN'; % .mat files folder

mkdir(path_save)

local_db_files = dir([db_folder '/*.mat']); % list of all mat files

% pre-processing parameters
fs = 1000.0; % post-decimation sampling frequency


for m = 1:length(local_db_files)

    tic
    disp(m)
    % load data and rpeaks
    in_fname = local_db_files(m).name;
    in_fname_full = [db_folder '/' in_fname];
    dat = load(in_fname_full); % load the data file

    % decimate the ECG and PCG channel to fs
    ecg_raw_org = decimate(dat.ECG, round(dat.fs/fs));

    % ecg_raw_org = ecg_raw_org+0.02*sin(2*pi*50/fs*(1:length(ecg_raw_org))); % powerline noise

    % NOTCH FILTERING THE ECG
    fc = 50.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = filtfilt(b, a, ecg_raw_org); % zero-phase non-causal filtering

    Qfactor = 25; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = filtfilt(b, a, ecg_raw); % zero-phase non-causal filtering

    % NOTCH FILTERING THE ECG
    fc = 100.0; % powerline frequency
    Qfactor = 5; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = filtfilt(b, a, ecg_raw); % zero-phase non-causal filtering

    ecg_raw = ecg_raw - movmean(movmedian(ecg_raw,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);

    % INTERPOLATE THE HEART RATE SEQUENCES
    t_second = (0 : length(ecg_raw)-1)/fs;
    tmin = t_second/60;% time in minute unit

    ecg_rpeaks_table = readtable([rpeak_folder,'/',in_fname(1:end-4),'_rpeaks.csv'],"FileType","spreadsheet");% read Rpeak index
    ecg_rpeaks_index = ecg_rpeaks_table.R_peak_indexes_ms;

    if contains(ecg_rpeaks_table.ECGNotes{1},'Disconnected')
        save_fname_full = [path_save '/' in_fname];
        ecg_denoised = ecg_raw;
        save( save_fname_full, 'ecg_denoised','fs');
        continue;
    end

    rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index);
    avg_intervals_ecg = movmean(rr_intervals_ecg,[30,30]);
    avg_intervals_ecg = [avg_intervals_ecg;avg_intervals_ecg(end)];
    rr_intervals_ecg = [rr_intervals_ecg;rr_intervals_ecg(end)];


    ecg_denoised = 0*ecg_raw;
    sample_70ms = round(fs*0.07);
    sample_200ms = round(fs*0.2);
    sample_300ms = round(fs*0.3);
    sample_350ms = round(fs*0.35);
    sample_rem = ceil(prctile(rr_intervals_ecg,99)) - sample_350ms;

    pqrs_bloks = zeros(length(ecg_rpeaks_index)-2,sample_200ms+sample_70ms+1);
    t_bloks = zeros(length(ecg_rpeaks_index)-2,round(sample_350ms*min(2,max(max(1,avg_intervals_ecg/(2*sample_350ms)))))-sample_70ms);
    remained_bloks = zeros(length(ecg_rpeaks_index)-2,sample_rem);

    for p = 2:length(ecg_rpeaks_index)-1
        this_qrs_index = ecg_rpeaks_index(p)-min(sample_200ms,floor(0.3*rr_intervals_ecg(p-1))):ecg_rpeaks_index(p)+sample_70ms;
        if any(this_qrs_index<1)||any(this_qrs_index>length(ecg_raw))
            continue;
        end
        pqrs_bloks(p-1,end-length(this_qrs_index)+1:end) = ecg_raw(this_qrs_index) ;%- linspace(mean(ecg_raw(this_qrs_index(1:fs*0.01))), mean(ecg_raw(this_qrs_index(end-round(fs*0.01):end))), 2*sample_70ms+1);

        %this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+min(sample_300ms,floor(0.7*rr_intervals_ecg(p)));
        sample_T = round( min(sample_350ms*min(2,max(1,avg_intervals_ecg(p)/(2*sample_350ms))),floor(0.7*rr_intervals_ecg(p))));
        this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+sample_T;
        if any(this_t_index>length(ecg_raw))
            continue;
        end

        t_bloks(p-1,1:length(this_t_index)) =  ecg_raw(this_t_index);% - linspace(init_mn, mean(ecg_raw(this_t_index(end-round(fs*0.015):end))), length(this_t_index));

        this_rem_index = ecg_rpeaks_index(p)+sample_T+1:ecg_rpeaks_index(p)+min([sample_rem+sample_350ms,rr_intervals_ecg(p)-min(sample_200ms,round(0.3*rr_intervals_ecg(p)))-1]);
        if any(this_rem_index>length(ecg_raw))
            continue;
        end
        remained_bloks (p-1,1:length(this_rem_index)) =  ecg_raw(this_rem_index) ;

    end

    ecg_blocks = [pqrs_bloks,t_bloks,remained_bloks];

    qrs_bloks_mv = movmean(pqrs_bloks,[10,10],1,'omitmissing');
    t_bloks_mv = movmean(t_bloks,[10,10],1,'omitmissing');


    max_val = max(prctile(pqrs_bloks(:),99.5),prctile(t_bloks(:),99.5));
    min_val = min(prctile(pqrs_bloks(:),0.5),prctile(t_bloks(:),0.5));

    ecg_blocks(ecg_blocks>max_val) = max_val;
    ecg_blocks(ecg_blocks<min_val) = min_val;
    ecg_blocks = (ecg_blocks - min_val) / (max_val-min_val);
    ecg_blocks_mv = movmean(ecg_blocks(:,1:sample_300ms+sample_70ms),[10,10],1,'omitmissing');


    ecg_noise = ecg_blocks(:,1:sample_300ms+sample_70ms) - ecg_blocks_mv;
    noise_power = movmean( sum(ecg_noise.^2,2),[30,30]);
    signal_power = movmean( sum(ecg_blocks_mv(:,1:sample_300ms+sample_70ms).^2,2),[30,30]);

    snr_signal_db = 10*log10(signal_power./noise_power);

    combining_weights = [exp((snr_signal_db-20)/10),exp((snr_signal_db-20)/10),exp((snr_signal_db-20)/10),exp((20-snr_signal_db)/10)];
    combining_weights = combining_weights./sum(combining_weights,2);

    ecg_noise = [ecg_noise,circshift(ecg_noise(:,end:-1:1),1,1),circshift(ecg_noise,1,5),circshift(ecg_noise(:,end:-1:1),1,3),circshift(ecg_noise,1,3)];
    ecg_noise = [ecg_noise,ecg_noise];
    ecg_noise = ecg_noise(:,1:size(ecg_blocks,2));

    PSD = abs(fft2(ecg_noise, size(ecg_noise,1), size(ecg_noise,2))).^2;
    PSD = movmean(PSD,[3,3]);
    PSD = movmean(PSD,[3,3],2); % noise power spectrum

    y_est_bm3d = BM3D(ecg_blocks, PSD ,'refilter');
    y_est_bm3d = y_est_bm3d*(max_val-min_val) + min_val;

    y_est_curvelet = ACT_filter(ecg_blocks, PSD, 'ksigma');
    y_est_curvelet = y_est_curvelet*(max_val-min_val) + min_val;

    y_est_adap_curvelet = ACT_filter(ecg_blocks, PSD);
    y_est_adap_curvelet = y_est_adap_curvelet*(max_val-min_val) + min_val;

    % RF3D parameters
    profile           = 'np';  % Complexity profile:
    transform_2D_name = 'dct'; % 2D block transform used in RF3D (it also defines the PSDs)
    do_wiener         = false;  % Wiener Filtering (true or false).
    erf3d             = true;  % Enhanced FPN suppression (true or false).
    sigma_rnd = 1;            % Scaling factor of random noise (RND)
    sigma_fpn = 1;            % Scaling factor of fixed-pattern noise (FPN)

    batch_size = 1000;
    num_batch = floor(size(pqrs_bloks_on,2)/batch_size);
    y_est_rf3d = [];

    for b = 1:num_batch
        if b~=num_batch
            ind_filter = (b-1)*batch_size+1:b*batch_size;
        else
            ind_filter = (b-1)*batch_size+1:size(pqrs_bloks_on,2);
        end
        y_est_rf3d1 = RF3D(ecg_blocks(ind_filter,:), sigma_rnd, sigma_fpn, [], [], transform_2D_name, erf3d, do_wiener, profile);
        y_est_rf3d = [y_est_rf3d;y_est_rf3d1];
    end

    % if (size(ecg_blocks,1)>3200 && m==60) ||  m>=61
    %     y_est_rf3d1 = RF3D(ecg_blocks(1:1000,:), sigma_rnd, sigma_fpn, [], [], transform_2D_name, erf3d, do_wiener, profile);
    %     y_est_rf3d2 = RF3D(ecg_blocks(1001:end,:), sigma_rnd, sigma_fpn, [], [], transform_2D_name, erf3d, do_wiener, profile);
    %     y_est_rf3d = [y_est_rf3d1;y_est_rf3d2];
    % else
    %     y_est_rf3d = RF3D(ecg_blocks, sigma_rnd, sigma_fpn, [], [], transform_2D_name, erf3d, do_wiener, profile);
    % end

    y_est_rf3d = y_est_rf3d*(max_val-min_val) + min_val;

    % final combination
    y_est = combining_weights(:,1).*y_est_bm3d + combining_weights(:,2).*y_est_curvelet +  combining_weights(:,3).*y_est_adap_curvelet +  combining_weights(:,4).*y_est_rf3d;

    qrs_bloks_denoise = y_est(:,1:size(pqrs_bloks,2));
    t_bloks_denoise = y_est(:,size(pqrs_bloks,2)+1:size(pqrs_bloks,2)+size(t_bloks,2));
    remained_bloks_denoise = y_est(:,size(pqrs_bloks,2)+size(t_bloks,2)+1:end);

    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_rpeaks_index(p)-min(sample_200ms,floor(0.3*rr_intervals_ecg(p-1))):ecg_rpeaks_index(p)+sample_70ms;
        if any(this_qrs_index<1)||any(this_qrs_index>length(ecg_raw))
            continue;
        end
        ecg_denoised(this_qrs_index) = qrs_bloks_denoise(p-1,end-length(this_qrs_index)+1:end);

        %this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+min(sample_300ms,floor(0.7*rr_intervals_ecg(p)));
        this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:ecg_rpeaks_index(p)+min(sample_350ms*min(2,max(1,avg_intervals_ecg(p)/(2*sample_350ms))),floor(0.7*rr_intervals_ecg(p)));

        if any(this_t_index>length(ecg_raw))
            continue;
        end
        ecg_denoised(this_t_index) = t_bloks_denoise(p-1,1:length(this_t_index));

        this_rem_index = ecg_rpeaks_index(p)+sample_300ms+1:ecg_rpeaks_index(p)+min([sample_rem+sample_300ms,rr_intervals_ecg(p)-min(sample_200ms,round(0.3*rr_intervals_ecg(p)))-1]);
        if any(this_rem_index>length(ecg_raw))
            continue;
        end

        ecg_denoised(this_rem_index) = remained_bloks_denoise (p-1,1:length(this_rem_index));

    end

    close all

    a = figure('Position', [130 130 1500 800]);
    lg = {};
    plot(t_second,ecg_raw) ;lg = cat(2, lg, {'ECG'});
    hold on
    plot(t_second,ecg_denoised,LineWidth=1.5) ;lg = cat(2, lg, {'ECG-env'});
    plot(t_second,ecg_raw_org) ;lg = cat(2, lg, {'ECG'});
    grid on
    legend(lg,'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
    xlabel('time (sec)',Interpreter='latex',FontSize=14)

    toc
    pause(1)
    save_fname_full = [path_save '/' in_fname];

    save( save_fname_full, 'ecg_denoised','fs');

end


%%

for m = 1:length(local_db_files)

    tic
    disp(m)
    % load data and rpeaks
    in_fname = local_db_files(m).name;
    in_fname_full = [db_folder '/' in_fname];
    dat = load(in_fname_full); % load the data file

    save_fname_full = [path_save '/' in_fname];
    load( save_fname_full, 'ecg_denoised');

    % decimate the ECG and PCG channel to fs
    ecg_raw_org = decimate(dat.ECG, round(dat.fs/fs));

    % ecg_raw_org = ecg_raw_org+0.02*sin(2*pi*50/fs*(1:length(ecg_raw_org))); % powerline noise

    % NOTCH FILTERING THE ECG
    fc = 50.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = filtfilt(b, a, ecg_raw_org); % zero-phase non-causal filtering

    % NOTCH FILTERING THE ECG
    fc = 100.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = filtfilt(b, a, ecg_raw); % zero-phase non-causal filtering

    ecg_raw = ecg_raw - movmean(movmedian(ecg_raw,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);

    t_second = (0 : length(ecg_raw)-1)/fs;
    tmin = t_second/60;% time in minute unit

    close all

    a = figure('Position', [130 130 1500 800]);
    lg = {};
    plot(t_second,ecg_raw) ;lg = cat(2, lg, {'ECG'});
    hold on
    plot(t_second,ecg_denoised,LineWidth=1.5) ;lg = cat(2, lg, {'ECG-env'});
    plot(t_second,ecg_raw_org) ;lg = cat(2, lg, {'ECG'});
    grid on
    legend(lg,'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
    xlabel('time (sec)',Interpreter='latex',FontSize=14)

    toc
    pause()


end

