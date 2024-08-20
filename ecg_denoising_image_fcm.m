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
    ecg_raw = ecg_raw_org;
    [Pxx,freq] = pwelch(ecg_raw,[],[],[],fs);
    fc = 50.0; % powerline frequency
    if sum(Pxx(freq>fc-1 & freq<fc+1 ))/sum(Pxx(freq<fc-1&freq>1))> 0.05 || sum(Pxx(freq>2*fc-1 & freq<2*fc+1 ))/sum(Pxx(freq<fc-1&freq>1))> 0.05
        % NOTCH FILTERING THE ECG
        Qfactor = 25; % Q-factor of the notch filter
        Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
        [b,a] = iirnotch(Wo, BW); % design the notch filter
        ecg_raw = filtfilt(b, a, ecg_raw); % zero-phase non-causal filtering

        % NOTCH FILTERING THE ECG
        if 2*fc<fs/2
            Qfactor = 25; % Q-factor of the notch filter
            Wo = 2*fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
            [b,a] = iirnotch(Wo, BW); % design the notch filter
            ecg_raw = filtfilt(b, a, ecg_raw); % zero-phase non-causal filtering
        end

    end

    ecg_raw = ecg_raw - movmean(movmedian(ecg_raw,[floor(0.3*fs),floor(0.3*fs)]),[floor(0.15*fs),floor(0.15*fs)]);
    ecg_raw = lp_filter_zero_phase(ecg_raw, 45/fs);

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
    sample_10ms = round(fs*0.01);
    sample_250ms = round(fs*0.25);
    sample_300ms = round(fs*0.3);
    sample_350ms = round(fs*0.35);
    sample_rem = ceil(prctile(rr_intervals_ecg,99)) - sample_350ms;

    pqrs_bloks = zeros(length(ecg_rpeaks_index)-2,sample_250ms+sample_70ms+1);
    t_bloks = zeros(length(ecg_rpeaks_index)-2,floor(sample_350ms*min(2,max(max(1,avg_intervals_ecg/(2*sample_350ms)))))-sample_70ms);
    remained_bloks = zeros(length(ecg_rpeaks_index)-2,sample_rem);
    size_T = size(t_bloks,2);
    size_R = size(remained_bloks,2);
    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_rpeaks_index(p)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1))):ecg_rpeaks_index(p)+sample_70ms;
        next_qrs_index = ecg_rpeaks_index(p+1)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p))):ecg_rpeaks_index(p+1)+sample_70ms;

        sample_T = size_T+sample_70ms;
        this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:min(next_qrs_index(1)-1,ecg_rpeaks_index(p)+sample_T);

        this_rem_index = this_t_index(end)+1:min(next_qrs_index(1)-1,this_t_index(end)+size_R);



        if any(this_qrs_index<1)||any(this_qrs_index>length(ecg_raw))
            continue;
        end

        if any(this_t_index>length(ecg_raw))
            continue;
        end

        if any(this_rem_index>length(ecg_raw))
            continue;
        end

        pqrs_bloks(p-1,end-length(this_qrs_index)+1:end) = ecg_raw(this_qrs_index) ;%- linspace(mean(ecg_raw(this_qrs_index(1:fs*0.01))), mean(ecg_raw(this_qrs_index(end-floor(fs*0.01):end))), 2*sample_70ms+1);
        pqrs_bloks(p-1,1:end-length(this_qrs_index)) = ecg_raw(this_qrs_index(1));

        temp = ecg_raw(this_t_index);
        if size_T>length(this_t_index) && isempty(this_rem_index)
            t_bloks(p-1,:) =  [temp,repelem(temp(end),1,size_T-length(this_t_index))];% - linspace(init_mn, mean(ecg_raw(this_t_index(end-floor(fs*0.015):end))), length(this_t_index));
            remained_bloks (p-1,:) = temp(end);
        elseif isempty(this_rem_index)
            t_bloks(p-1,:) =  temp;
            remained_bloks (p-1,:) = temp(end);
        else
            t_bloks(p-1,:) =  temp;
        end

        temp = ecg_raw(this_rem_index);
        if size_R>length(this_rem_index) && ~isempty(this_rem_index)
            remained_bloks(p-1,:) =  [temp,repelem(temp(end),1,size_R-length(this_rem_index))];% - linspace(init_mn, mean(ecg_raw(this_t_index(end-floor(fs*0.015):end))), length(this_t_index));
        elseif ~isempty(this_rem_index)
            remained_bloks (p-1,:) = temp;
        end


    end

    ecg_blocks = [pqrs_bloks,t_bloks,remained_bloks];

    qrs_bloks_mv = movmean(pqrs_bloks,[10,10],1,'omitmissing');
    t_bloks_mv = movmean(t_bloks,[10,10],1,'omitmissing');


    max_val = max(prctile(pqrs_bloks(:),99.5),prctile(t_bloks(:),99.5));
    min_val = min(prctile(pqrs_bloks(:),0.5),prctile(t_bloks(:),0.5));

    ecg_blocks(ecg_blocks>max_val) = max_val;
    ecg_blocks(ecg_blocks<min_val) = min_val;
    ecg_blocks = (ecg_blocks - min_val) / (max_val-min_val);

    fcm_features = [ecg_blocks(:,sample_70ms:2*sample_300ms),repmat(avg_intervals_ecg(2:end-1),1,3*sample_10ms),repmat(linspace(-1,1,size(ecg_blocks,1))',1,3*sample_10ms)];
    fcm_features = fillmissing(fcm_features,"linear");
    fcm_features = zscore(fcm_features);

    fcm_options = fcmOptions( MaxNumIteration=25, Exponent=1.1);
    [fcm_centers, fcm_part_mat] = fcm(fcm_features,fcm_options);
    [~,cluster_fcm] = max(fcm_part_mat);

    num_cls = size(fcm_centers,1);
    index_clustering = [];
    count_cls = zeros(num_cls,1);
    for c = 1:num_cls
        index_temp = find(cluster_fcm==c);
        count_cls(c) = length(index_temp);
        index_clustering = [index_clustering;[index_temp(:),c*ones(length(index_temp),1)]];
    end

    if length(cluster_fcm)<60
        index_clustering(:,2)=1;
    else
        index_clustering_org = index_clustering;
        for c = find(count_cls(:)'<30)
            [~,cls_max] = max(count_cls);
            index_clustering(index_clustering_org(:,2)==c,2) = cls_max;
        end
    end

    num_cls = unique(index_clustering(:,2));

    ecg_blocks_org = ecg_blocks;

    ecg_blocks = ecg_blocks(index_clustering(:,1),:);
    ecg_blocks_mv = movmean(ecg_blocks(:,sample_70ms:2*sample_300ms),[10,10],1,'omitmissing');


    ecg_noise = ecg_blocks(:,sample_70ms:2*sample_300ms) - ecg_blocks_mv;
    noise_power = movmean( sum(ecg_noise.^2,2),[30,30]);
    signal_power = movmean( sum(ecg_blocks_mv(:,1:sample_300ms+sample_70ms).^2,2),[30,30]);

    snr_signal_db = 10*log10(signal_power./noise_power);

    combining_weights = [exp((snr_signal_db-20)/10),exp((snr_signal_db-20)/10),exp((snr_signal_db-20)/10),exp((20-snr_signal_db)/10)];
    combining_weights = combining_weights./sum(combining_weights,2);

    ecg_noise = [ecg_noise,circshift(ecg_noise(:,end:-1:1),1,1),circshift(ecg_noise,1,5),circshift(ecg_noise(:,end:-1:1),1,3),circshift(ecg_noise,1,3)];
    ecg_noise = [ecg_noise,ecg_noise];
    ecg_noise = ecg_noise(:,1:size(ecg_blocks,2));

    y_est_bm3d = nan*ecg_blocks;
    y_est_curvelet = nan*ecg_blocks;
    y_est_adap_curvelet = nan*ecg_blocks;
    y_est_rf3d = nan*ecg_blocks;

    for c = num_cls(:)'

        ind_c = find(index_clustering(:,2)==c) ;
        PSD = abs(fft2(ecg_noise(ind_c,:), size(ecg_noise(ind_c,:),1), size(ecg_noise(ind_c,:),2))).^2;
        PSD = movmean(PSD,[3,3]);
        PSD = movmean(PSD,[3,3],2); % noise power spectrum

        y_est_bm3d(ind_c,:) = BM3D(ecg_blocks(ind_c,:), PSD ,'refilter');

        y_est_curvelet(ind_c,:) = ACT_filter(ecg_blocks(ind_c,:), PSD, 'ksigma');

        y_est_adap_curvelet(ind_c,:) = ACT_filter(ecg_blocks(ind_c,:), PSD);

        % RF3D parameters
        profile           = 'lc';  % Complexity profile:
        transform_2D_name = 'dct'; % 2D block transform used in RF3D (it also defines the PSDs)
        do_wiener         = false;  % Wiener Filtering (true or false).
        erf3d             = true;  % Enhanced FPN suppression (true or false).
        sigma_rnd = 1;            % Scaling factor of random noise (RND)
        sigma_fpn = 1;            % Scaling factor of fixed-pattern noise (FPN)

        y_est_rf3d(ind_c,:) = RF3D(ecg_blocks(ind_c,:), sigma_rnd, sigma_fpn, [], [], transform_2D_name, erf3d, do_wiener, profile);

    end

    y_est_adap_curvelet = y_est_adap_curvelet*(max_val-min_val) + min_val;
    y_est_curvelet = y_est_curvelet*(max_val-min_val) + min_val;
    y_est_bm3d = y_est_bm3d*(max_val-min_val) + min_val;
    y_est_rf3d = y_est_rf3d*(max_val-min_val) + min_val;

    % final combination
    y_est = combining_weights(:,1).*y_est_bm3d + combining_weights(:,2).*y_est_curvelet +  combining_weights(:,3).*y_est_adap_curvelet +  combining_weights(:,4).*y_est_rf3d;


    ecg_blocks(index_clustering(:,1),:) = ecg_blocks;
    y_est(index_clustering(:,1),:) = y_est;

    qrs_bloks_denoise = y_est(:,1:size(pqrs_bloks,2));
    t_bloks_denoise = y_est(:,size(pqrs_bloks,2)+1:size(pqrs_bloks,2)+size(t_bloks,2));
    remained_bloks_denoise = y_est(:,size(pqrs_bloks,2)+size(t_bloks,2)+1:end);
    %%
    P = size(qrs_bloks_denoise,1);
    for p = 2:length(ecg_rpeaks_index)-1

        this_qrs_index = ecg_rpeaks_index(p)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p-1))):ecg_rpeaks_index(p)+sample_70ms;
        next_qrs_index = ecg_rpeaks_index(p+1)-min(sample_250ms,floor(0.3*rr_intervals_ecg(p))):ecg_rpeaks_index(p+1)+sample_70ms;

        sample_T = size_T+sample_70ms;
        this_t_index = ecg_rpeaks_index(p)+sample_70ms+1:min(next_qrs_index(1)-1,ecg_rpeaks_index(p)+sample_T);

        this_rem_index = this_t_index(end)+1:min(next_qrs_index(1)-1,this_t_index(end)+size_R);



        if any(this_qrs_index<1)||any(this_qrs_index>length(ecg_raw))
            continue;
        end

        if any(this_t_index>length(ecg_raw))
            continue;
        end

        if any(this_rem_index>length(ecg_raw))
            continue;
        end

        ecg_denoised(this_qrs_index) = qrs_bloks_denoise(p-1,end-length(this_qrs_index)+1:end);

        ecg_denoised(this_t_index) = t_bloks_denoise(p-1,1:length(this_t_index));

        if ~isempty(this_rem_index)
            temp = remained_bloks_denoise(p-1,1:length(this_rem_index));
            ecg_denoised(this_rem_index) = temp - linspace(0,temp(end)-qrs_bloks_denoise(min(p,P),end-length(next_qrs_index)+1),length(temp));
        elseif length(this_rem_index) > size_R && ~isempty(this_rem_index)

            numrep = length(this_rem_index) - size_R;
            temp = [remained_bloks_denoise(p-1,:),repelem(remained_bloks_denoise(p-1,end),1,numrep)];
            ecg_denoised(this_rem_index) = temp - linspace(0,temp(end)-qrs_bloks_denoise(min(p,P),end-length(next_qrs_index)+1),length(temp));
        elseif isempty(this_rem_index)
            temp = t_bloks_denoise(p-1,1:length(this_t_index));
            ecg_denoised(this_t_index) = temp - linspace(0,temp(end)-qrs_bloks_denoise(min(p,P),end-length(next_qrs_index)+1),length(temp));
        end

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

