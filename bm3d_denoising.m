clear;
close all;
clc;


bm3d_path = './BM3D';% enter the path of BM3D toolbox
addpath(genpath(bm3d_path))

db_folder = '.\localLU_noisy';
local_db_files = dir([db_folder '/*.mat']); % list of all mat files

f50 = 50;

for m = 1:length(local_db_files)

    clc
    close all

    disp(m)

    % load data and rpeaks
    in_fname = local_db_files(m).name(1:end-4);

    load([db_folder,'/',in_fname,'.mat']);

    % ecg_noisy = ecg;
    for ch = 1:12

        index_R =  true_position(ch).R;
        ecg_rpeaks_index = index_R;

        ecg_raw = ecg_noisy(:,ch)';

ecg_pwelch = ecg_raw - movmean(movmedian(ecg_raw,[floor(0.3*fs),floor(0.3*fs)]),[floor(0.15*fs),floor(0.15*fs)]);
        [Pxx,freq] = pwelch(ecg_pwelch,[],[],[],fs);
        fc = 50.0; % powerline frequency
        if sum(Pxx(freq>fc-1 & freq<fc+1 ))/sum(Pxx(freq<fc-1&freq>2))> 0.05 || sum(Pxx(freq>2*fc-1 & freq<2*fc+1 ))/sum(Pxx(freq<fc-1&freq>1))> 0.05 || sum(Pxx(freq>fc-1 & freq<fc+1 ))/sum(Pxx(freq<fc-1&freq>20))> 0.1
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



        peak_detector_params.RETURN_SIGNAL_PEAKS = true; % return signal peaks or energy envelope peaks
        peak_detector_params.PLOT_RESULTS = false; % plot the results using the internal plot function of peak_det_likelihood or not
        peak_detector_params.PLOT_DIAGNOSTIC = false; % diagnostic mode (do not activate unless diving deep into the code! run only on short segments, since many figures are created)
        peak_detector_params.verbose = false; % reports all the default values for the internal parameters of peak_det_likelihood, which can be modified through this data structure if needed.
        peak_detector_params.REFINE_PEAKS = true;

        overlap_time = 1.0; % overlap between segments for continuity (1.0-2.0 seconds is enough)
        seg_len_time = 10.0; % segment length in seconds

        [peaks, ecg_rpeaks_index, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(lp_filter_zero_phase(ecg_raw, 30/fs), fs, seg_len_time, overlap_time, peak_detector_params);
        ecg_rpeaks_index = peak_indexes_consensus;
        ecg_rpeaks_index =ecg_rpeaks_index(:);

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

            temp_block = ecg_blocks(ind_c,:);
            temp_noise = ecg_noise(ind_c,:);

            while size(temp_block,1)<8
                temp_block = [temp_block;temp_block(randperm(length(ind_c)),:)];
                temp_noise = [temp_noise;temp_noise(randperm(length(ind_c)),:)];
            end

            PSD = abs(fft2(temp_noise, size(temp_noise,1), size(temp_noise,2))).^2;
            PSD = movmean(PSD,[3,3]);
            PSD = movmean(PSD,[3,3],2); % noise power spectrum


            temp = BM3D(temp_block, PSD ,'refilter');
            temp = BM3D(temp, PSD ,'refilter');
            y_est_bm3d(ind_c,:) = temp(1:length(ind_c),:);

            % % % temp =  ACT_filter(temp_block, PSD, 'ksigma');
            % % % y_est_curvelet(ind_c,:) = temp(1:length(ind_c),:);
            % % % temp =  ACT_filter(temp_block, PSD);
            % % % y_est_adap_curvelet(ind_c,:) = temp(1:length(ind_c),:);
            % % % 
            % % % % RF3D parameters
            % % % profile           = 'lc';  % Complexity profile:
            % % % transform_2D_name = 'dct'; % 2D block transform used in RF3D (it also defines the PSDs)
            % % % do_wiener         = false;  % Wiener Filtering (true or false).
            % % % erf3d             = true;  % Enhanced FPN suppression (true or false).
            % % % sigma_rnd = 1;            % Scaling factor of random noise (RND)
            % % % sigma_fpn = 1;            % Scaling factor of fixed-pattern noise (FPN)
            % % % 
            % % % temp = RF3D(temp_block, sigma_rnd, sigma_fpn, [], [], transform_2D_name, erf3d, do_wiener, profile);
            % % % y_est_rf3d(ind_c,:) = temp(1:length(ind_c),:);

        end

        % % % y_est_adap_curvelet = y_est_adap_curvelet*(max_val-min_val) + min_val;
        % % % y_est_curvelet = y_est_curvelet*(max_val-min_val) + min_val;
        y_est_bm3d = y_est_bm3d*(max_val-min_val) + min_val;
        % % % y_est_rf3d = y_est_rf3d*(max_val-min_val) + min_val;

        % final combination
        % % % y_est = combining_weights(:,1).*y_est_bm3d + combining_weights(:,2).*y_est_curvelet +  combining_weights(:,3).*y_est_adap_curvelet +  combining_weights(:,4).*y_est_rf3d;

        y_est = y_est_bm3d;

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

        ecg_clean = ecg(:,ch)' ;
ecg_recomb = ecg_denoised;

        figure
        plot(t_second,ecg_noisy(:,ch))
        hold on
        plot(t_second,ecg_recomb,LineWidth=1.5)
        plot(t_second,ecg_clean,LineWidth=1.5)
        grid on
        legend({'Noisy','Denoised','Clean'},'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
        xlabel('time (sec)',Interpreter='latex',FontSize=14)

        figure
        plot(t_second,ecg_clean,LineWidth=1.5,Color='#EDB120')
        hold on
        plot(t_second,ecg_recomb,LineWidth=1.5)
        grid on
        legend({'Clean','Denoised'},'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
        xlabel('time (sec)',Interpreter='latex',FontSize=14)
        xlim([6,7.5])

    end


end
