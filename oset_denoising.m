clear;
close all;
clc;


oset_path = 'D:\projects\toolboxes\OSET'; 'D:\my_projects\Alphanumeric\multimodal-cardiac-biomarkers\OSET';% enter the path of OSET toolbox
addpath(genpath(oset_path))

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

        ecg_in = ecg_noisy(:,ch)';


        % ecg_in = ecg_den_lti_smoother(ecg_in,2,20);
        % ecg_in = ecg_in - movmean(movmedian(ecg_in,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);


        params = [];
        ecg_rpeaks = zeros(1,length(ecg_in));
        ecg_rpeaks(ecg_rpeaks_index) = 1;
        % % [ecg_in, data_prior_est, n_var] = ecg_den_phase_domain_gp(ecg_in, ecg_rpeaks, params);
        % % [ecg_in, data_prior_est, n_var] = ecg_den_time_domain_gp(ecg_in, ecg_rpeaks, params);
        % % ecg_in = ecg_in - movmean(movmedian(ecg_in,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);


        % [ecg_in, x_filtered1] = ecg_den_seg_wise_smoother(ecg_in,2, 4*fs, 10);
        % ecg_in = ecg_in - movmean(movmedian(ecg_in,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);


        params.TPTR =  {'sqtwolog'};% {'rigrsure', 'heursure', 'sqtwolog', 'minimaxi'};
        params.SCAL = {'mln'}; % {'one', 'sln', 'mln'};
        params.SORH = {'s'};
        params.WLEVELS = [5,6];
        [y_med, ecg_in] = ecg_den_wavelet(ecg_in, params);
        [y_med, ecg_in] = ecg_den_wavelet(ecg_in, params);
        ecg_in = ecg_in - movmean(movmedian(ecg_in,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);

        ecg_clean = ecg(:,ch)' ;
        ecg_recomb = ecg_in;

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
