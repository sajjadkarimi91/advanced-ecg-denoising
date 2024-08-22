
clear;
close all;
clc;


oset_path = 'D:\projects\toolboxes\OSET'; 'D:\my_projects\Alphanumeric\multimodal-cardiac-biomarkers\OSET';% enter the path of OSET toolbox
addpath(genpath(oset_path))

db_folder = '.\localLU';
local_db_files = dir([db_folder '/*.mat']); % list of all mat files
db_folder_save = '.\localLU_noisy';
f50 = 50;

for m = 1:length(local_db_files)

    clc
    close all

    disp(m)

    % load data and rpeaks
    in_fname = local_db_files(m).name(1:end-4);

    load([db_folder,'/',in_fname,'.mat']);

    ecg_noisy = ecg;
    for ch = 1:12

        index_R =  true_position(ch).R;
        ecg_rpeaks_index = index_R;

        ecg_denoised = ecg(:,ch)';
        ecg_denoised = ecg_denoised - movmean(movmedian(ecg_denoised,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);

        baseline_wander = randn(1,length(ecg_denoised)+5*fs);
        baseline_wander = lp_filter_zero_phase(baseline_wander, 0.1/fs);
        baseline_wander = baseline_wander(5*fs+1:end);
        baseline_wander = baseline_wander * std(ecg_denoised)/std(baseline_wander);

        white_noise = randn(1,length(ecg_denoised));
        white_noise = white_noise * std(ecg_denoised)/std(white_noise);

        power_line = sin(2*pi*f50/fs*(1:length(ecg_denoised)+2*fs)+2*pi*rand(1));
        baseline_50 = randn(1,length(ecg_denoised)+2*fs);
        baseline_50 = lp_filter_zero_phase(baseline_50, 0.5/fs);
        power_line = power_line.* abs(baseline_50);
        power_line = power_line(2*fs+1:end);
        power_line = power_line * std(ecg_denoised)/std(power_line);


        emg_signal = randn(1,length(ecg_denoised)+2*fs);
        emg_signal = sjk_eeg_filter(emg_signal, fs,20+randi(10),100+randi(50));

        % emg_signal = lp_filter_zero_phase(emg_signal, 100/fs);
        % emg_signal = emg_signal - lp_filter_zero_phase(emg_signal, 25/fs);
        baseline_50 = randn(1,length(ecg_denoised)+2*fs);
        baseline_50 = lp_filter_zero_phase(baseline_50, 0.2/fs);
        emg_signal = emg_signal.* abs(baseline_50);
        emg_signal = emg_signal(2*fs+1:end);
        emg_signal = emg_signal * std(ecg_denoised)/std(emg_signal);

        ecg_noisy(:,ch) = ecg_denoised + 0.25*emg_signal + 0.25*power_line + 0.1*white_noise + baseline_wander;
        ecg(:,ch) = ecg_denoised;

    end

    save( [db_folder_save,'\',in_fname,'.mat'], 'ecg_noisy','ecg','t_second', 'true_position','fs')
end



