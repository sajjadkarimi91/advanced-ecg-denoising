
clear;
close all;
clc;


oset_path = 'D:\projects\toolboxes\OSET'; 'D:\my_projects\Alphanumeric\multimodal-cardiac-biomarkers\OSET';% enter the path of OSET toolbox
addpath(genpath(oset_path))

db_folder = '.\localLU';
local_db_files = dir([db_folder '/*.mat']); % list of all mat files

for m = 1:length(local_db_files)

    clc
    close all

    disp(m)

    % load data and rpeaks
    in_fname = local_db_files(m).name(1:end-4);

    load([db_folder,'/',in_fname,'.mat']);

    a = figure('Position', [130 130 1500 800]);
    for i = 1:12
        subplot(6,2,i)
        plot(t_second,ecg(:,i),LineWidth=1.5)
        grid on
        if i==11 || i==12
            xlabel('time (sec)',Interpreter='latex',FontSize=16)
        end
    end

    for ch = 1:12

        index_R =  true_position(ch).R;
        ecg_rpeaks_index = index_R;

        ecg_denoised = ecg(:,ch)';
        ecg_denoised = ecg_denoised - movmean(movmedian(ecg_denoised,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);
        ecg_denoised = lp_filter_zero_phase(ecg_denoised, 30/fs);

        ecg_rpeaks_index_p = index_R;
        ecg_rpeaks_index_p(isnan(ecg_rpeaks_index_p))=[];


        a = figure('Position', [130 130 1500 800]);
        lg = {};
        plot(t_second,ecg_denoised,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
        hold on
        ecg_plot = ecg_denoised;
        plot(t_second(ecg_rpeaks_index_p),ecg_plot(ecg_rpeaks_index_p),'*',MarkerSize=6,LineWidth=1);lg = cat(2, lg, {'R'});

        grid on
        legend(lg,'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
        xlabel('time (sec)',Interpreter='latex',FontSize=14)

        pause()
    end


end



