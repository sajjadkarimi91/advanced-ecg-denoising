
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



    for ch = 1:12

        index_R =  true_position(ch).R;
        ecg_rpeaks_index = index_R;

        ecg_denoised = ecg(:,ch)';
        ecg_denoised = ecg_denoised - movmean(movmedian(ecg_denoised,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);
        ecg_denoised = lp_filter_zero_phase(ecg_denoised, 30/fs);

        ecg_rpeaks_index_p = index_R;
        ecg_rpeaks_index_p(isnan(ecg_rpeaks_index_p))=[];

        ecg_QRSon_true = true_position(ch).QRSon; ecg_QRSon_wavdet_p = ecg_QRSon_true;  ecg_QRSon_wavdet_p(isnan(ecg_QRSon_wavdet_p))=[];


        ecg_QRSon_true = true_position(ch).QRSon; ecg_QRSon_wavdet_p = ecg_QRSon_true;  ecg_QRSon_wavdet_p(isnan(ecg_QRSon_wavdet_p))=[];

        ecg_QRSoff_true = true_position(ch).QRSoff; ecg_QRSoff_wavedet_p = ecg_QRSoff_true;  ecg_QRSoff_wavedet_p(isnan(ecg_QRSoff_wavedet_p))=[];

        ecg_Pon_index = true_position(ch).Pon;  ecg_Pon_index_p = ecg_Pon_index; ecg_Pon_index_p(isnan(ecg_Pon_index_p))=[];
        ecg_P_index = true_position(ch).P;  ecg_P_index_p = ecg_P_index; ecg_P_index_p(isnan(ecg_P_index_p))=[];
        ecg_Poff_index = true_position(ch).Poff;  ecg_Poff_index_p = ecg_Poff_index; ecg_Poff_index_p(isnan(ecg_Poff_index_p))=[];

        ecg_Ton_true = true_position(ch).Ton';     ecg_Ton_wavedet_p = ecg_Ton_true;  ecg_Ton_wavedet_p(isnan(ecg_Ton_wavedet_p))=[];

        ecg_T_true = true_position(ch).T';     ecg_T_wavedet_p = ecg_T_true;  ecg_T_wavedet_p(isnan(ecg_T_wavedet_p))=[];

        ecg_Toff_true = true_position(ch).Toff';  ecg_Toff_wavedet_p = ecg_Toff_true;  ecg_Toff_wavedet_p(isnan(ecg_Toff_wavedet_p))=[];


        a = figure('Position', [130 130 1500 800]);
        lg = {};
        plot(t_second,ecg_denoised,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
        hold on
        ecg_plot = ecg_denoised;
        plot(t_second(ecg_rpeaks_index_p),ecg_plot(ecg_rpeaks_index_p),'*',MarkerSize=6,LineWidth=1);lg = cat(2, lg, {'R'});

        if ~isempty(ecg_Ton_wavedet_p)
            plot(t_second(ecg_Ton_wavedet_p),ecg_plot(ecg_Ton_wavedet_p),'r+',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-Ton'});
        end
        if ~isempty(ecg_T_wavedet_p)
            plot(t_second(ecg_T_wavedet_p),ecg_plot(ecg_T_wavedet_p),'r*',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-T'});
        end

        if ~isempty(ecg_Toff_wavedet_p)
            plot(t_second(ecg_Toff_wavedet_p),ecg_plot(ecg_Toff_wavedet_p),'rx',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-Toff'});
        end

        if ~isempty(ecg_Pon_index_p)
            plot(t_second(ecg_Pon_index_p),ecg_plot(ecg_Pon_index_p),'k+',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'Pon'});
            plot(t_second(ecg_P_index_p),ecg_plot(ecg_P_index_p),'k*',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'P'});
            plot(t_second(ecg_Poff_index_p),ecg_plot(ecg_Poff_index_p),'kx',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'Poff'});
        end

        plot(t_second(ecg_QRSon_wavdet_p),ecg_plot(ecg_QRSon_wavdet_p),'go',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-QRSon'});
        plot(t_second(ecg_QRSoff_wavedet_p),ecg_plot(ecg_QRSoff_wavedet_p),'go',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-QRSoff'});

        grid on
        legend(lg,'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
        xlabel('time (sec)',Interpreter='latex',FontSize=14)

        pause()
    end
end



