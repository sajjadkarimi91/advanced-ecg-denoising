clear;
close all;
clc;

% https://github.com/shloked/Denoise-ECG-using-EMD
emd_path = '.\ECG-EMD-master';% enter the path of EMD toolbox
addpath(genpath(emd_path))

db_folder = '.\localLU_noisy';
local_db_files = dir([db_folder '/*.mat']); % list of all mat files


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


        %% Find the IMFs

        sig_power=(var(ecg_in))/2;
        SNR=0;%dB SNR=5
        Nstd=sqrt(sig_power/(10^(SNR/10)));
        T = length(ecg_in)/fs;
        NR = 80;
        MaxIter = 500;
        tic
        [modes, ~]=ceemdan_v2014(ecg_in,Nstd,NR,MaxIter,2);
        % modes=emd(ecg_sig);

        % ZCR
        zcr=zeros(size(modes,1),1);
        for i=1:size(modes,1)
            zc=abs(diff(sign(modes(i,:))));
            zc_nos=length(zc(zc==2));
            zcr(i)=zc_nos/T;
        end

        % recombine
        ecg_recomb=zeros(1,length(ecg_in));
        bw_est=zeros(1,length(ecg_in));
        bw_noise=zeros(1,length(ecg_in));
        for i=1:size(modes,1)
            if zcr(i)>1.5 && zcr(i)<50
                ecg_recomb=ecg_recomb+modes(i,:);
            elseif zcr(i)<1.5
                bw_est=bw_est+modes(i,:);
            elseif zcr(i)>=50
                bw_noise=bw_noise+modes(i,:);
            end
        end


        ecg_clean = ecg(:,ch)' ;

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
