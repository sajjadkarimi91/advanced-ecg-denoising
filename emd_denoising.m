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


%% Find the IMFs
NR = 80;
MaxIter = 500;
tic
[modes, ~]=ceemdan_v2014(ecg_sig,Nstd,NR,MaxIter,2);
% modes=emd(ecg_sig);





%% ZCR
zcr=zeros(size(modes,1),1);
for i=1:size(modes,1)
    zc=abs(diff(sign(modes(i,:))));
    zc_nos=length(zc(zc==2));
    zcr(i)=zc_nos/T;
end

%% recombine
ecg_recomb=zeros(1,length(ecg_sig_uncorr));
bw_est=zeros(1,length(ecg_sig_uncorr));
for i=1:size(modes,1)
    if zcr(i)>1.5
        ecg_recomb=ecg_recomb+modes(i,:);
    else
        bw_est=bw_est+modes(i,:);
    end
end


        ecg_denoised = ecg(:,ch)' ;

        figure
        plot(t_second,ecg_noisy(:,ch))
        hold on
        plot(t_second,ecg_in)
        plot(t_second,ecg_denoised)

    end


end
