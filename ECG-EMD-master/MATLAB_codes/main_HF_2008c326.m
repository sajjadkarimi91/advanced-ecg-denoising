close all
%load ecg signal
load('m103.mat');
Fs=360;
gain=200;
T=10; % 10 seconds
ecg_sig=(val(1:T*Fs)-1024);
sig_len=length(ecg_sig);
% t=linspace(0,T,T*Fs);
t=1:sig_len;


%% Generate noise
% sig_power=(1/sig_len)*(norm(ecg_sig).^2);
sig_power=var(ecg_sig);
SNR=10;
Nstd=sqrt(sig_power/(10^(SNR/10)));

SNR=10;
Nstd_emg=sqrt(sig_power/(10^(SNR/10)));
emg_syn=Nstd_emg*randn([1,sig_len]);
emg_syn(600:1200)=0;
emg_syn(2400:3000)=0;

%% Construct noisy signal
ecg_in=ecg_sig+emg_syn;
% ecg_in=ecg_sig;


%% R peak detection (fiducial point)

amp_thresh=0.4*max(ecg_in);


sig_clip=zeros(1,sig_len);

for i=1:sig_len
    if ecg_in(i)>amp_thresh
        sig_clip(i)=ecg_in(i);
    else
        sig_clip(i)=amp_thresh;
    end
end
sig_clip=sig_clip/gain;
sig_diff= diff(sig_clip);


qrs_win_len=floor(Fs*0.12); %120ms qrs size;
deriv_thresh=0.02*(360/Fs);

qrs_window=zeros(1,sig_len);
qrs_onset_list=[];
r_index=[];
i=1;
j=1;
while i<sig_len
    if sig_diff(i)>deriv_thresh ...
            && sig_diff(i+1)>deriv_thresh ...
            %             && sig_diff(i+2)>deriv_thresh
        qrs_onset_list(j)=i;
        
        ind_till=min(i+qrs_win_len,sig_len);
        qrs_window(i:ind_till)=1;
        [~,r_index(j)]=max(ecg_in(i:ind_till));
        r_index(j)=r_index(j)+i-1;
        
        j=j+1;
        i=i+qrs_win_len;
    else
        i=i+1;
    end
end

num_qrs=numel(r_index);
% plot(sig_clip)
r_positions=zeros(1,sig_len);
r_positions(r_index)=400;

figure
plot(t,400*qrs_window,t,ecg_sig,t,r_positions);
legend('qrs window','ecg signal','R positions')
% plot(ecg_sig)


%% Find the IMFs
NR = 80;
MaxIter = 500;
tic
[modes ~]=ceemdan_v2014(ecg_in,Nstd,NR,MaxIter,2);
% modes=emd(ecg_in);

t1=toc;
t=1:length(ecg_in);

[a b]=size(modes);

figure;
subplot(a+1,1,1);
plot(t,ecg_in);% the ECG signal is in the first row of the subplot
% title('CEEMDAN Decomposition for ECG signal')
ylabel('ECG')
set(gca,'xtick',[])
ylim([min(ecg_in),max(ecg_in)])
set(gca,'ytick',[])
axis tight;

for i=2:a
    subplot(a+1,1,i);
    plot(t,modes(i-1,:));
    ylabel ([num2str(i-1)]);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([1 length(ecg_in)])   
    ylim([min(modes(i-1,:)),max(modes(i-1,:))])
end;

subplot(a+1,1,a+1)
plot(t,modes(a,:))
ylabel([num2str(a)])
xlim([1 length(ecg_in)])
ylim([min(modes(a,:)),max(modes(a,:))])
set(gca,'ytick',[])
%% find QRS

first3imfSum=sum(modes(1:3,:),1);
% plot(first3imfSum)
%%
right_zc=zeros(1,num_qrs);
left_zc=zeros(1,num_qrs);

zero_crossings=diff(sign(first3imfSum));
% plot(t,[zero_crossings,0],t,first3imfSum/200);

k=1;
for i= r_index
    j=0;
    while zero_crossings(i+j)<=0
        j=j+1;
    end
    right_zc(k)=i+j;
    k=k+1;
end

k=1;
for i= r_index
    j=0;
    while zero_crossings(i-j)>=0
        j=j+1;
    end
    left_zc(k)=i-j;
    k=k+1;
end

qrs_box=zeros(1,sig_len);
for i=1:num_qrs
    qrs_box(left_zc(i):right_zc(i))=1;
end

% plot(t,(max(ecg_sig)+abs(min(ecg_sig)))*qrs_box+min(ecg_sig),t,first3imfSum);
figure
plot(t,first3imfSum,t,ecg_sig);
legend('sum of c_1 to c_3', 'ECG signal')
grid on
% xL = xlim;
% yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis
%% t-test
alpha=0.1;
P=t_test(alpha,modes)% P=filter order

%%
qrs_box_len=right_zc-left_zc;
qrs_centre=round((right_zc+left_zc)/2);
P=4;
a=0.1;
low_order_sum=zeros(1,sig_len);
low_order_sum_cmpl=zeros(1,sig_len);
for j=1:P
    beta=j*0.3;
    qrs_tukey = apply_tukey(right_zc,left_zc,beta,sig_len);
    qrs_tukey_cmpl=1-qrs_tukey;
    low_order_sum=low_order_sum+ qrs_tukey.*modes(j,:);
    low_order_sum_cmpl=low_order_sum_cmpl+ a*qrs_tukey_cmpl.*modes(j,:);
end

high_order_sum=sum(modes(P+1:end,:),1);

recons_sig=low_order_sum+low_order_sum_cmpl+high_order_sum;

% Plots
% figure
% subplot(3,1,1)
% plot(t,ecg_sig)
% 
% subplot(3,1,2)
% plot(t,recons_sig)
% 
% subplot(3,1,3)
% plot(t,ecg_in)

%% SER calculation

err_power=norm(ecg_sig-ecg_in)^2;
sig_power= norm(ecg_sig)^2;
SER=sig_power/err_power;
SER_dB=10*log10(SER)

%%
qrs_tukey = apply_tukey(right_zc,left_zc,0.3,sig_len);
% plot(t,(max(ecg_sig)+abs(min(ecg_sig)))*qrs_box+min(ecg_sig),t,ecg_sig);
% plot(t,qrs_tukey)
% ylim([-0.2,1.2])
ceemdan_recons = recons_sig;
% save ceemdan_recons

%%
% Plots
figure
subplot(3,1,1)
plot(t,ecg_in)
title('(a) ECG with intermittent noise (10 dB SNR)')
ylim([min(ecg_in), max(ecg_in)]);
xlim([1,length(ecg_in)]);


subplot(3,1,2)
plot(t,ceemdan_recons)
title('(b) Reconstruction using improved CEEMDAN')
ylim([min(ceemdan_recons), max(ceemdan_recons)]);
xlim([1,length(ceemdan_recons)]);

subplot(3,1,3)
plot(t,emd_recons)
title('(c) Reconstruction using EMD')
ylim([min(emd_recons), max(emd_recons)]);
xlim([1,length(emd_recons)]);

% plot(t,ecg_sig)












