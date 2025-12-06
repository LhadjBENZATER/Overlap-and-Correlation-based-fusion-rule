clear; clc; close all; set(0, 'DefaultAxesFontSize', 14); rng(0)
N = 1024; M = round(0.25*N); a=0.7; bk=0.65; Psi = fft(eye(N)) / sqrt(N);     
Phi = Phi_selection(N,M,2,2e7); A=Phi/Psi; At=A';th=0.1;
%[x, xf, ~, f, K, index, ~] = WFG_Re(9e8, 2e7, N, 20); close figure 1;
[x,xf,index,K,f,t,ind_shft]=LFM(2.9e9,4e7,N,1);%close figure 1;
MC=1; L=8; Pd=zeros(5,L-2); Pf=Pd; RMSE=Pd; Y=Phi*x; tic;ncor=zeros(1,L-2);
for l = 3:L
disp(['Network size (L)', ' = ', num2str(l)]);
Xf_mat = zeros(N,2*l+1); P_d = zeros(5, 1); P_f = P_d;RMSE_loc=P_d;
for frame = 1:MC
%[x1,x2,x3,x4]=Read_USRP(frame);[x5,x6,x7,x8]=Read_USRP(frame+300);
[x1,x2,x3,x4,x5,x6,x7,x8]=Read_LFM(frame);
x=[x1 x2 x3 x4 x5 x6 x7 x8];x = awgn(x, 10, 'measured');
x=fct_syn(x); x=fct_Phase(x);
y=Phi*x; YL=Virtual(y(:,1:l));          
for i = 1:2 * l + 1
    Xf_mat(:,i) = fftshift(abs(fct_MRMP(A,At,YL(:,i),K,bk,a)));
end
%% Geometric Mean
xf_gm=prod(Xf_mat(:,l+1:2*l),2).^(1/l);
mx = max(xf_gm);mx(mx==0)=1;xf_gm=abs(xf_gm)./mx;
P_d(1)=P_d(1)+sum(xf_gm(index)>th)/K;
P_f(1)=P_f(1)+(sum(xf_gm>th)-sum(xf_gm(index)>th))/(N-K);
RMSE_loc(1)=RMSE_loc(1)+sqrt(mean((xf-xf_gm).^2, 1)); 
%% Arithmetic Mean
xf_m=mean(Xf_mat,2);xf_m=xf_m-median(xf_m); xf_m=xf_m./max(xf_m);
P_d(2)=P_d(2)+sum(xf_m(index)>th)/K;
P_f(2)=P_f(2)+(sum(xf_m>th)-sum(xf_m(index)>th))/(N-K);
RMSE_loc(2)=RMSE_loc(2)+sqrt(mean((xf-xf_m).^2, 1)); 
%% Combined by Corr th=0.1
ncor(l-2)=sum(xf_m.*xf_gm)/(norm(xf_m)*norm(xf_gm));
xf_comb = (1-ncor(l-2))*xf_gm+ncor(l-2)*xf_m;%xf_comb=xf_comb-median(xf_comb); 
xf_comb=xf_comb/max(xf_comb);
P_d(3)=P_d(3)+sum(xf_comb(index)>th)/K;
P_f(3)=P_f(3)+(sum(xf_comb>th)-sum(xf_comb(index)>th))/(N-K);
RMSE_loc(3)=RMSE_loc(3)+sqrt(mean((xf-xf_comb).^2,1)); 
%s=sort(xf_comb,'descend');

%% Combined overlapping
xf_comb = comb(xf_m,xf_gm,th);
P_d(4)=P_d(4)+sum(xf_comb(index)>th)/K;
P_f(4)=P_f(4)+(sum(xf_comb>th)-sum(xf_comb(index)>th))/(N-K);
RMSE_loc(4)=RMSE_loc(4)+sqrt(mean((xf-xf_comb).^2,1)); 

end
Pd(:,l-2)=P_d; Pf(:,l-2)=P_f; RMSE(:,l-2)=RMSE_loc;
end    
Pd=Pd/MC; Pf=Pf/MC; RMSE=RMSE/MC;
disp(['Time = ', datestr(now, 'HH:MM:SS'), ', Duration(min)', ' = ', num2str(toc / 60)]);
subplot(2, 2, 1); 
plot(3:L,ncor,'LineWidth',1.2);grid on;title('Normalized correlation');
xlabel('Network size');ylabel('Normalized correlation');
subplot(2,2,2); plot(3:8,Pd,'LineWidth',1.5); grid on;ylabel('Q_d');title('Q_d');
legend('Geometric Mean','Arhithmetic Mean','Corr Comb','Overlap Comb');xlabel('Network size'); 
subplot(2,2,3); semilogy(3:8,Pf,'LineWidth',1.5); grid on;title('Q_f');ylabel('Q_f');
legend('Geometric Mean','Arhithmetic Mean','Corr Comb','Overlap Comb');xlabel('Network size'); 
subplot(2,2,4); semilogy(3:8,RMSE,'LineWidth',1.5); grid on;ylabel('RMSE');title('RMSE');
legend('Geometric Mean','Arhithmetic Mean','Corr Comb','Overlap Comb');xlabel('Network size'); 

figure
xf1=fftshift(abs(fft(x1)))/max(fftshift(abs(fft(x1))));
subplot(2,2,1);plot(f,xf1,'LineWidth',1.5,'LineStyle','--'); grid on;
title('Original Spectrum');xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');
subplot(2,2,2);plot(f,xf_m,'LineWidth',1.5,'LineStyle','--'); grid on;
title('Arithmetic Avg Spectrum');xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');
subplot(2,2,3);plot(f,xf_gm,'LineWidth',1.5,'LineStyle','-.'); grid on;
title('Geometric Avg Spectrum');xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');
subplot(2,2,4);plot(f,xf_comb,'LineWidth',1.5,'LineStyle',':'); grid on;
title('Combined Avg Spectrum');xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');


