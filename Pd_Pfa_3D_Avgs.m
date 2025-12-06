%% Methods weighting
clear; clc; close all; set(0,'DefaultAxesFontSize',14); %rng(0)
N=1024; M=round(0.25*N); a=0.7; bk=0.65; Psi=fft(eye(N))/sqrt(N);
%[x, xf, ~, f,K,index,ind_shft] = WFG_Re(2.4e9,4e7,N,20); 
W=[0 1 0 1 1 0 0];[x,xf,index,K,f,t,ind_shft]=Sig(2.9e9,4e7,N,1,W);
%[x,xf,index,K,f,t,ind_shft]=LFM(2.9e9,4e7,N,1);
%Doub_to_bin_DAT(real(x),'x_LFM.dat'); %Already saved
Phi=Phi_selection(N,M,2,2e7); A=Phi/Psi; At=A'; Y = Phi*x;
th=0.05; L=[4,8]; SNR=-20:40:20; MC=1;alphas = 0:0.01:1;
iter = length(alphas);
l=length(L); snr=length(SNR); Pd=zeros(3,3,l,snr); Pf=Pd; RMSE=Pd;  
tic;ncor=zeros(l,snr);  
for s = 1:snr
fprintf('SNR(dB) = %d\n', SNR(s));  
Pd_loc=zeros(3,3,l); Pf_loc=Pd_loc; RMSE_loc=Pd_loc;
 for j = 1:l
  for frame = 1:MC
  y = Virtual(diff_chn(Y, L(j), SNR(s), 0));     
%  y = diff_chn(Y, L(j), SNR(s), 0);     
   Xf_mat = zeros(N, 2*L(j)+1);
   for i = 1:2*L(j)+1
    Xf_mat(:,i)=fftshift(abs(fct_MRMP(A, At, y(:,i),K,bk,a)));
   end
%Xf_mat(:,1:L(j))=fftshift(abs(fct_MRMP(A, At, y(:,1:L(j)),K,bk,a)));
%Xf_mat(:,L(j)+1:2*L(j)+1)=fftshift(abs(fct_MRMP(A, At, y(:,L(j)+1:2*L(j)+1),K,bk,a)));
%fct_MRMP_Virt
n = max(Xf_mat,[],1);n(n==0)=1;Xf_mat=abs(Xf_mat)./n;
w = [ones(1, L(j)),sqrt(L(j)-1)*ones(1, L(j)+1)];w=w/max(w);

%% Geometric Mean
xf_gm=[prod(Xf_mat(:,L(j)+1:2*L(j)),2).^(1/L(j)),...
   prod(Xf_mat,2).^(1/(2*L(j)+1)) prod(Xf_mat(:,1:L(j)),2).^(1/L(j))];
mx = max(xf_gm,[],1);mx(mx==0)=1;xf_gm=abs(xf_gm)./mx;
Pd_loc(1,:,j)=Pd_loc(1,:,j)+sum(xf_gm(index,:)>th)/K;
Pf_loc(1,:,j)=Pf_loc(1,:,j)+(sum(xf_gm>th)-sum(xf_gm(index,:)>th))/(N-K);
RMSE_loc(1,:,j)=RMSE_loc(1,:,j)+sqrt(mean((xf-xf_gm).^2, 1)); 

%% Mean
xf_m=[mean(Xf_mat(:,L(j)+1:2*L(j)),2) mean(Xf_mat,2) Xf_mat*w'];
xf_m=xf_m-median(xf_m);xf_m(xf_m<0)=0;xf_m=xf_m./ max(xf_m,[],1);
Pd_loc(2,:,j)=Pd_loc(2,:,j)+sum(xf_m(index,:)>th)/K;
Pf_loc(2,:,j)=Pf_loc(2,:,j)+(sum(xf_m>th)-sum(xf_m(index,:)>th))/(N-K);
RMSE_loc(2,:,j)=RMSE_loc(2,:,j)+sqrt(mean((xf-xf_m).^2, 1));
%% Alpha comb
Pd_comb=zeros(iter,1);Pf_comb=Pd_comb;RMSE_comb=Pd_comb;
for idx=1:iter
alpha=alphas(idx);
xf_comb=alpha*xf_gm(:,1)+(1-alpha)*xf_m(:,3); xf_comb=xf_comb/max(xf_comb);
Pd_comb(idx)=sum(xf_comb(index)>th)/K; 
Pf_comb(idx)=(sum(xf_comb>th)-sum(xf_comb(index)>th))/(N-K); 
RMSE_comb(idx)=sqrt(mean((xf-xf_comb).^2)); 
end
Obj=Pd_comb-10*(Pf_comb+RMSE_comb);[~,opt_idx]= max(Obj);opt_alpha = alphas(opt_idx);
xf_comb = opt_alpha*xf_gm(:,1)+(1-opt_alpha)*xf_m(:,3); xf_comb=xf_comb/max(xf_comb);
Pd_loc(3,1,j)=Pd_loc(3,1,j)+sum(xf_comb(index)>th)/K;
Pf_loc(3,1,j)=Pf_loc(3,1,j)+(sum(xf_comb>th)-sum(xf_comb(index)>th))/(N-K);
RMSE_loc(3,1,j)=RMSE_loc(3,1,j)+sqrt(mean((xf-xf_comb).^2));
%% Corr Comb
Xf_gm = max(xf_gm(:,[1,3]),[],2);
ncor(j,s)=sum(xf_m(:,2).*Xf_gm)/(norm(xf_m(:,2))*norm(Xf_gm));
xf_comb2 = (1-ncor(j,s))*Xf_gm+ncor(j,s)*xf_m(:,3);%xf_comb2=xf_comb2-median(xf_comb2); 
xf_comb2=xf_comb2/max(xf_comb2);
Pd_loc(3,2,j)=Pd_loc(3,2,j)+sum(xf_comb2(index)>th)/K;
Pf_loc(3,2,j)=Pf_loc(3,2,j)+(sum(xf_comb2>th)-sum(xf_comb2(index)>th))/(N-K);
RMSE_loc(3,2,j)=RMSE_loc(3,2,j)+sqrt(mean((xf-xf_comb2).^2));
%% Overlap Comb
xf_comb3 = comb(xf_m(:,3),Xf_gm, th);
Pd_loc(3,3,j)=Pd_loc(3,3,j)+sum(xf_comb3(index)>th)/K;
Pf_loc(3,3,j)=Pf_loc(3,3,j)+(sum(xf_comb3>th)-sum(xf_comb3(index)>th))/(N-K);
RMSE_loc(3,3,j)=RMSE_loc(3,3,j)+sqrt(mean((xf-xf_comb3).^2));
 end
subplot(2, 1, 2); plot(SNR(1:s), ncor(:,1:s),'LineWidth',1.2);grid on;
title('Normalized correlation');xlabel('SNR (dB)');ylabel('Normalized correlation');
subplot(2,4,1);plot(f,xf_gm(:,1));grid on;xlabel('Frequency (GHz)');ylabel('Amplitude');
title(['Geo, SNR =',num2str(SNR(s)),'dB, L=',num2str(L(j))]);
subplot(2,4,2);plot(f,xf_m(:,3));grid on;xlabel('Frequency (GHz)');ylabel('Amplitude');
title('Arith');
subplot(2,4,3);plot(f,xf_comb);grid on;xlabel('Frequency (GHz)');ylabel('Amplitude');
title('Optimized');
subplot(2,4,4);plot(f,xf_comb3);grid on;xlabel('Frequency (GHz)');ylabel('Amplitude');
title('Ovelapping');drawnow;
 end
Pd(:,:,:,s)=Pd_loc/MC; Pf(:,:,:,s)=Pf_loc/MC; RMSE(:,:,:,s)=RMSE_loc/MC;
end
leg = cell(1,l); 
for j = 1:l
    leg{j} = ['L = ',num2str(L(j))]; 
end
disp(['Time=', datestr(now, 'HH:MM:SS'), ', Duration(min) = ', num2str(toc / 60)]);
Qd_V_geo(:,:)=Pd(1,1,:,:); Qf_V_geo(:,:)=Pf(1,1,:,:); 
Qd_V_mean(:,:)= Pd(2,1,:,:);Qf_V_mean(:,:)= Pf(2,1,:,:);
Qd_w_geo(:,:)=Pd(1,3,:,:); Qf_w_geo(:,:)=Pf(1,3,:,:); 
Qd_w_mean(:,:)= Pd(2,3,:,:);Qf_w_mean(:,:)= Pf(2,3,:,:);   
Qd_all_geo(:,:)=Pd(1,2,:,:); Qf_all_geo(:,:)=Pf(1,2,:,:); 
Qd_all_mean(:,:)= Pd(2,2,:,:);Qf_all_mean(:,:)= Pf(2,2,:,:);
RMSE_V_geo(:,:)=RMSE(1,1,:,:); RMSE_all_geo(:,:)=RMSE(1,2,:,:); 
RMSE_V_mean(:,:)= RMSE(2,1,:,:);RMSE_all_mean(:,:)= RMSE(2,2,:,:);
RMSE_w_geo(:,:)=RMSE(1,3,:,:); RMSE_w_mean(:,:)= RMSE(2,3,:,:);
Qd_comb(:,:)=Pd(3,1,:,:);Qf_comb(:,:)=Pf(3,1,:,:);RMSE_com(:,:)=RMSE(3,1,:,:);
Qd_comb2(:,:)=Pd(3,2,:,:);Qf_comb2(:,:)=Pf(3,2,:,:);RMSE_com2(:,:)=RMSE(3,2,:,:);
Qd_comb3(:,:)=Pd(3,3,:,:);Qf_comb3(:,:)=Pf(3,3,:,:);RMSE_com3(:,:)=RMSE(3,3,:,:);
[L_grid,SNR_grid] = meshgrid(L,SNR);methods = {'geo', 'mean','comb'};
titles_latex = {
'$|\mathbf{\hat{X}_g(f)}|=\left(\prod_{i=1}^L w_i\cdot|\mathbf{ \hat{X}_i}|\right)^{\frac{1}{L}}$', ...  
'$|\mathbf{\bar{X}(f)}|=\sum_{i=1}^L w_i |\mathbf{\hat{X}_i}|$', ...  
'$|\mathbf{\hat{X}_{comb}(f)}|=\rho \cdot|\mathbf{\bar{X}(f)}|+(1-\rho)\cdot |\mathbf{\hat{X}_g(f)}|$'};
figure;legend(leg, 'Location', 'best');

for i = 1:2
Qd_V = eval(['Qd_V_' methods{i}]); Qf_V = eval(['Qf_V_' methods{i}]);
RMSE_V = eval(['RMSE_V_' methods{i}]); Qd_all = eval(['Qd_all_' methods{i}]);
Qf_all = eval(['Qf_all_' methods{i}]); RMSE_all = eval(['RMSE_all_' methods{i}]); 
Qd_w = eval(['Qd_w_' methods{i}]); Qf_w = eval(['Qf_w_' methods{i}]); 
RMSE_w = eval(['RMSE_w_' methods{i}]); 

subplot(3, 3, i);    
surf(L_grid, SNR_grid, Qd_all', 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeColor', 'blue', 'LineWidth', 1.3);hold on; grid on;    
surf(L_grid, SNR_grid, Qd_V', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'red', 'LineWidth', 1.3);
surf(L_grid, SNR_grid, Qd_w', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'green','LineStyle','--', 'LineWidth', 1.3);
legend('All Ch', 'Virt Ch', 'Weighted Ch', 'Location', 'best'); 

zlabel('Recovered Support'); ylabel('SNR (dB)'); 
xlabel('L'); title(titles_latex{i}, 'Interpreter', 'latex'); 
subplot(3, 3, 3 + i);
surf(L_grid, SNR_grid, Qf_all', 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeColor', 'blue', 'LineWidth', 1.3);hold on; grid on;
surf(L_grid, SNR_grid, Qf_V', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'red', 'LineWidth', 1.3);
surf(L_grid, SNR_grid, Qf_w', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'green','LineStyle','--', 'LineWidth', 1.3);
title(['Qf - ' methods{i}]);set(gca,'ZScale','log');
zlabel('Erroneous Support');    ylabel('SNR (dB)'); xlabel('L');
subplot(3, 3, 2 * 3 + i);
surf(L_grid, SNR_grid, RMSE_all', 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeColor', 'blue', 'LineWidth', 1.3);hold on; grid on;
surf(L_grid, SNR_grid, RMSE_V', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'red', 'LineWidth', 1.3);
surf(L_grid, SNR_grid, RMSE_w', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'green','LineStyle','--', 'LineWidth', 1.3);
zlabel('RMSE'); ylabel('SNR (dB)'); xlabel('L'); set(gca, 'ZScale', 'log');set(gca, 'XScale', 'log');
title(['RMSE - ' methods{i}]);
end

subplot(3, 3, 3);surf(L_grid, SNR_grid, Qd_comb', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'LineWidth', 1.3);
grid on; hold on; surf(L_grid, SNR_grid, Qd_comb2', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'cyan', 'LineWidth', 1.3);
surf(L_grid, SNR_grid, Qd_comb3', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'magenta', 'LineWidth', 1.3);
legend('Opt Comb', 'Corr Comb','Overlap Comb', 'Location', 'best'); 
title(titles_latex{3}, 'Interpreter', 'latex');set(gca,'XScale','log');
zlabel('Recovered Support'); ylabel('SNR (dB)'); xlabel('L');
subplot(3, 3, 6);surf(L_grid, SNR_grid, Qf_comb', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'LineWidth', 1.3);
grid on; hold on; surf(L_grid, SNR_grid, Qf_comb2', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'cyan', 'LineWidth', 1.3);
surf(L_grid, SNR_grid, Qf_comb3', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'magenta', 'LineWidth', 1.3);
set(gca,'ZScale','log');set(gca,'XScale','log');
zlabel('Erroneous Support');    ylabel('SNR (dB)'); xlabel('L');

subplot(3, 3, 9);surf(L_grid, SNR_grid, RMSE_com', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'black', 'LineWidth', 1.3);
grid on; hold on; surf(L_grid, SNR_grid, RMSE_com2', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'cyan', 'LineWidth', 1.3);
surf(L_grid, SNR_grid, RMSE_com3', 'FaceColor', 'none', 'FaceAlpha', 0.3, 'EdgeColor', 'magenta', 'LineWidth', 1.3);
zlabel('RMSE'); ylabel('SNR (dB)'); xlabel('L'); set(gca, 'ZScale', 'log');
set(gca,'XScale','log');
curves(SNR,L,Qd_w_geo,Qf_w_geo,RMSE_w_geo,Qd_V_geo,Qd_w_mean,Qd_comb,Qd_comb2,Qd_comb3,Qf_V_geo,Qf_w_mean,Qf_comb,Qf_comb2,Qf_comb3,RMSE_V_geo,RMSE_w_mean,RMSE_com,RMSE_com2,RMSE_com3)

tightfig;
