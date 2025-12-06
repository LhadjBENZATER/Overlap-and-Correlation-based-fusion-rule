clear; clc; close all; set(0, 'DefaultAxesFontSize', 14); %rng(0)
N = 1024; M = round(0.25*N); a=0.7; bk=0.65; Psi = fft(eye(N)) / sqrt(N);     
Phi = Phi_selection(N,M,2,2e7); A=Phi/Psi; At=A';th=0.1;
[x,xf,index,K,f,t,ind_shft]=LFM(2.9e9,4e7,N,1);
MC=1; l=8;Y=Phi*x; tic;
Xf_mat = zeros(N,2*l+1); 
[x1,x2,x3,x4,x5,x6,x7,x8]=Read_LFM(20);
x=[x1 x2 x3 x4 x5 x6 x7 x8];x = awgn(x,-2, 'measured');
x=fct_syn(x); x=fct_Phase(x);
y=Phi*x; YL=Virtual(y(:,1:l));          
for i = 1:2 * l + 1
    Xf_mat(:,i) = fftshift(abs(fct_MRMP(A,At,YL(:,i),K,bk,a)));
end
%% Geometric Mean
xf_gm=prod(Xf_mat(:,l+1:2*l),2).^(1/l);
mx = max(xf_gm);mx(mx==0)=1;xf_gm=abs(xf_gm)./mx;
%% Arithmetic Mean
xf_m=mean(Xf_mat,2);xf_m=xf_m-median(xf_m); xf_m=xf_m./max(xf_m);
%% Combined by Corr th=0.1
ncor=sum(xf_m.*xf_gm)/(norm(xf_m)*norm(xf_gm));
xf_comb = (1-ncor)*xf_gm+ncor*xf_m;
xf_cor=xf_comb/max(xf_comb);
%% Combined overlapping
xf_comb = comb(xf_m,xf_gm,th);
%% Joint sparse recovery
Xf_MRMP = prod(fftshift(abs(fct_MRMP(A, At, YL, K, 0.65, 0.7))), 2).^(1/8);
%% Original Spectrum
xf1 = fftshift(abs(fft(x(:,1))))/max(abs(fft(x(:,1))));

% 1. Original Spectrum
figure;plot(f, xf1, 'm', 'LineWidth', 2); grid on;
title('Original Spectrum');
xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');

% 2. Arithmetic Mean
figure;plot(f, xf_m, 'r', 'LineWidth', 1.5); grid on;
title('Arithmetic Mean');
xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');

% 3. Geometric Mean
figure;plot(f, xf_gm, 'b', 'LineWidth', 1.5); grid on;
title('Geometric Mean');
xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');

% 4. Correlation-based
figure;plot(f, xf_cor, 'Color',[1 0.5 0], 'LineWidth', 1.5); grid on;
title('Correlation-based');
xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');

% 5. Overlapping Combined
figure;plot(f, xf_comb, 'g', 'LineWidth', 1.5); grid on;
title('Overlapping Combined');
xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');

% 6. Joint Sparse Recovery (MRMP)
figure;plot(f, Xf_MRMP/max(Xf_MRMP), 'k', 'LineWidth', 1.5); grid on;
title('Joint Sparse Recovery (MRMP)');
xlabel('Frequency (GHz)'); ylabel('Normalized Spectrum');
