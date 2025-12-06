function curves(SNR,L,Qd_geo,Qf_geo,RMSE_geo,Qd_1,Qd_2,Qd_3,Qd_4,Qd_5,Qf1,Qf_2,Qf_3,Qf_4,Qf_5,...
RMSE_1,RMSE_2,RMSE_3,RMSE_4,RMSE_5)
figure;subplot(1,3,1); colors = lines(6); 
for i = [1, length(L)] 
    if i== 1
        lineStyle = '-'; 
    else  
        lineStyle = '--';
    end    
semilogy(SNR, Qd_1(i,:),'Color',colors(1,:),'LineStyle', lineStyle, 'LineWidth', 2, 'DisplayName', ['L = ' num2str(L(i)) ', Geo Virt']); hold on; grid on;
plot(SNR,Qd_geo(i,:),'Color',colors(6,:),'LineStyle', lineStyle, 'LineWidth', 2, 'DisplayName', ['L = ' num2str(L(i)) ', Geo Real']);
plot(SNR,Qd_2(i,:),'Color',colors(2,:),'LineStyle', lineStyle, 'LineWidth', 2, 'DisplayName', ['L = ' num2str(L(i)) ', Arith Weight']);
plot(SNR,Qd_3(i,:),'Color',colors(3,:),'LineStyle', lineStyle, 'LineWidth', 2, 'DisplayName', ['L = ' num2str(L(i)) ', Opt Comb']);
plot(SNR,Qd_4(i,:),'Color',colors(4,:),'LineStyle', lineStyle, 'LineWidth', 2, 'DisplayName', ['L = ' num2str(L(i)) ', Corr Comb']);
plot(SNR,Qd_5(i,:),'Color',colors(5,:),'LineStyle', lineStyle, 'LineWidth', 2, 'DisplayName', ['L = ' num2str(L(i)) ', Overlap Comb']);
end
xlabel('SNR (dB)');ylabel('Recovered Support');
title('Q_d vs. SNR for L = 4 and L = 16');legend('Location', 'best');

%% Qf
subplot(1,3,2)
for i = [1, length(L)] 
    if i== 1
        lineStyle = '-'; 
    else  
        lineStyle = '--';
    end   
    
semilogy(SNR, Qf1(i, :), 'Color', colors(1,:),'LineStyle',lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Geo Virt']); hold on; grid on;
plot(SNR,Qf_geo(i,:),'Color',colors(6,:),'LineStyle', lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Geo Real']);
plot(SNR, Qf_2(i, :), 'Color', colors(2,:),'LineStyle',lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Arith Weight']);
plot(SNR, Qf_3(i, :), 'Color', colors(3,:),'LineStyle',lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Opt Comb']);
plot(SNR, Qf_4(i, :), 'Color', colors(4,:),'LineStyle',lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Corr Comb']);
plot(SNR, Qf_5(i, :), 'Color', colors(5,:),'LineStyle',lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Overlap Comb']);
end
xlabel('SNR (dB)');ylabel('Erroneous Support');
title('Q_f vs. SNR for L = 3 and L = 18');

%% RMSE
subplot(1,3,3)
for i = [1, length(L)] 
    if i== 1
        lineStyle = '-'; 
    else  
        lineStyle = '--';
    end    
semilogy(SNR, RMSE_1(i, :), 'Color', colors(1, :), 'LineStyle', lineStyle, 'LineWidth', 2,'DisplayName', ['L = ' num2str(L(i)) ', Geo Virt']); hold on; grid on;
plot(SNR,RMSE_geo(i,:),'Color',colors(6,:),'LineStyle', lineStyle,'LineWidth',2,'DisplayName', ['L = ' num2str(L(i)) ', Geo Real']);
plot(SNR, RMSE_2(i, :), 'Color', colors(2, :), 'LineStyle', lineStyle, 'LineWidth', 2,'DisplayName', ['L = ' num2str(L(i)) ', Arith Weight']);
plot(SNR, RMSE_3(i, :), 'Color', colors(3, :), 'LineStyle', lineStyle, 'LineWidth', 2,'DisplayName', ['L = ' num2str(L(i)) ', Opt Comb']);
plot(SNR, RMSE_4(i, :), 'Color', colors(4, :), 'LineStyle', lineStyle, 'LineWidth', 2,'DisplayName', ['L = ' num2str(L(i)) ', Corr Comb']);
plot(SNR, RMSE_5(i, :), 'Color', colors(5, :), 'LineStyle', lineStyle, 'LineWidth', 2,'DisplayName', ['L = ' num2str(L(i)) ', Overlap Comb']);
end

xlabel('SNR (dB)');ylabel('RMSE');
title('RMSE vs. SNR for L = 3 and L = 18');
