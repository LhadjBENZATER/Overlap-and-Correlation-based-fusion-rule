function result = demo_fusion(makePlot)
% Demonstrate the repository's comb.m on generated nonnegative spectra.
if nargin<1,makePlot=true;end
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'examples'));
rng(20261121,'twister');N=256;L=4;threshold=0.55;
truth=false(N,1);truth([35:49 140:157])=true;
amplitude=zeros(N,1);amplitude(truth)=1;
spectra=abs(repmat(amplitude,1,L)+0.3*(randn(N,L)+1i*randn(N,L)));
arithmetic=mean(spectra,2);geometric=exp(mean(log(max(spectra,realmin)),2));
combined=comb(arithmetic,geometric,threshold);
curves=[arithmetic,geometric,combined];values=zeros(3,3);
for method=1:3,values(method,:)=support_metrics(truth,curves(:,method)>threshold);end
result=struct('metrics',values,'spectra',curves,'threshold',threshold,'seed',20261121);
out=fullfile(root,'results');if ~exist(out,'dir'),mkdir(out);end
write_demo_csv(fullfile(out,'demo_metrics.csv'),'method,pd,pf,f1',[(1:3)',values]);
save(fullfile(out,'demo_results.mat'),'result');
fprintf('Fusion demo: fixed amplitude threshold %.2f; no threshold calibration.\n',threshold);
disp(array2table(values,'VariableNames',{'Pd','Pf','F1'},'RowNames',{'Arithmetic','Geometric','Overlap_comb'}));
if makePlot
 fig=figure('Visible','off','Color','w','Position',[100 100 1100 600]);
 subplot(2,1,1);plot(curves,'LineWidth',1.3);hold on;plot([1 N],[threshold threshold],'k--');grid on;
 xlabel('Bin');ylabel('Amplitude');legend('Arithmetic mean','Geometric mean','Overlap comb','Threshold','Location','best');
 subplot(2,1,2);bar(values);ylim([0 1]);grid on;set(gca,'XTickLabel',{'Arithmetic','Geometric','Overlap comb'});legend('P_d','P_f','F1','Location','best');
 sgtitle('Synthetic amplitude-fusion example using the existing comb.m');
 print(fig,fullfile(root,'docs','preview.svg'),'-dsvg');close(fig);
end
end
