clc;close all;

load('F:\PW基本处理\两级传递函数\result\幅值+相角\TF.mat');
load('F:\PW基本处理\两级传递函数\result\幅值+相角\sig.mat');
load('F:\PW基本处理\两级传递函数\result\幅值+相角\feature.mat');
load('F:\PW基本处理\两级传递函数\result\幅值+相角\para.mat');
load('F:\PW基本处理\两级传递函数\result\幅值+相角\result.mat');

figure;
subplot(3,1,1)
plot(sig.align{1, 2}.aorta,'linewidth',0.8);hold on;plot(sig.normalize{1, 2}.aorta_base,'k','linewidth',0.8);hold off; ylabel({'Pressure';'(mmHg)'});set(gca,'Fontsize',8,'FontName','Times New Roman','ytick',[]) 
axis off;
subplot(3,1,2)
plot(sig.align{1, 2}.aorta-sig.normalize{1, 2}.aorta_base,'linewidth',0.8);ylabel({'Pressure' ;'(mmHg)'});set(gca,'Fontsize',8,'FontName','Times New Roman','ytick',[]) 
axis off;
subplot(3,1,3)
plot(sig.normalize{1, 2}.aorta,'linewidth',0.8);ylabel('Amplitude');xlabel({'Sampling points'});set(gca,'Fontsize',8,'FontName','Times New Roman','ytick',[]) 
axis off;


for si=30
    figure;
    m=4;
%     f(1)=plot(sig_MP.sig.align{1, si}.aorta(feature_MP.feature.align{1, si}.aorta.bottom(m,1):feature.align{1, si}.aorta.bottom(m+1,1)),'k*');hold on;

    f(1)=plot(sig.align{1, si}.aorta(feature.align{1, si}.aorta.bottom(m,1):feature.align{1, si}.aorta.bottom(m+1,1)),'k','linewidth',2);hold on;
%     ylabel('Aortic Pressure (mmHg)');xlabel('Sampling points');
%     axis off;
    f(2)=plot(sig.est.general_cv.align{1, si}(feature.interp.est.general_cv.align{1, si}.bottom(m,1):feature.interp.est.general_cv.align{1, si}.bottom(m+1,1)),'b:','linewidth',0.8);hold on;
%     axis off;
    f(3)=plot(sig.est.general_cv.calibrate{1, si}(feature.interp.est.general_cv.calibrate{1, si}.bottom(m,1):feature.interp.est.general_cv.calibrate{1, si}.bottom(m+1,1)),'g-.','linewidth',0.8);hold on;
%     axis off;
    f(4)=plot(sig.est.adaptive_cv.align{1, si}(feature.interp.est.adaptive_cv.align{1, si}.bottom(m,1):feature.interp.est.adaptive_cv.align{1, si}.bottom(m+1,1)),'m--','linewidth',0.8);hold on;
%     axis off;
    f(5)=plot(sig.est.adaptive_cv.calibrate{1, si}(feature.interp.est.adaptive_cv.calibrate{1, si}.bottom(m,1):feature.interp.est.adaptive_cv.calibrate{1, si}.bottom(m+1,1)),'r','linewidth',0.8);hold on;
%     axis off;
    
    legend([f(1),f(2),f(3),f(4),f(5)],'Measured','GTF','Two-level GTF strategy','ATF','Two-level ATF strategy');%%长8宽4
    legend('boxoff') 
    set(gca,'Fontsize',8,'FontName','Times New Roman')
    xlabel({'Sampling points'});ylabel('Estimated APWs (mmHg)')
end
% sig_multiple.sig.est.adaptive_cv.align{1, si}(feature_multiple.feature.interp.est.adaptive_cv.calibrate{1, si}.bottom(1,1):feature_multiple.feature.interp.est.adaptive_cv.calibrate{1, si}.bottom(2,1))-...
% sig.est.adaptive_cv.align{1, si}(feature.interp.est.adaptive_cv.calibrate{1, si}.bottom(1,1):feature.interp.est.adaptive_cv.calibrate{1, si}.bottom(2,1))    

% for si=1:34
%     figure;
%     m=4;
% %     f(1)=plot(sig_MP.sig.align{1, si}.aorta(feature_MP.feature.align{1, si}.aorta.bottom(m,1):feature.align{1, si}.aorta.bottom(m+1,1)),'k*');hold on;
% 
%     f(1)=plot(sig.align{1, si}.aorta(50:300),'k','linewidth',2);hold on;
% %     ylabel('Aortic Pressure (mmHg)');xlabel('Sampling points');
% %     axis off;
%     f(2)=plot(sig.est.general_cv.align{1, si}(50:300),'b:','linewidth',0.8);hold on;
% %     axis off;
%     f(3)=plot(sig.est.general_cv.calibrate{1, si}(50:300),'g-.','linewidth',0.8);hold on;
% %     axis off;
%     f(4)=plot(sig.est.adaptive_cv.align{1, si}(50:300),'m--','linewidth',0.8);hold on;
% %     axis off;
%     f(5)=plot(sig.est.adaptive_cv.calibrate{1, si}(50:300),'r','linewidth',0.8);hold on;
% %     axis off;
%     
%     legend([f(1),f(2),f(3),f(4),f(5)],'Measured','GTF','Two-level GTF strategy','ATF','Two-level ATF strategy');%%长8宽4
%     legend('boxoff') 
%     set(gca,'Fontsize',8,'FontName','Times New Roman')
%     xlabel({'Sampling points';['Subject ',num2str(si)]});ylabel('Estimated APWs (mmHg)')
% end
% 

True=[]; GTF1=[]; ATF1=[];GTF2=[]; ATF2=[];
for si=1:34
    True=[True;sig.align{1, si}.aorta];
    GTF1=[GTF1;sig.est.general_cv.align{1, si}];
    GTF2=[GTF2;sig.est.general_cv.calibrate{1, si}];
    ATF1=[ATF1;sig.est.adaptive_cv.align{1, si}];
    ATF2=[ATF2;sig.est.adaptive_cv.calibrate{1, si}];

end
[r1,p1]=corrcoef(True,GTF1);
r1(2,1).^2
[r2,p2]=corrcoef(True,GTF2);
r2(2,1).^2
[r3,p3]=corrcoef(True,ATF1);
r3(2,1).^2
[r4,p4]=corrcoef(True,ATF2);
r4(2,1).^2



figure;
subplot(2,2,1)
plot(True,GTF1,'b+');hold on;plot(True,True,'r');hold off;
gtext(['R^2=',num2str(round(r1(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'Measured aortic pressure (mmHg)'})
ylabel({'Estimated aortic pressure using';' GTF (mmHg)'} )
xlim([min(True),max(True)]);
ylim([min(True),max(True)]);
set(gca,'Fontsize',8,'FontName','Times New Roman')
subplot(2,2,2)
plot(True,GTF2,'b+');hold on;plot(True,True,'r');hold off;
gtext(['R^2=',num2str(round(r2(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'Measured aortic pressure (mmHg)'})
ylabel({'Estimated aortic pressure using';'two-leavel GTF stragy (mmHg)'} )
xlim([min(True),max(True)]);
ylim([min(True),max(True)]);
set(gca,'Fontsize',8,'FontName','Times New Roman')
subplot(2,2,3)
plot(True,ATF1,'b+');hold on;plot(True,True,'r');hold off;
gtext(['R^2=',num2str(round(r3(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'Measured aortic pressure (mmHg)'})
ylabel({'Estimated aortic pressure using';'ATF (mmHg)'} )
xlim([min(True),max(True)]);
ylim([min(True),max(True)]);
set(gca,'Fontsize',8,'FontName','Times New Roman')
subplot(2,2,4)
plot(True,ATF2,'b+');hold on;plot(True,True,'r');hold off;
gtext(['R^2=',num2str(round(r4(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'Measured aortic pressure (mmHg)'})
ylabel({'Estimated aortic pressure using';'two-leavel ATF stragy (mmHg)'} )
xlim([min(True),max(True)]);
ylim([min(True),max(True)]);
set(gca,'Fontsize',8,'FontName','Times New Roman')





% TFw=[GTF';ATF'];
% altman(TFw,True',{'Average of the estimated';'and measured aortic pressure (mmHg)'},{'Estimated error of'; 'aortic pressure (mmHg)'})




figure,%%output of estimated peak_f
[r1,p1]=corrcoef(TF.individual.align.peak_f(:,1),TF.adaptive_cv.align.peak_f(:,1));
subplot(2,2,1)
% plot(TF_single.TF.individual.align.peak_f(:,1),TF_single.TF.adaptive_cv.align.peak_f(:,1),'r+');hold on;
plot(TF.individual.align.peak_f(:,1),TF.adaptive_cv.align.peak_f(:,1),'b+');hold on;
% plot(TF.general_cv.align.peak_f(:,1),'*');hold on;
plot(TF.individual.align.peak_f(:,1),TF.individual.align.peak_f(:,1),'r');hold off;
gtext(['R^2=',num2str(round(r1(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
% legend('MATF1','SATF1')
xlabel({'True F_g_m_1 (Hz)'})
ylabel('Estimated F_g_m_1 (Hz)')
xlim([min(TF.individual.align.peak_f(:,1)),max(TF.individual.align.peak_f(:,1))]);
ylim([min(TF.individual.align.peak_f(:,1)),max(TF.individual.align.peak_f(:,1))]);
set(gca,'Fontsize',8,'FontName','Times New Roman')

subplot(2,2,2)
[r1,p1]=corrcoef(TF.individual.align.peak_f(:,2),TF.adaptive_cv.align.peak_f(:,2));
% plot(TF_single.TF.individual.align.peak_f(:,1),TF_single.TF.adaptive_cv.align.peak_f(:,1),'r+');hold on;
plot(TF.individual.align.peak_f(:,2),TF.adaptive_cv.align.peak_f(:,2),'b+');hold on;
% plot(TF.general_cv.align.peak_f(:,1),'*');hold on;
plot(TF.individual.align.peak_f(:,2),TF.individual.align.peak_f(:,2),'r');hold off;
% legend('MATF1','SATF1')
gtext(['R^2=',num2str(round(r1(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'True F_p_m_1 (Hz)'})
ylabel('Estimated F_p_m_1 (Hz)')
xlim([min(TF.individual.align.peak_f(:,2)),max(TF.individual.align.peak_f(:,2))]);
ylim([min(TF.individual.align.peak_f(:,2)),max(TF.individual.align.peak_f(:,2))]);
set(gca,'Fontsize',8,'FontName','Times New Roman')

subplot(2,2,3)
[r1,p1]=corrcoef(TF.individual.normalize.peak_f(:,1),TF.adaptive_cv.normalize.peak_f(:,1));
% plot(TF_single.TF.individual.align.peak_f(:,1),TF_single.TF.adaptive_cv.align.peak_f(:,1),'r+');hold on;
plot(TF.individual.normalize.peak_f(:,1),TF.adaptive_cv.normalize.peak_f(:,1),'b+');hold on;
% plot(TF.general_cv.align.peak_f(:,1),'*');hold on;
plot(TF.individual.normalize.peak_f(:,1),TF.individual.normalize.peak_f(:,1),'r');hold off;
% legend('MATF1','SATF1')
gtext(['R^2=',num2str(round(r1(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'True F_g_m_2 (Hz)'})
ylabel('Estimated F_g_m_2 (Hz)')
xlim([min(TF.individual.normalize.peak_f(:,1)),max(TF.individual.normalize.peak_f(:,1))]);
ylim([min(TF.individual.normalize.peak_f(:,1)),max(TF.individual.normalize.peak_f(:,1))]);
set(gca,'Fontsize',8,'FontName','Times New Roman')

subplot(2,2,4)
[r1,p1]=corrcoef(TF.individual.normalize.peak_f(:,2),TF.adaptive_cv.normalize.peak_f(:,2));
% plot(TF_single.TF.individual.align.peak_f(:,1),TF_single.TF.adaptive_cv.align.peak_f(:,1),'r+');hold on;
plot(TF.individual.normalize.peak_f(:,2),TF.adaptive_cv.normalize.peak_f(:,2),'b+');hold on;
% plot(TF.general_cv.align.peak_f(:,1),'*');hold on;
plot(TF.individual.normalize.peak_f(:,2),TF.individual.normalize.peak_f(:,2),'r');hold off;
% legend('MATF1','SATF1')
gtext(['R^2=',num2str(round(r1(2,1).^2,2))],'Fontsize',8,'FontName','Times New Roman');
gtext(['\itP','\rm<0.05'],'Fontsize',8,'FontName','Times New Roman');
xlabel({'True F_p_m_2 (Hz)'})
ylabel('Estimated F_p_m_2 (Hz)')
xlim([min(TF.individual.normalize.peak_f(:,2)),max(TF.individual.normalize.peak_f(:,2))]);
ylim([min(TF.individual.normalize.peak_f(:,2)),max(TF.individual.normalize.peak_f(:,2))]);
set(gca,'Fontsize',8,'FontName','Times New Roman')





figure;
subplot(2,1,1)
[mu,sigma,muci,sigmaci]=normfit(abs(TF.individual.align.H));
upper_align=muci(1,:); 
lower_align=muci(2,:); 
f(1)=plot(TF.general.align.f(1,:),mu,'b','linewidth',2);hold on;
% upper_align=abs(TF.general.align.H(1,:))+1.96*std(abs(TF.individual.align.H));
% lower_align=abs(TF.general.align.H(1,:))-1.96*std(abs(TF.individual.align.H));
f(2)=plot(TF.general.align.f(1,:),upper_align,'b:');hold on;
f(3)=plot(TF.general.align.f(1,:),lower_align,'b:');hold on;
[mu,sigma,muci,sigmaci]=normfit(abs(TF.individual.normalize.H));
upper_normalize=muci(1,:); 
lower_normalize=muci(2,:); 
f(4)=plot(TF.general.normalize.f(1,:),mu,'r','linewidth',2);hold on;
% upper_normalize=abs(TF.general.normalize.H(1,:))+1.96*std(abs(TF_single.TF.general.normalize.H));
% lower_normalize=abs(TF.general.normalize.H(1,:))-1.96*std(abs(TF_single.TF.general.normalize.H));
f(5)=plot(TF.general.normalize.f(1,:),upper_normalize,'r:');hold on;
f(6)=plot(TF.general.normalize.f(1,:),lower_normalize,'r:');hold on;
% f(7)=plot(TF.general.normalize.f(1,:),ones(length(TF.general.normalize.f(1,:)),1),'k--');hold off;

legend([f(1),f(4)],'GTF1','GTF2');
ylabel('Magnitude')
xlabel('Frequency (Hz)')
xlim([0 15])
% ylim([0 2.5])
set(gca,'Fontsize',8,'FontName','Times New Roman')  %%长7宽7

for subjIdx=1:34
    Phase1(subjIdx,:)=phase(TF.individual.align.H(subjIdx,:));
    Phase2(subjIdx,:)=phase(TF.individual.normalize.H(subjIdx,:));
end
subplot(2,1,2)
[mu,sigma,muci,sigmaci]=normfit(Phase1);
upper_align=muci(1,:); 
lower_align=muci(2,:); 
f(1)=plot(TF.general.align.f(1,:),mu,'b','linewidth',2);hold on
f(2)=plot(TF.general.align.f(1,:),upper_align,'b:');hold on;
f(3)=plot(TF.general.align.f(1,:),lower_align,'b:');hold on;
[mu,sigma,muci,sigmaci]=normfit(Phase2);
upper_normalize=muci(1,:); 
lower_normalize=muci(2,:); 
f(4)=plot(TF.general.normalize.f(1,:),mu,'r','linewidth',2);hold on;
f(5)=plot(TF.general.normalize.f(1,:),upper_normalize,'r:');hold on;
f(6)=plot(TF.general.normalize.f(1,:),lower_normalize,'r:');hold on;
% f(7)=plot(TF.general.normalize.f(1,:),zeros(length(TF.general.normalize.f(1,:)),1),'k--');
hold off;
legend([f(1),f(4)],'GTF1','GTF2');
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')
xlim([0 15])
% ylim([0 2.5])
set(gca,'Fontsize',8,'FontName','Times New Roman')


% figure;
% boxplot([result.linear.subjAverage.general_cv.align.RMSE,result.linear.subjAverage.general_cv.calibrate.RMSE,...
%     result.linear.subjAverage.adaptive_cv.align.RMSE,result.linear.subjAverage.adaptive_cv.calibrate.RMSE])
% figure;
% boxplot([result.linear.subjAverage.general_cv.align.SBP,result.linear.subjAverage.general_cv.calibrate.SBP,...
%     result.linear.subjAverage.adaptive_cv.align.SBP,result.linear.subjAverage.adaptive_cv.calibrate.SBP])
% % figure;
% % boxplot([result.linear.subjAverage.general_cv.align.DBP,result.linear.subjAverage.general_cv.calibrate.DBP,...
% %     result.linear.subjAverage.adaptive_cv.align.DBP,result.linear.subjAverage.adaptive_cv.calibrate.DBP])
% figure;
% boxplot([result.linear.subjAverage.general_cv.align.PP,result.linear.subjAverage.general_cv.calibrate.PP,...
%     result.linear.subjAverage.adaptive_cv.align.PP,result.linear.subjAverage.adaptive_cv.calibrate.PP])
% figure;
% boxplot([result.linear.subjAverage.general_cv.align.AI,result.linear.subjAverage.general_cv.calibrate.AI,...
%     result.linear.subjAverage.adaptive_cv.align.AI,result.linear.subjAverage.adaptive_cv.calibrate.AI])

average=[mean(result.linear.subjAverage.general_cv.align.RMSE),mean(result.linear.subjAverage.general_cv.calibrate.RMSE),...
    mean(result.linear.subjAverage.adaptive_cv.align.RMSE),mean(result.linear.subjAverage.adaptive_cv.calibrate.RMSE)];
variance=[std(result.linear.subjAverage.general_cv.align.RMSE),std(result.linear.subjAverage.general_cv.calibrate.RMSE),...
    std(result.linear.subjAverage.adaptive_cv.align.RMSE),std(result.linear.subjAverage.adaptive_cv.calibrate.RMSE)];
figure;
subplot(1,5,1)
b=bar(diag(average),'stack');
% xlabel('(a)','Fontsize',8,'FontName','Times New Roman')
ylabel('Errors of TW (mmHg)','Fontsize',8,'FontName','Times New Roman')
set(b(1),'FaceColor','y');set(b(2),'FaceColor','b');set(b(3),'FaceColor','g');set(b(4),'FaceColor','r');
set(gca,'XTickLabel',{'GTF ','Two-level GTF strategy','ATF ','Two-level ATF strategy'},'Fontsize',8,'FontName','Times New Roman');
% set(gca,'XTickLabel',{' ',' ',' ',' '},'Fontsize',8,'FontName','Times New Roman');
% gtext('GTF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'GTF stragy'},'Fontsize',8,'FontName','Times New Roman');
% gtext('ATF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'ATF stragy'},'Fontsize',8,'FontName','Times New Roman');
set(gca,'XTickLabelRotation',46,'Fontsize',8,'FontName','Times New Roman');
hold on;
e=errorbar(average,variance,'Linestyle','None');
set(e,'Color','k')

average=[mean(result.linear.subjAverage.general_cv.align.SBP),mean(result.linear.subjAverage.general_cv.calibrate.SBP),...
    mean(result.linear.subjAverage.adaptive_cv.align.SBP),mean(result.linear.subjAverage.adaptive_cv.calibrate.SBP)];
variance=[std(result.linear.subjAverage.general_cv.align.SBP),std(result.linear.subjAverage.general_cv.calibrate.SBP),...
    std(result.linear.subjAverage.adaptive_cv.align.SBP),std(result.linear.subjAverage.adaptive_cv.calibrate.SBP)];
subplot(1,5,2)
b=bar(diag(average),'stack');
% xlabel('(b)','Fontsize',8,'FontName','Times New Roman')
ylabel('Errors of SBP (mmHg)','Fontsize',8,'FontName','Times New Roman')
set(b(1),'FaceColor','y');set(b(2),'FaceColor','b');set(b(3),'FaceColor','g');set(b(4),'FaceColor','r');
set(gca,'XTickLabel',{'GTF ','Two-level GTF strategy','ATF ','Two-level ATF strategy'},'Fontsize',8,'FontName','Times New Roman');
% set(gca,'XTickLabel',{' ',' ',' ',' '},'Fontsize',8,'FontName','Times New Roman');
% gtext('GTF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'GTF stragy'},'Fontsize',8,'FontName','Times New Roman');
% gtext('ATF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'ATF stragy'},'Fontsize',8,'FontName','Times New Roman');
set(gca,'XTickLabelRotation',46,'Fontsize',8,'FontName','Times New Roman');
hold on;
e=errorbar(average,variance,'Linestyle','None');
set(e,'Color','k')

average=[mean(result.linear.subjAverage.general_cv.align.DBP),mean(result.linear.subjAverage.general_cv.calibrate.DBP),...
    mean(result.linear.subjAverage.adaptive_cv.align.DBP),mean(result.linear.subjAverage.adaptive_cv.calibrate.DBP)];
variance=[std(result.linear.subjAverage.general_cv.align.DBP),std(result.linear.subjAverage.general_cv.calibrate.DBP),...
    std(result.linear.subjAverage.adaptive_cv.align.DBP),std(result.linear.subjAverage.adaptive_cv.calibrate.DBP)];
subplot(1,5,3)
b=bar(diag(average),'stack');
% xlabel('(c)','Fontsize',8,'FontName','Times New Roman')
ylabel('Errors of DBP (mmHg)','Fontsize',8,'FontName','Times New Roman')
set(b(1),'FaceColor','y');set(b(2),'FaceColor','b');set(b(3),'FaceColor','g');set(b(4),'FaceColor','r');
set(gca,'XTickLabel',{'GTF ','Two-level GTF strategy','ATF ','Two-level ATF strategy'},'Fontsize',8,'FontName','Times New Roman');
% set(gca,'XTickLabel',{' ',' ',' ',' '},'Fontsize',8,'FontName','Times New Roman');
% gtext('GTF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'GTF stragy'},'Fontsize',8,'FontName','Times New Roman');
% gtext('ATF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'ATF stragy'},'Fontsize',8,'FontName','Times New Roman');
set(gca,'XTickLabelRotation',46,'Fontsize',8,'FontName','Times New Roman');
hold on;
e=errorbar(average,variance,'Linestyle','None');
set(e,'Color','k')


average=[mean(result.linear.subjAverage.general_cv.align.PP),mean(result.linear.subjAverage.general_cv.calibrate.PP),...
    mean(result.linear.subjAverage.adaptive_cv.align.PP),mean(result.linear.subjAverage.adaptive_cv.calibrate.PP)];
variance=[std(result.linear.subjAverage.general_cv.align.PP),std(result.linear.subjAverage.general_cv.calibrate.PP),...
    std(result.linear.subjAverage.adaptive_cv.align.PP),std(result.linear.subjAverage.adaptive_cv.calibrate.PP)];
subplot(1,5,4)
b=bar(diag(average),'stack');
% xlabel('(d)','Fontsize',8,'FontName','Times New Roman')
ylabel('Errors of PP (mmHg)','Fontsize',8,'FontName','Times New Roman')
set(b(1),'FaceColor','y');set(b(2),'FaceColor','b');set(b(3),'FaceColor','g');set(b(4),'FaceColor','r');
set(gca,'XTickLabel',{'GTF ','Two-level GTF strategy','ATF ','Two-level ATF strategy'},'Fontsize',8,'FontName','Times New Roman');
% set(gca,'XTickLabel',{' ',' ',' ',' '},'Fontsize',8,'FontName','Times New Roman');
% gtext('GTF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'GTF stragy'},'Fontsize',8,'FontName','Times New Roman');
% gtext('ATF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'ATF stragy'},'Fontsize',8,'FontName','Times New Roman');
set(gca,'XTickLabelRotation',46,'Fontsize',8,'FontName','Times New Roman');
hold on;
e=errorbar(average,variance,'Linestyle','None');
set(e,'Color','k')

average=[mean(result.linear.subjAverage.general_cv.align.AI),mean(result.linear.subjAverage.general_cv.calibrate.AI),...
    mean(result.linear.subjAverage.adaptive_cv.align.AI),mean(result.linear.subjAverage.adaptive_cv.calibrate.AI)];
variance=[std(result.linear.subjAverage.general_cv.align.AI),std(result.linear.subjAverage.general_cv.calibrate.AI),...
    std(result.linear.subjAverage.adaptive_cv.align.AI),std(result.linear.subjAverage.adaptive_cv.calibrate.AI)];
subplot(1,5,5)
b=bar(diag(average),'stack');
% xlabel('(e)','Fontsize',8,'FontName','Times New Roman')
ylabel('Errors of AIx (%)','Fontsize',8,'FontName','Times New Roman')
set(b(1),'FaceColor','y');set(b(2),'FaceColor','b');set(b(3),'FaceColor','g');set(b(4),'FaceColor','r');
set(gca,'XTickLabel',{'GTF ','Two-level GTF strategy','ATF ','Two-level ATF strategy'},'Fontsize',8,'FontName','Times New Roman');
% set(gca,'XTickLabel',{' ',' ',' ',' '},'Fontsize',8,'FontName','Times New Roman');
% gtext('GTF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'GTF stragy'},'Fontsize',8,'FontName','Times New Roman');
% gtext('ATF','Fontsize',8,'FontName','Times New Roman');
% gtext({'Two-level';'ATF stragy'},'Fontsize',8,'FontName','Times New Roman');
set(gca,'XTickLabelRotation',46,'Fontsize',8,'FontName','Times New Roman');
hold on;
e=errorbar(average,variance,'Linestyle','None');
set(e,'Color','k')





% for subjIdx=1:para.Nsubj
% SBP_est_g1(subjIdx)=mean(feature.interp.est.general_cv.align{1, subjIdx}.SBP);
% DBP_est_g1(subjIdx)=mean(feature.interp.est.general_cv.align{1, subjIdx}.DBP);
% PP_est_g1(subjIdx)=mean(feature.interp.est.general_cv.align{1, subjIdx}.PP);
% MBP_est_g1(subjIdx)=mean(feature.interp.est.general_cv.align{1, subjIdx}.MBP);
% ED_est_g1(subjIdx)=mean(feature.interp.est.general_cv.align{1, subjIdx}.ED);
% AI_est_g1(subjIdx)=mean(feature.interp.est.general_cv.align{1, subjIdx}.AI);
% 
% SBP_est_g2(subjIdx)=mean(feature.interp.est.general_cv.calibrate{1, subjIdx}.SBP);
% DBP_est_g2(subjIdx)=mean(feature.interp.est.general_cv.calibrate{1, subjIdx}.DBP);
% PP_est_g2(subjIdx)=mean(feature.interp.est.general_cv.calibrate{1, subjIdx}.PP);
% MBP_est_g2(subjIdx)=mean(feature.interp.est.general_cv.calibrate{1, subjIdx}.MBP);
% ED_est_g2(subjIdx)=mean(feature.interp.est.general_cv.calibrate{1, subjIdx}.ED);
% AI_est_g2(subjIdx)=mean(feature.interp.est.general_cv.calibrate{1, subjIdx}.AI);
% 
% 
% SBP_est_MP1(subjIdx)=mean(feature.interp.est.adaptive_cv.align{1, subjIdx}.SBP);
% DBP_est_MP1(subjIdx)=mean(feature.interp.est.adaptive_cv.align{1, subjIdx}.DBP);
% PP_est_MP1(subjIdx)=mean(feature.interp.est.adaptive_cv.align{1, subjIdx}.PP);
% MBP_est_MP1(subjIdx)=mean(feature.interp.est.adaptive_cv.align{1, subjIdx}.MBP);
% ED_est_MP1(subjIdx)=mean(feature.interp.est.adaptive_cv.align{1, subjIdx}.ED);
% AI_est_MP1(subjIdx)=mean(feature.interp.est.adaptive_cv.align{1, subjIdx}.AI);
% 
% 
% SBP_est_MP2(subjIdx)=mean(feature.interp.est.adaptive_cv.calibrate{1, subjIdx}.SBP);
% DBP_est_MP2(subjIdx)=mean(feature.interp.est.adaptive_cv.calibrate{1, subjIdx}.DBP);
% PP_est_MP2(subjIdx)=mean(feature.interp.est.adaptive_cv.calibrate{1, subjIdx}.PP);
% MBP_est_MP2(subjIdx)=mean(feature.interp.est.adaptive_cv.calibrate{1, subjIdx}.MBP);
% ED_est_MP2(subjIdx)=mean(feature.interp.est.adaptive_cv.calibrate{1, subjIdx}.ED);
% AI_est_MP2(subjIdx)=mean(feature.interp.est.adaptive_cv.calibrate{1, subjIdx}.AI);
% 
% end
% SBP_mea=feature.subjAverage.align.aorta.SBP';
% DBP_mea=feature.subjAverage.align.aorta.DBP';
% PP_mea= feature.subjAverage.align.aorta.PP';
% MBP_mea=feature.subjAverage.align.aorta.MBP';
% ED_mea=feature.subjAverage.align.aorta.ED';
% AI_mea=feature.subjAverage.align.aorta.AI';
% 
% 
% SBP=[SBP_est_g1;SBP_est_g2;SBP_est_MP1;SBP_est_MP2];
% DBP=[DBP_est_g1;DBP_est_g2;DBP_est_MP1;DBP_est_MP2];
% PP=[PP_est_g1;PP_est_g2;PP_est_MP1;PP_est_MP2];
% MBP=[MBP_est_g1;MBP_est_g2;MBP_est_MP1;MBP_est_MP2];
% AI=[AI_est_g1;AI_est_g2;AI_est_MP1;AI_est_MP2];
% 
% 
% 
% altman(SBP,SBP_mea,{'Average of the measured';'and estimated SBP (mmHg)'},'Estimated error of SBP (mmHg)')
% altman(DBP,DBP_mea,{'Average of the measured';'and estimated DBP (mmHg)'},'Estimated error of DBP (mmHg)')
% % altman(MBP,MBP_mea,{'Average of the estimated';'and measured MBP (mmHg)'},'Estimated error of MBP (mmHg)')
% altman(PP,PP_mea,{'Average of the measured';'and estimated PP (mmHg)'},'Estimated error of PP (mmHg)')
% altman(AI,AI_mea,{'Average of the measured';'and estimated AIx (%)'},'Estimated error of AIx (%)')


