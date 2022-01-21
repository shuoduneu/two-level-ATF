function f=featureCalc_brachial(S,bottom_ref_x,fs)
% This function calculates features from the brachial pulse wave
% bottom_ref_x: the reference for finding the foot.

dS1=diff5(S,1,1);
dS2=diff5(S,2,1);
N=length(bottom_ref_x);

% find the bottom
bottom=nan(length(bottom_ref_x),2);
win=round(0.05*fs);
for ii=1:N
    [~,bottom_x]=min(S(bottom_ref_x(ii)-win:bottom_ref_x(ii)+win));
    bottom_x=bottom_ref_x(ii)-win+bottom_x-1;
    bottom(ii,:)=[bottom_x,S(bottom_x)];
end

% find peak
win=0.5*fs;
peak=nan(size(bottom));
for ii=1:N
        [~,peak_x]=max(S(bottom(ii,1):bottom(ii,1)+win));
    peak_x=bottom(ii,1)+peak_x-1;
    peak(ii,:)=[peak_x,S(peak_x)];
end

% find the onset
onset=nan(size(peak));
for ii=1:N
    [maxdS1_y,maxdS1_x]=max(dS1(bottom(ii,1):peak(ii,1)));
    maxdS1_x=maxdS1_x+bottom(ii,1)-1;
    onset_x=maxdS1_x-round((S(maxdS1_x)-bottom(ii,2))/maxdS1_y);
    if onset_x<bottom(ii,:)
        onset_x=bottom(ii,:);
        warning('Foot calculated using intersection tangent appears before the bottom')
    end
    onset(ii,:)=[onset_x,S(onset_x)];
end

% find the notch
notch=nan(size(peak));
for ii=1:N
    [~,notch_x]=max(dS2(peak(ii,1):peak(ii,1)+round(0.30*fs)));
    notch_x=notch_x+peak(ii,1)-1+1;
    notch(ii,:)=[notch_x,S(notch_x)];
end

% find the instability point
instable=nan(size(peak));
for ii=1:N
    [~,instable_x]=min(dS1(peak(ii,1):notch(ii,1)));
    instable_x=instable_x+peak(ii,1)-1;
    instable(ii,:)=[instable_x,S(instable_x)];
end

% figure;plot(S); hold on;plot(instable(:,1),instable(:,2),'*');hold off;

% calculate the parameters
T=bottom(2:end,1)-bottom(1:end-1,1);T=[T;nan];T=T./fs;
SBP=peak(:,2);
DBP=bottom(:,2);
MBP=nan(size(SBP));
for ii=1:N-1
    MBP(ii)=mean(S(bottom(ii,1):bottom(ii+1,1)));
    dP_min(ii)=min(diff(S(bottom(ii,1):bottom(ii+1,1))));
    dP_max(ii)=max(diff(S(bottom(ii,1):bottom(ii+1,1))));
end
PP=SBP-DBP;
FF=(MBP-DBP)./PP*100;
ED=(notch(:,1)-bottom(:,1))./fs*1000;
SV=nan(size(SBP));
for ii=1:N-1
    A=sum(S(onset(ii,1):notch(ii,1)))/fs;
    K=90./MBP(ii);
    Z=((peak(ii,2)-onset(ii,2))/(peak(ii,1)-onset(ii,1))+notch(ii,2)/(onset(ii+1,1)-notch(ii,1))-instable(ii,2)/(onset(ii+1,1)-instable(ii,1)))*fs.*K; % condition 1: no visible second peak
    SV(ii,1)=A/Z;
end
% SV=0.28*(SBP(2,:)-DBP(2,:)).*T/(FF.^2);  % ml
CO=SV.*60./T;                % ml/min
TPR=MBP./CO;
% K=(MBP-DBP)./PP;
% pct=(notch(:,2)-bottom(:,2))./PP*100;
pct=notch(:,2);
onsetp=onset(:,2);
HR=floor(60./T);
% HR=60./T;
dP_min=dP_min';
dP_max=dP_max';

% save all in structure f
f.peak=peak;
f.bottom=bottom;
f.onset=onset;
f.notch=notch;
% f.T=T(1:end-1);
% f.F=F(1:end-1);
f.SBP=SBP(1:end-1);
f.DBP=DBP(1:end-1);
f.MBP=MBP(1:end-1);
f.PP=PP(1:end-1);
% f.FF=FF(1:end-1);
f.ED=ED(1:end-1);
f.SV=SV(1:end-1);
f.CO=CO(1:end-1);
% f.K=K(1:end-1);
% f.pct=pct(1:end-1);
% f.onsetp=onsetp(1:end-1);
f.HR=HR(1:end-1);
f.dP_max=dP_max(1:end-1)/(1/fs);
f.dP_min=dP_min(1:end-1)/(1/fs);
% f.TPR=TPR(1:end-1);