function f=featureCalc_aorta(S,varargin)
% This function calculates features from the aortic pulse wave
% Time axis of foot points calculated from the aortic pulse wave in 'align'
% step is considerd reference for the extraction of foot points in others
% steps.
% The number of input in this function is 3, when a reference is available,and 2 when not. 

dS1=diff5(S,1,1);
dS2=diff5(S,2,1);

switch nargin
    case 2
        fs=varargin{1};
        % find the peak
        thr=(max(S)-min(S))*0.7;
        count=0;
        peak_y=[];
        peak_x=[];
        win=0.5*fs;
        for ii=win+1:length(S)-win
            if abs(max(S)-S(ii))<thr&&S(ii)==max(S(ii-win:ii+win))
                if count==0
                    peak_y(1,1)=S(ii);
                    peak_x(1,1)=ii;
                    count=1;
                end
                if count~=0&&ii-peak_x(count)>round(0.3*fs)
                    peak_y(count+1,1)=S(ii);
                    peak_x(count+1,1)=ii;
                    count=count+1;
                elseif count~=0&&ii-peak_x(count)<=round(0.3*fs)
                    fault=[count;peak_x];
                end
            end
        end
        peak=[peak_x,peak_y];
        
        % find the bottom
        N=size(peak,1);
        bottom=nan(size(peak));
        for ii=1:N
            [~,bottom_x]=min(S(max(peak(ii,1)-round(0.5*fs),1):peak(ii,1)));
            bottom_x=max(peak(ii,1)-round(0.5*fs),1)+bottom_x-1;
            bottom(ii,:)=[bottom_x,S(bottom_x)];
        end
    case 3
        bottom_ref_x=varargin{1};
        fs=varargin{2};
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
    onset(ii,:)=[onset_x;S(onset_x)];
end

% figure;plot(S); hold on;plot(onset(:,1),onset(:,2),'*');hold off;

% find the inflection point
inflect=nan(size(peak));
for ii=1:N
    [~,mindS2_x]=min(dS2(onset(ii,1):onset(ii,1)+round(0.08*fs)));
    mindS2_x=mindS2_x+onset(ii,1)-1+1;
    [~,inflect_x]=max(dS2(mindS2_x:peak(ii,1)));%Inf: inflection point
    inflect_x=inflect_x+mindS2_x-1+1;
    inflect(ii,:)=[inflect_x;S(inflect_x)];
end

% find the notch
notch=nan(size(peak));
for ii=1:N
    [~,notch_x]=max(dS2(peak(ii,1):peak(ii,1)+round(0.30*fs)));
    notch_x=notch_x+peak(ii,1)-1+1;
    notch(ii,:)=[notch_x;S(notch_x)];
end

% calculate parameters
SBP=peak(:,2);
DBP=bottom(:,2);
MBP=nan(size(SBP));
for ii=1:N-1
    MBP(ii)=mean(S(bottom(ii,1):bottom(ii+1,1)));
    dP_min(ii)=min(diff(S(bottom(ii,1):bottom(ii+1,1))));
    dP_max(ii)=max(diff(S(bottom(ii,1):bottom(ii+1,1))));
end
PP=SBP-DBP;
ED=notch(:,1)-bottom(:,1)/fs*1000;
AI=(SBP-inflect(:,2))./PP*100;
% AI=inflect(:,2);
% K=(MBP-DBP)./PP;
FF=(MBP-DBP)./PP*100;
% pct=(notch(:,2)-bottom(:,2))./PP*100;
pct=(notch(:,2)-bottom(:,2));
% HR=floor(size(DBP,2))*60*fs/(DBP(end)-DBP(1));
dP_min=dP_min';
dP_max=dP_max';





% save all in structure f
f.peak=peak;
f.bottom=bottom;
f.onset=onset;
f.inflect=inflect;
f.notch=notch;
f.SBP=SBP(1:end-1);
f.DBP=DBP(1:end-1);
f.MBP=MBP(1:end-1);
f.PP=PP(1:end-1);
f.ED=ED(1:end-1);
f.AI=AI(1:end-1);
f.FF=FF(1:end-1);
f.pct=pct(1:end-1);
% f.HR=HR;
f.dP_max=dP_max(1:end-1);
f.dP_min=dP_min(1:end-1);


