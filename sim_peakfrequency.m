function y=sim_peakfrequency(S,A,B,fp,fp0,fs)
% This function calculates the adaptive TF and estimates the aortic pulse
% waves.
% S: input signal
% A,B: A and B parameters of GTF
% fp: the peak frequency of the GTF
% fp0: the target peak frequency

n=1024; % number of sample points frequency response
fp_x=round(fp/fs*n);
fp0_x=round(fp0/fs*n);


X=freqz(B,A,n,'whole');
H=abs(X); % magnitude
theta=phase(X); % phase

% temp=interp(H(1:fp_x(1)),fp0_x(1));
% Hhead=temp(1:fp_x(1):end);
% temp=interp(H(fp_x(1)+1:n/2+1),n/2+1-fp0_x(1));
% Hend=temp(1:(n/2+1-fp_x(1)):end);
% H_ATF=[Hhead;Hend];
% H_ATF=[H_ATF;H_ATF(end-1:-1:2)];

Hhead_xx=linspace(1,fp_x(1),fp0_x(1));
Hhead_yy=interp1(1:fp_x(1),H(1:fp_x(1)),Hhead_xx,'linear');
Hend_xx=linspace(fp_x(1)+1,n/2+1,n/2+1-fp0_x(1)); 
Hend_yy=interp1(fp_x(1)+1:n/2+1,H(fp_x(1)+1:n/2+1),Hend_xx,'linear');
H_ATF=[Hhead_yy,Hend_yy];
H_ATF=[H_ATF,H_ATF(end-1:-1:2)];
H_ATF=H_ATF';


if length(fp0_x)==2
% temp=interp(theta(1:fp_x(2)),fp0_x(2));
% Hhead=temp(1:fp_x(2):end);
% temp=interp(theta(fp_x(2)+1:n/2+1),n/2+1-fp0_x(2));
% Hend=temp(1:(n/2+1-fp_x(2)):end);
% theta_ATF=[Hhead;Hend];
% theta_ATF=[theta_ATF;-1*theta_ATF(end-1:-1:2)];

Hhead_xx=linspace(1,fp_x(2),fp0_x(2));
Hhead_yy=interp1(1:fp_x(2),theta(1:fp_x(2)),Hhead_xx,'linear');
Hend_xx=linspace(fp_x(2)+1,n/2+1,n/2+1-fp0_x(2)); 
Hend_yy=interp1(fp_x(2)+1:n/2+1,theta(fp_x(2)+1:n/2+1),Hend_xx,'linear');
theta_ATF=[Hhead_yy,Hend_yy];
theta_ATF=[theta_ATF,-1*theta_ATF(end-1:-1:2)];
theta_ATF=theta_ATF';
end

% Hhead_xx=linspace(1,fp_x,fp0_x);
% Hhead_yy=interp1(1:fp_x,theta(1:fp_x),Hhead_xx,'spline');
% Hend_xx=linspace(fp_x+1,n/2+1,n/2+1-fp0_x); 
% Hend_yy=interp1(fp_x+1:n/2+1,theta(fp_x+1:n/2+1),Hend_xx,'spline');
% theta_ATF=[Hhead_yy,Hend_yy];
% theta_ATF=[theta_ATF,-1*theta_ATF(end-1:-1:2)];
% theta_ATF=theta_ATF';


% figure;plot(theta_ATF);hold on;plot(theta);hold off;
% figure;plot(H_ATF);hold on;plot(H);hold off;



% %%%%%%%%%%%%%%%%%%%%能量校准%%%%%%%%%%%%%%%%%%
% H_ATF_cal=sum(H_ATF.^2)/sum(H.^2)*H;
% H_ATF=H_ATF_cal;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temp=interp(theta(1:fp_x),fp0_x);
% thetahead=temp(1:fp_x:end);
% temp=interp(theta(fp_x+1:n/2+1),n/2+1-fp0_x);
% thetaend=temp(1:(n/2+1-fp_x):end);
% theta_ATF=[thetahead;thetaend];
% theta_ATF=[theta_ATF;-1*theta_ATF(end-1:-1:2)];
% X_ATF=H_ATF.*exp(1i.*theta_ATF);

if length(fp0_x)==2
    X_ATF=H_ATF.*exp(1i.*theta_ATF);
else
    X_ATF=H_ATF.*exp(1i.*theta);
end




% calculate output as the product of the frequency response with the input
% signal in frequency domain.
Xa=fft([S;zeros(length(X)-1,1)]);
X_ATF=fft([ifft(X_ATF);zeros(length(S)-1,1)]);
Xy=X_ATF.*Xa;
y=abs(ifft(Xy));
% y=real(ifft(Xy));
index_neg=find(real(ifft(Xy))<0);
y(index_neg)=y(index_neg)*-1;
y=y(1:length(S));


