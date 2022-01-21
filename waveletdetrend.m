function [s_detrend,base]=waveletdetrend(s)

[C,L]=wavedec(s,10,'db7');
s_detrend=s-wrcoef('a',C,L,'db7',10);
base=wrcoef('a',C,L,'db7',10);

% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure,hold on,
% subplot(211),
% plot((1:length(s))/100,s/10,'k'),
% plot((1:length(base))/100,base/10,'k','LineWidth',2);
% subplot(212),
% plot((1:length(s_detrend))/100,s_detrend/10,'k');
% ylabel('Amplitude (mmHg)');
% xlabel('Time (s)');

