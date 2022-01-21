function [sig,feature]=filtering(sig,para)
% This function is to filter the raw signals using a ButterWorth filter
% Features are calculated here for the next step 'align'.
switch para.butterfilter.switch
    case 1 % raw signal go through the filter
        [b,a]=butter(para.butterfilter.order,para.butterfilter.Wn);
        for subjIdx=1:para.Nsubj
            temp=filter(b,a,sig.raw{subjIdx}.aorta);
            sig.butterfilter{subjIdx}.aorta=temp(30:end,1);
            temp=filter(b,a,sig.raw{subjIdx}.brachial);
            sig.butterfilter{subjIdx}.brachial=temp(30:end,1);
        end
    case 0 % raw signal does not go through the filter
        sig.butterfilter=sig.raw;
end
% calculate features of all signals
feature=featureCalc(sig,para);

if para.butterfilter.outputFlag && para.butterfilter.switch
    figure,freqz(b,a,(0:0.05:20),para.fs);
    figure, 
    subplot(211),hold on
    plot(sig.butterfilter{1}.aorta,'k');
    plot(sig.raw{1}.aorta,'r');
    title('ButterWorth Filtering: Aortic Signal');
    legend('filtered signal','raw signal');
    subplot(212),hold on
    plot(sig.butterfilter{1}.brachial,'k');
    plot(sig.raw{1}.brachial,'r');
    title('ButterWorth Filtering: Brachial Signal');
    legend('filtered signal','raw signal');
end
