function [sig,feature]=alignsig(sig,para,feature)
% This function aligns all aortic and brachial pulse wave signals. Also,
% all the features points calculated from the 'butterfilter' step are
% shifted in time axis to fit the aligned signals.
% All aligned signals are saved in streucture 'sig.align'.
% All shifted feature points and other features are saved in 'feature.align'

feature.align=feature.butterfilter;
feature.subjAverage.align=feature.subjAverage.butterfilter;
sig.align=sig.butterfilter;
if para.align.switch
for subjIdx=1:para.Nsubj
    deltT=(feature.butterfilter{subjIdx}.aorta.bottom(:,1)-feature.butterfilter{subjIdx}.brachial.bottom(:,1));
%     deltT=(feature.butterfilter{subjIdx}.aorta.onset(:,1)-feature.butterfilter{subjIdx}.brachial.onset(:,1));
    lag=round(mean(deltT));
%     lag=8;
    feature.align{subjIdx}.lag=lag; % time lag between the aortic and brachial pulse waves
%     lag=0;
    if lag>0 
        % aortic pulse wave goes behind the brachial pulse wave
        % cut the head of aortic and the tail of brachial pulse wave 
        sig.align{subjIdx}.aorta=sig.butterfilter{subjIdx}.aorta(lag+1:end);
        sig.align{subjIdx}.brachial=sig.butterfilter{subjIdx}.brachial(1:end-lag);
        featureNames=fields(feature.butterfilter{subjIdx}.aorta);
        for ii=1:length(featureNames)
            if size(feature.align{subjIdx}.aorta.(featureNames{ii}),2)>1
                feature.align{subjIdx}.aorta.(featureNames{ii})(:,1)=feature.align{subjIdx}.aorta.(featureNames{ii})(:,1)-lag;
            end
        end
        feature.ref{subjIdx}=feature.ref{subjIdx}-lag;
    else
        % aortic pulse wave goes beyond the brachial pulse wave
        % cut the head of brachial and the tail of aortic pulse wave 
        lag=abs(lag);
        sig.align{subjIdx}.aorta=sig.butterfilter{subjIdx}.aorta(1:end-lag);
        sig.align{subjIdx}.brachial=sig.butterfilter{subjIdx}.brachial(lag+1:end);
        featureNames=fields(feature.butterfilter{subjIdx}.brachial);
        for ii=1:length(featureNames)
            if size(feature.align{subjIdx}.brachial.(featureNames{ii}),2)>1
                feature.align{subjIdx}.brachial.(featureNames{ii})(:,1)=feature.align{subjIdx}.brachial.(featureNames{ii})(:,1)-lag;
            end
        end
    end
end
end

if para.align.switch && para.align.outputFlag
    subjIdx=para.align.outputSubjIdx;
    figure,
    h1=subplot(211);hold on
    plot(sig.align{subjIdx}.aorta(110:310),'k');axis off;
%     title('Measured Aortic Signal');
    h2=subplot(212);hold on
    plot(sig.align{subjIdx}.brachial(110:310),'k');axis off;
%     title('Measured Brachial Signal');
    linkaxes([h1,h2],'x')
    figure,
    h1=subplot(211);hold on
    plot(sig.align{subjIdx}.aorta,'k');
    scatter(feature.align{subjIdx}.aorta.bottom(:,1),feature.align{subjIdx}.aorta.bottom(:,2),'r');
    title('Alignment: Aortic Signal');
    h2=subplot(212);hold on
    plot(sig.align{subjIdx}.brachial,'k');
    scatter(feature.align{subjIdx}.brachial.bottom(:,1),feature.align{subjIdx}.brachial.bottom(:,2),'r');
    title('Alignment: Brachial Signal');
    linkaxes([h1,h2],'x')
    figure,
    for subjIdx=1:para.Nsubj
        lagMat(subjIdx,1)=feature.align{subjIdx}.lag;
    end
    plot(lagMat/para.fs,'--o')
    ylabel('Lag b/w Pulse Waveforms, s')
    xlabel('Subject')
end