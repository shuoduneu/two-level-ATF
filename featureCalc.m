function feature=featureCalc(sig,para,varargin)
% This function would be used in two steps for the calculation of features:
% 'butterfilter': results saved in: 'feature.butterfilter'
% 'est': results saved in: 'feature.interp.est'

if nargin==3
    feature=varargin{1};
end
currentstep=para.currentstep;

switch currentstep
    case 'butterfilter'  % calculate features from aortic and brachial pulse waveforms for 'butterfilter' step
        for subjIdx=1:para.Nsubj
            feature.butterfilter{subjIdx}.aorta=featureCalc_aorta(sig.butterfilter{subjIdx}.aorta,para.fs);
            feature.ref{subjIdx}=feature.butterfilter{subjIdx}.aorta.bottom(:,1);
            feature.butterfilter{subjIdx}.brachial=featureCalc_brachial(sig.butterfilter{subjIdx}.brachial,...
                feature.ref{subjIdx},para.fs);
            featureNames=fields(feature.butterfilter{subjIdx}.aorta);
            for ii=1:length(featureNames)
                if size(feature.butterfilter{subjIdx}.aorta.(featureNames{ii}),2)==1
                    feature.subjAverage.butterfilter.aorta.(featureNames{ii})(subjIdx,1)=...
                        mean(feature.butterfilter{subjIdx}.aorta.(featureNames{ii}));
                end
            end
            featureNames=fields(feature.butterfilter{subjIdx}.brachial);
            for ii=1:length(featureNames)
                if size(feature.butterfilter{subjIdx}.brachial.(featureNames{ii}),2)==1
                    feature.subjAverage.butterfilter.brachial.(featureNames{ii})(subjIdx,1)=...
                        mean(feature.butterfilter{subjIdx}.brachial.(featureNames{ii}));
                end
            end
        end
        
    case 'est' % calculate features from only the aortic pulse waveform for 'est' step
        % interpolate signals to over 1000Hz
        if para.featureCalc.interp
            para.result.interpFactor=ceil(10/para.load.interpFactor); % increase signal sampling rate to over 1000Hz.
            para.result.fs=para.fs*para.result.interpFactor;
            for subjIdx=1:para.Nsubj
                feature.interp.ref{subjIdx}=feature.ref{subjIdx}*para.result.interpFactor;
            end
            [sig_interp,para]=interpsig(sig,para);
        else 
            sig_interp=sig;
            para.result.fs=para.fs;
            for subjIdx=1:para.Nsubj
                feature.interp.ref{subjIdx}=feature.ref{subjIdx};
            end
        end
        
        % calculate features for 'align' step
        for subjIdx=1:para.Nsubj
            feature.interp.align{subjIdx}.aorta=featureCalc_aorta(sig_interp.align{subjIdx}.aorta,...
                feature.interp.ref{subjIdx},para.result.fs);
            feature.interp.align{subjIdx}.brachial=featureCalc_brachial(sig_interp.align{subjIdx}.brachial,...
                feature.interp.ref{subjIdx},para.result.fs);
        end
        
        % calculate features for 'est' step
        TFs=fields(sig_interp.est);
        for tfIdx=1:length(TFs)
            currentTF=TFs{tfIdx};
            steps=fields(sig_interp.est.(currentTF));
            for stepIdx=1:length(steps)
                currentstep=steps{stepIdx};
                if strcmp(currentstep,'calibratefactor'), continue, end
                for subjIdx=1:para.Nsubj
                    feature.interp.est.(currentTF).(currentstep){subjIdx}=...
                        featureCalc_aorta(sig_interp.est.(currentTF).(currentstep){subjIdx},...
                        feature.interp.align{subjIdx}.aorta.bottom(:,1),para.result.fs);
                    
                    feature.interp.est_nolinear_1.(currentTF).(currentstep){subjIdx}=...
                        featureCalc_aorta(sig_interp.est_nolinear_1.(currentTF).(currentstep){subjIdx},...
                        feature.interp.align{subjIdx}.aorta.bottom(:,1),para.result.fs);
                    
                    feature.interp.est_nolinear_2.(currentTF).(currentstep){subjIdx}=...
                        featureCalc_aorta(sig_interp.est_nolinear_2.(currentTF).(currentstep){subjIdx},...
                        feature.interp.align{subjIdx}.aorta.bottom(:,1),para.result.fs);
                    
                    feature.interp.est_nolinear_3.(currentTF).(currentstep){subjIdx}=...
                        featureCalc_aorta(sig_interp.est_nolinear_3.(currentTF).(currentstep){subjIdx},...
                        feature.interp.align{subjIdx}.aorta.bottom(:,1),para.result.fs);
                    
%                     feature.interp.est_nolinear_4.(currentTF).(currentstep){subjIdx}=...
%                         featureCalc_aorta(sig_interp.est_nolinear_4.(currentTF).(currentstep){subjIdx},...
%                         feature.interp.align{subjIdx}.aorta.bottom(:,1),para.result.fs);
                    
                end
            end
        end
end

