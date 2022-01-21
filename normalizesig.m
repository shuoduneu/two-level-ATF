function sig=normalizesig(sig,para,feature)
% This function 
% (1)normalizes the aortic and brachial pulse signals in 'align' step.
% (Results stored in sig.normalize)
% (2)normalizes the estimated aortic pulse signals that were not normalized
% to derive the baseline, foot-to-foot line, and the scaling factor.
% (Results stored in 'sig.est.(currentTF).calibratefactor')
% 
currentstep=para.currentstep;
normalizemethod=para.normalizemethod;
scalemethod=para.scalemethod;
switch currentstep
    case 'normalize'
        for subjIdx=1:para.Nsubj
            aorta=sig.align{subjIdx}.aorta;
            brachial=sig.align{subjIdx}.brachial;
            
            % detrend and calculate the base
            [aorta,sig.normalize{subjIdx}.aorta_base]=waveletdetrend(aorta);
            [brachial,sig.normalize{subjIdx}.brachial_base]=waveletdetrend(brachial);
%             [aorta,sig.normalize{subjIdx}.aorta_base]=emd_detrend(aorta);
%             [brachial,sig.normalize{subjIdx}.brachial_base]=emd_detrend(brachial);
            
            % calculate features
            featureTemp.aorta=featureCalc_aorta(aorta,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
            featureTemp.brachial=featureCalc_brachial(brachial,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
            
            switch normalizemethod
                case 'f2f'
            % calculate the foot-to-foot line
            [aorta,sig.normalize{subjIdx}.aorta_f2fline]=f2flineCalc(aorta,featureTemp.aorta.bottom);
            [brachial,sig.normalize{subjIdx}.brachial_f2fline]=f2flineCalc(brachial,featureTemp.brachial.bottom);  
                case 'p2p'
            % calculate the foot-to-foot line
            [aorta,sig.normalize{subjIdx}.aorta_f2fline]=f2flineCalc(aorta,featureTemp.aorta.peak);
            [brachial,sig.normalize{subjIdx}.brachial_f2fline]=f2flineCalc(brachial,featureTemp.brachial.peak);
            end
            
            switch scalemethod
                case 'scale'
            % calculate the scaling factor
            [aorta,sig.normalize{subjIdx}.aorta_scalefactor]=scaleCalc(aorta,featureTemp.aorta.bottom);
            [brachial,sig.normalize{subjIdx}.brachial_scalefactor]=scaleCalc(brachial,featureTemp.brachial.bottom);
                case 'noscale'
                    sig.normalize{subjIdx}.aorta_scalefactor=1;
                    sig.normalize{subjIdx}.brachial_scalefactor=1;                    
            end
            
            % calculate the normalized signal
            sig.normalize{subjIdx}.aorta=aorta;
            sig.normalize{subjIdx}.brachial=brachial;
        end
    case 'calibrate'
        % the normalization step here is totally for the calibration.
        TFs=fields(sig.est);
        for tfIdx=1:length(TFs)
            currentTF=TFs{tfIdx};
            for subjIdx=1:para.Nsubj
                aorta=sig.est.(currentTF).align{subjIdx};
                aorta_nolinear_1=sig.est_nolinear_1.(currentTF).align{subjIdx};
                aorta_nolinear_2=sig.est_nolinear_2.(currentTF).align{subjIdx};
                aorta_nolinear_3=sig.est_nolinear_3.(currentTF).align{subjIdx};
%                 
                % detrend and calculate the base
                [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta);
%                 [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.base]=emd_detrend(aorta);
%                 [aorta,sig.normalize{subjIdx}.aorta_base]=emd_detrend(aorta);
                [aorta_nolinear_1,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta_nolinear_1);
                [aorta_nolinear_2,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta_nolinear_2);
                [aorta_nolinear_3,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta_nolinear_3);
%                 
                
                % calculate features
                featureTemp.aorta=featureCalc_aorta(aorta,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
                featureTemp.aorta_nolinear_1=featureCalc_aorta(aorta_nolinear_1,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
                featureTemp.aorta_nolinear_2=featureCalc_aorta(aorta_nolinear_2,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
                featureTemp.aorta_nolinear_3=featureCalc_aorta(aorta_nolinear_3,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
                
                switch normalizemethod
                    case 'f2f'
                % calculate the foot-to-foot line
                [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta,featureTemp.aorta.bottom);
                [aorta_nolinear_1,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_1,featureTemp.aorta_nolinear_1.bottom);
                [aorta_nolinear_2,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_2,featureTemp.aorta_nolinear_2.bottom);
                [aorta_nolinear_3,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_3,featureTemp.aorta_nolinear_3.bottom);
                    case 'p2p'
                % calculate the foot-to-foot line
                [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.f2fline]=p2plineCalc(aorta,featureTemp.aorta.peak);
%                 [aorta_nolinear_1,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_1,featureTemp.aorta_nolinear_1.peak);
%                 [aorta_nolinear_2,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_2,featureTemp.aorta_nolinear_2.peak);
%                 [aorta_nolinear_3,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_3,featureTemp.aorta_nolinear_3.peak);
                end
                
                % calculate the scaling factor
%                 [~,sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta,featureTemp.aorta.bottom);
                switch scalemethod
                    case 'scale'
                % calculate the scaling factor
                [~,sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta,featureTemp.aorta.bottom);
                [~,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta_nolinear_1,featureTemp.aorta_nolinear_1.bottom);
                [~,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta_nolinear_2,featureTemp.aorta_nolinear_2.bottom);
                [~,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta_nolinear_3,featureTemp.aorta_nolinear_3.bottom);
                    case 'noscale'
                sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor=1;
                sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor=1; 
                sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor=1; 
                sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor=1; 
                end
            
                
                % normalized signals are not needed for the calibration.
                

            end
        end
end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,f2fline]=f2flineCalc(s,bottom)
%calculate the foot-to-foot line (set 'Y' of foot to 0)
f2fline=nan(size(s));
bottom_x=bottom(:,1);
bottom_y=bottom(:,2);
f2fline(1:bottom_x(1))=bottom_y(1);
f2fline(bottom_x(end):end)=bottom_y(end);
for ii=1:length(bottom_x)-1
    N=bottom_x(ii+1)-bottom_x(ii)+1;
    f2fline(bottom_x(ii):bottom_x(ii+1))=...
        linspace(bottom_y(ii),bottom_y(ii+1),N);
% X=[bottom_x(ii),bottom_x(ii+1)]; Y=[bottom_y(ii),bottom_y(ii+1)];
%     f2fline(bottom_x(ii):bottom_x(ii+1))=...
%         spline(X,Y,N);
end
% figure;plot(s);hold on;plot(f2fline,'r');hold off;
s=s-f2fline;
% figure;plot(s);
end

function [s,p2pline]=p2plineCalc(s,peak)
%calculate the foot-to-foot line (set 'Y' of foot to 0)
p2pline=nan(size(s));
peak_x=peak(:,1);
peak_y=peak(:,2);
p2pline(1:peak_x(1))=peak_y(1);
p2pline(peak_x(end):end)=peak_y(end);
for ii=1:length(peak_x)-1
    N=peak_x(ii+1)-peak_x(ii)+1;
    p2pline(peak_x(ii):peak_x(ii+1))=...
        linspace(peak_y(ii),peak_y(ii+1),N);
end
s=p2pline-s;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,scalefactor]=scaleCalc(s,bottom)
% calculate scaling factor 'f'
% set mean of each period to 0, foot to -1
bottom_x=bottom(:,1);
scalefactor=nan(size(s));
scalefactor(1:bottom_x(1))=mean(s(bottom_x(1):bottom_x(2)));
scalefactor(bottom_x(end):end)=mean(s(bottom_x(end-1):bottom_x(end)));
% scalefactor(1:bottom_x(1))=mean(s(1:bottom_x(2)));
% scalefactor(bottom_x(end):end)=mean(s(bottom_x(end):s(end)));
for ii=1:length(bottom_x)-1
    scalefactor(bottom_x(ii):bottom_x(ii+1))=mean(s(bottom_x(ii):bottom_x(ii+1)));
end
s=s./scalefactor-1;
% figure;plot(s);
end
% function sig=normalizesig(sig,para,feature)
% % This function 
% % (1)normalizes the aortic and brachial pulse signals in 'align' step.
% % (Results stored in sig.normalize)
% % (2)normalizes the estimated aortic pulse signals that were not normalized
% % to derive the baseline, foot-to-foot line, and the scaling factor.
% % (Results stored in 'sig.est.(currentTF).calibratefactor')
% % 
% currentstep=para.currentstep;
% normalizemethod=para.normalizemethod;
% scalemethod=para.scalemethod;
% switch currentstep
%     case 'normalize'
%         for subjIdx=1:para.Nsubj
%             aorta=sig.align{subjIdx}.aorta;
%             brachial=sig.align{subjIdx}.brachial;
%             
%             % detrend and calculate the base
%             [aorta,sig.normalize{subjIdx}.aorta_base]=waveletdetrend(aorta);
%             [brachial,sig.normalize{subjIdx}.brachial_base]=waveletdetrend(brachial);
% %             [aorta,sig.normalize{subjIdx}.aorta_base]=emd_detrend(aorta);
% %             [brachial,sig.normalize{subjIdx}.brachial_base]=emd_detrend(brachial);
%             
% %             % calculate features
% %             featureTemp.aorta=featureCalc_aorta(aorta,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
% %             featureTemp.brachial=featureCalc_brachial(brachial,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
%             
%             switch normalizemethod
%                 case 'f2f'
%             % calculate the foot-to-foot line
% %             [aorta,sig.normalize{subjIdx}.aorta_f2fline]=f2flineCalc(aorta,featureTemp.aorta.bottom);
% %             [brachial,sig.normalize{subjIdx}.brachial_f2fline]=f2flineCalc(brachial,featureTemp.brachial.bottom); 
%             sig.normalize{subjIdx}.aorta_f2fline=0;
%             sig.normalize{subjIdx}.brachial_f2fline=0;
%             
%                 case 'p2p'
%             % calculate the foot-to-foot line
%             [aorta,sig.normalize{subjIdx}.aorta_f2fline]=f2flineCalc(aorta,featureTemp.aorta.peak);
%             [brachial,sig.normalize{subjIdx}.brachial_f2fline]=f2flineCalc(brachial,featureTemp.brachial.peak);
%             end
%             
%             switch scalemethod
%                 case 'scale'
%             % calculate the scaling factor
% %                     [aorta,sig.normalize{subjIdx}.aorta_scalefactor]=scaleCalc(aorta,featureTemp.aorta.bottom);
% %                     [brachial,sig.normalize{subjIdx}.brachial_scalefactor]=scaleCalc(brachial,featureTemp.brachial.bottom);
%                     sig.normalize{subjIdx}.aorta_scalefactor=1;
%                     sig.normalize{subjIdx}.brachial_scalefactor=1;
%                 case 'noscale'
%                     sig.normalize{subjIdx}.aorta_scalefactor=1;
%                     sig.normalize{subjIdx}.brachial_scalefactor=1;                    
%             end
%             
%             % calculate the normalized signal
%             sig.normalize{subjIdx}.aorta=aorta;
%             sig.normalize{subjIdx}.brachial=brachial;
%         end
%     case 'calibrate'
%         % the normalization step here is totally for the calibration.
%         TFs=fields(sig.est);
%         for tfIdx=1:length(TFs)
%             currentTF=TFs{tfIdx};
%             for subjIdx=1:para.Nsubj
%                 aorta=sig.est.(currentTF).align{subjIdx};
% %                 aorta_nolinear_1=sig.est_nolinear_1.(currentTF).align{subjIdx};
% %                 aorta_nolinear_2=sig.est_nolinear_2.(currentTF).align{subjIdx};
% %                 aorta_nolinear_3=sig.est_nolinear_3.(currentTF).align{subjIdx};
% %                 
%                 % detrend and calculate the base
%                 [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta);
% %                 [aorta_nolinear_1,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta_nolinear_1);
% %                 [aorta_nolinear_2,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta_nolinear_2);
% %                 [aorta_nolinear_3,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base]=waveletdetrend(aorta_nolinear_3);
% %                 
%                 
%                 % calculate features
% %                 featureTemp.aorta=featureCalc_aorta(aorta,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
% %                 featureTemp.aorta_nolinear_1=featureCalc_aorta(aorta_nolinear_1,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
% %                 featureTemp.aorta_nolinear_2=featureCalc_aorta(aorta_nolinear_2,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
% %                 featureTemp.aorta_nolinear_3=featureCalc_aorta(aorta_nolinear_3,feature.align{subjIdx}.aorta.bottom(:,1),para.fs);
%                 
%                 switch normalizemethod
%                     case 'f2f'
%                 % calculate the foot-to-foot line
% %                 [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta,featureTemp.aorta.bottom);
%                 sig.est.(currentTF).calibratefactor{subjIdx}.f2fline=0;
% %                 [aorta_nolinear_1,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_1,featureTemp.aorta_nolinear_1.bottom);
% %                 [aorta_nolinear_2,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_2,featureTemp.aorta_nolinear_2.bottom);
% %                 [aorta_nolinear_3,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_3,featureTemp.aorta_nolinear_3.bottom);
%                     case 'p2p'
%                 % calculate the foot-to-foot line
%                 [aorta,sig.est.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta,featureTemp.aorta.peak);
% %                 [aorta_nolinear_1,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_1,featureTemp.aorta_nolinear_1.peak);
% %                 [aorta_nolinear_2,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_2,featureTemp.aorta_nolinear_2.peak);
% %                 [aorta_nolinear_3,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline]=f2flineCalc(aorta_nolinear_3,featureTemp.aorta_nolinear_3.peak);
%                 end
%                 
%                 % calculate the scaling factor
% %                 [~,sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta,featureTemp.aorta.bottom);
%                 switch scalemethod
%                     case 'scale'
%                 % calculate the scaling factor
% %                 [~,sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta,featureTemp.aorta.bottom);
%                 sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor=1;
% %                 [~,sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta_nolinear_1,featureTemp.aorta_nolinear_1.bottom);
% %                 [~,sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta_nolinear_2,featureTemp.aorta_nolinear_2.bottom);
% %                 [~,sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor]=scaleCalc(aorta_nolinear_3,featureTemp.aorta_nolinear_3.bottom);
%                     case 'noscale'
%                 sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor=1;
% %                 sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor=1; 
% %                 sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor=1; 
% %                 sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor=1; 
%                 end
%             
%                 
%                 % normalized signals are not needed for the calibration.
%                 
% 
%             end
%         end
% end
% end
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [s,f2fline]=f2flineCalc(s,bottom)
% %calculate the foot-to-foot line (set 'Y' of foot to 0)
% f2fline=nan(size(s));
% bottom_x=bottom(:,1);
% bottom_y=bottom(:,2);
% f2fline(1:bottom_x(1))=bottom_y(1);
% f2fline(bottom_x(end):end)=bottom_y(end);
% for ii=1:length(bottom_x)-1
%     N=bottom_x(ii+1)-bottom_x(ii)+1;
%     f2fline(bottom_x(ii):bottom_x(ii+1))=...
%         linspace(bottom_y(ii),bottom_y(ii+1),N);
% % X=[bottom_x(ii),bottom_x(ii+1)]; Y=[bottom_y(ii),bottom_y(ii+1)];
% %     f2fline(bottom_x(ii):bottom_x(ii+1))=...
% %         spline(X,Y,N);
% end
% % figure;plot(s);hold on;plot(f2fline,'r');hold off;
% s=s-f2fline;
% % figure;plot(s);
% end
% 
% function [s,p2pline]=p2plineCalc(s,peak)
% %calculate the foot-to-foot line (set 'Y' of foot to 0)
% p2pline=nan(size(s));
% peak_x=peak(:,1);
% peak_y=peak(:,2);
% p2pline(1:peak_x(1))=peak_y(1);
% p2pline(peak_x(end):end)=peak_y(end);
% for ii=1:length(peak_x)-1
%     N=peak_x(ii+1)-peak_x(ii)+1;
%     p2pline(peak_x(ii):peak_x(ii+1))=...
%         linspace(peak_y(ii),peak_y(ii+1),N);
% end
% s=s-p2pline;
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [s,scalefactor]=scaleCalc(s,bottom)
% % calculate scaling factor 'f'
% % set mean of each period to 0, foot to -1
% bottom_x=bottom(:,1);
% scalefactor=nan(size(s));
% scalefactor(1:bottom_x(1))=mean(s(bottom_x(1):bottom_x(2)));
% scalefactor(bottom_x(end):end)=mean(s(bottom_x(end-1):bottom_x(end)));
% for ii=1:length(bottom_x)-1
%     scalefactor(bottom_x(ii):bottom_x(ii+1))=mean(s(bottom_x(ii):bottom_x(ii+1)));
% end
% s=s./scalefactor-1;
% % figure;plot(s);
% end