function sig=calibratesig(sig,para,feature,TF)
% This function calibrate the normalized signals with unnormalized signals.

% normalize the estimated aortic pressure waveform and calculate the
% (1)base; (2)foot-to-foot line; (3)scaling factor
sig=normalizesig(sig,para,feature);

% calibrate the estimated normalized aortic pulse waveform
TFs=fields(sig.est);
for tfIdx=1:length(TFs)
    currentTF=TFs{tfIdx};
    for subjIdx=1:para.Nsubj
        switch currentTF
            case 'adaptive_cv'
            aorta=sig.est.(currentTF).normalize{subjIdx};
            base=sig.est.('adaptive_cv').calibratefactor{subjIdx}.base;
            f2fline=sig.est.('adaptive_cv').calibratefactor{subjIdx}.f2fline; 
            scalefactor=sig.est.('adaptive_cv').calibratefactor{subjIdx}.scalefactor;
            sig.est.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor+f2fline+base;
            
            
            base_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_1+f2fline_nolinear_1+base_nolinear_1;
            
            base_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_2+f2fline_nolinear_2+base_nolinear_2;
            
            base_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_3+f2fline_nolinear_3+base_nolinear_3;
            
            
            
            matrix=[(sig.est.(currentTF).calibrate{subjIdx}.^3)';(sig.est.(currentTF).calibrate{subjIdx}.^2)';...
            sig.est.(currentTF).calibrate{subjIdx}';ones(1,length(sig.est.(currentTF).calibrate{subjIdx}))];
        
                
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient1(subjIdx,:)*matrix(3:4,:);
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient2(subjIdx,:)*matrix(2:4,:);
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient3(subjIdx,:)*matrix(1:4,:);
%             sig.est_nolinear_4.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient4(subjIdx,1)*exp(TF.('general_cv').align.coefficient4(subjIdx,2)*matrix(3,:));
            
            case 'general_cv'
            
            aorta=sig.est.(currentTF).normalize{subjIdx};
            base=sig.est.(currentTF).calibratefactor{subjIdx}.base;
            f2fline=sig.est.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor=sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor+f2fline+base;
            
            base_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_1+f2fline_nolinear_1+base_nolinear_1;
            
            base_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_2+f2fline_nolinear_2+base_nolinear_2;
            
            base_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_3+f2fline_nolinear_3+base_nolinear_3;
            
            
            matrix=[(sig.est.(currentTF).calibrate{subjIdx}.^3)';(sig.est.(currentTF).calibrate{subjIdx}.^2)';...
            sig.est.(currentTF).calibrate{subjIdx}';ones(1,length(sig.est.(currentTF).calibrate{subjIdx}))];
        
                
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient1(subjIdx,:)*matrix(3:4,:);
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient2(subjIdx,:)*matrix(2:4,:);
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient3(subjIdx,:)*matrix(1:4,:);
%             sig.est_nolinear_4.(currentTF).calibrate{subjIdx}=TF.('general_cv').align.coefficient4(subjIdx,1)*exp(TF.('general_cv').align.coefficient4(subjIdx,2)*matrix(3,:));
            
            case 'individual'
            
            aorta=sig.est.(currentTF).normalize{subjIdx};
            base=sig.est.(currentTF).calibratefactor{subjIdx}.base;
            f2fline=sig.est.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor=sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor+f2fline+base;
            
            
            base_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_1+f2fline_nolinear_1+base_nolinear_1;
            
            base_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_2+f2fline_nolinear_2+base_nolinear_2;
            
            base_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_3+f2fline_nolinear_3+base_nolinear_3;
            
            
            matrix=[(sig.est.(currentTF).calibrate{subjIdx}.^3)';(sig.est.(currentTF).calibrate{subjIdx}.^2)';...
            sig.est.(currentTF).calibrate{subjIdx}';ones(1,length(sig.est.(currentTF).calibrate{subjIdx}))];
        
                
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=TF.('individual').align.coefficient1(subjIdx,:)*matrix(3:4,:);
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=TF.('individual').align.coefficient2(subjIdx,:)*matrix(2:4,:);
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=TF.('individual').align.coefficient3(subjIdx,:)*matrix(1:4,:);
%             sig.est_nolinear_4.(currentTF).calibrate{subjIdx}=TF.('individual').align.coefficient4(subjIdx,1)*exp(TF.('individual').align.coefficient4(subjIdx,2)*matrix(3,:));
            
            case 'general'
            
            aorta=sig.est.(currentTF).normalize{subjIdx};
            base=sig.est.(currentTF).calibratefactor{subjIdx}.base;
            f2fline=sig.est.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor=sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor+f2fline+base;
            
            base_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_1+f2fline_nolinear_1+base_nolinear_1;
            
            base_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_2+f2fline_nolinear_2+base_nolinear_2;
            
            base_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_3+f2fline_nolinear_3+base_nolinear_3;
            
            
            matrix=[(sig.est.(currentTF).calibrate{subjIdx}.^3)';(sig.est.(currentTF).calibrate{subjIdx}.^2)';...
            sig.est.(currentTF).calibrate{subjIdx}';ones(1,length(sig.est.(currentTF).calibrate{subjIdx}))];
%         
%                 
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient1(subjIdx,:)*matrix(3:4,:);
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient2(subjIdx,:)*matrix(2:4,:);
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient3(subjIdx,:)*matrix(1:4,:);
%             sig.est_nolinear_4.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient4(subjIdx,1)*exp(TF.('general').align.coefficient4(subjIdx,2)*matrix(3,:));
            
            case 'adaptive'
            
            aorta=sig.est.(currentTF).normalize{subjIdx};
            base=sig.est.(currentTF).calibratefactor{subjIdx}.base;
            f2fline=sig.est.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor=sig.est.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor+f2fline+base;
            
            base_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_1=sig.est_nolinear_1.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_1+f2fline_nolinear_1+base_nolinear_1;
            
            base_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_2=sig.est_nolinear_2.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_2+f2fline_nolinear_2+base_nolinear_2;
            
            base_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.base;
            f2fline_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.f2fline; 
            scalefactor_nolinear_3=sig.est_nolinear_3.(currentTF).calibratefactor{subjIdx}.scalefactor;
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=(aorta+1).*scalefactor_nolinear_3+f2fline_nolinear_3+base_nolinear_3;
            
            
            matrix=[(sig.est.(currentTF).calibrate{subjIdx}.^3)';(sig.est.(currentTF).calibrate{subjIdx}.^2)';...
            sig.est.(currentTF).calibrate{subjIdx}';ones(1,length(sig.est.(currentTF).calibrate{subjIdx}))];
        
                
            sig.est_nolinear_1.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient1(subjIdx,:)*matrix(3:4,:);
            sig.est_nolinear_2.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient2(subjIdx,:)*matrix(2:4,:);
            sig.est_nolinear_3.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient3(subjIdx,:)*matrix(1:4,:);
%             sig.est_nolinear_4.(currentTF).calibrate{subjIdx}=TF.('general').align.coefficient4(subjIdx,1)*exp(TF.('general').align.coefficient4(subjIdx,2)*matrix(3,:));
            
            
            
            
        end
        
              
                
        
        
    end
end




