function [sig,feature]=est(sig,para,feature,TF)
% This function estimates aortic pulse waves from normalized ('normalize'
% step) and unnormalized ('align' step) brachial pulse waves using all
% kinds of functions (individual, gerenal, general_cv, adaptive,
% adaptive_cv).

steps=para.est.steps;
TFs=para.est.TFs;
for TFIdx=1:length(TFs)
    currentTF=TFs{TFIdx};
    if any(strcmp(currentTF,{'individual','general','general_cv'}))
        for stepIdx=1:length(steps)
            currentstep=steps{stepIdx};
            for subjIdx=1:para.Nsubj
                model=idpoly(TF.(currentTF).(currentstep).A(subjIdx,:),TF.(currentTF).(currentstep).B(subjIdx,:),1,1,1,[],1/para.fs);
                sig.est.(currentTF).(currentstep){subjIdx}=sim(model,sig.(currentstep){subjIdx}.brachial);
                
                matrix=[(sig.est.(currentTF).(currentstep){subjIdx}.^3)';(sig.est.(currentTF).(currentstep){subjIdx}.^2)';...
                    sig.est.(currentTF).(currentstep){subjIdx}';ones(1,length(sig.est.(currentTF).(currentstep){subjIdx}))];
                
                sig.est_nolinear_1.(currentTF).(currentstep){subjIdx}=(TF.(currentTF).(currentstep).coefficient1(subjIdx,:)*matrix(3:4,:))';
                sig.est_nolinear_2.(currentTF).(currentstep){subjIdx}=(TF.(currentTF).(currentstep).coefficient2(subjIdx,:)*matrix(2:4,:))';
                sig.est_nolinear_3.(currentTF).(currentstep){subjIdx}=(TF.(currentTF).(currentstep).coefficient3(subjIdx,:)*matrix(1:4,:))';
%                 sig.est_nolinear_4.(currentTF).(currentstep){subjIdx}=TF.(currentTF).(currentstep).coefficient4(subjIdx,1)*exp(TF.(currentTF).(currentstep).coefficient4(subjIdx,2)*matrix(3,:));
            end
        end
    elseif any(strcmp(currentTF,{'adaptive','adaptive_cv'}))
        for stepIdx=1:length(steps)
            currentstep=steps{stepIdx};
            for subjIdx=1:para.Nsubj
                sig.est.(currentTF).(currentstep){subjIdx}=sim_peakfrequency(sig.(currentstep){subjIdx}.brachial,...
                    TF.general_cv.(currentstep).A(subjIdx,:),TF.general_cv.(currentstep).B(subjIdx,:),...
                    TF.general_cv.(currentstep).peak_f(subjIdx,:),TF.(currentTF).(currentstep).peak_f(subjIdx,:),para.fs);
             if strcmp(currentTF,'adaptive')
                matrix=[(sig.est.(currentTF).(currentstep){subjIdx}.^3)';(sig.est.(currentTF).(currentstep){subjIdx}.^2)';...
                    sig.est.(currentTF).(currentstep){subjIdx}';ones(1,length(sig.est.(currentTF).(currentstep){subjIdx}))];
                sig.est_nolinear_1.(currentTF).(currentstep){subjIdx}=(TF.('general').(currentstep).coefficient1(subjIdx,:)*matrix(3:4,:))';
                sig.est_nolinear_2.(currentTF).(currentstep){subjIdx}=(TF.('general').(currentstep).coefficient2(subjIdx,:)*matrix(2:4,:))';
                sig.est_nolinear_3.(currentTF).(currentstep){subjIdx}=(TF.('general').(currentstep).coefficient3(subjIdx,:)*matrix(1:4,:))';
%                 sig.est_nolinear_4.(currentTF).(currentstep){subjIdx}=TF.('general').(currentstep).coefficient4(subjIdx,1)*exp(TF.('general').(currentstep).coefficient4(subjIdx,2)*matrix(3,:));
             else
                 matrix=[(sig.est.(currentTF).(currentstep){subjIdx}.^3)';(sig.est.(currentTF).(currentstep){subjIdx}.^2)';...
                    sig.est.(currentTF).(currentstep){subjIdx}';ones(1,length(sig.est.(currentTF).(currentstep){subjIdx}))];
                sig.est_nolinear_1.(currentTF).(currentstep){subjIdx}=(TF.('general_cv').(currentstep).coefficient1(subjIdx,:)*matrix(3:4,:))';
                sig.est_nolinear_2.(currentTF).(currentstep){subjIdx}=(TF.('general_cv').(currentstep).coefficient2(subjIdx,:)*matrix(2:4,:))';
                sig.est_nolinear_3.(currentTF).(currentstep){subjIdx}=(TF.('general_cv').(currentstep).coefficient3(subjIdx,:)*matrix(1:4,:))';
%                 sig.est_nolinear_4.(currentTF).(currentstep){subjIdx}=TF.('general_cv').(currentstep).coefficient4(subjIdx,1)*exp(TF.('general_cv').(currentstep).coefficient4(subjIdx,2)*matrix(3,:));
                
             end
             
             
                
            end
        end
    end
end

% check for unstable systems
% TFs=fields(sig.est);
% for TFIdx=1:length(TFs)
%     currentTF=TFs{TFIdx};
%     for stepIdx=1:length(steps)
%         currentstep=steps{stepIdx};
%         for subjIdx=1:para.Nsubj
%             if max(sig.est.(currentTF).(currentstep){subjIdx})>=500
%                 error('system not stable')
%             end
%         end
%     end
end

