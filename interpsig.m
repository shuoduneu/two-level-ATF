function [sig_interp,para]=interpsig(sig,para)
interpFactor=para.result.interpFactor;
% interpolate signals in 'align' step
for subjIdx=1:para.Nsubj
    sig_interp.align{subjIdx}.aorta=interp(sig.align{subjIdx}.aorta,interpFactor);
    sig_interp.align{subjIdx}.brachial=interp(sig.align{subjIdx}.brachial,interpFactor);
end

% interpolate signals in 'est' step
TFs=fields(sig.est);
for tfIdx=1:length(TFs)
    currentTF=TFs{tfIdx};
    steps={'align','calibrate'};
    for stepIdx=1:length(steps)
        currentstep=steps{stepIdx};
        for subjIdx=1:para.Nsubj
            sig_interp.est.(currentTF).(currentstep){subjIdx}=interp(sig.est.(currentTF).(currentstep){subjIdx},interpFactor);
            sig_interp.est_nolinear_1.(currentTF).(currentstep){subjIdx}=interp(sig.est_nolinear_1.(currentTF).(currentstep){subjIdx},interpFactor);
            sig_interp.est_nolinear_2.(currentTF).(currentstep){subjIdx}=interp(sig.est_nolinear_2.(currentTF).(currentstep){subjIdx},interpFactor);
            sig_interp.est_nolinear_3.(currentTF).(currentstep){subjIdx}=interp(sig.est_nolinear_3.(currentTF).(currentstep){subjIdx},interpFactor);
%             sig_interp.est_nolinear_4.(currentTF).(currentstep){subjIdx}=interp(sig.est_nolinear_4.(currentTF).(currentstep){subjIdx},interpFactor);
        end
    end
end