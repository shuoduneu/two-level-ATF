% Main Function
% Structures generated from this analysis 
% 'sig': all signals
% 'para': all parameters for the whole analysis
% 'feature': all features extracted from the signals
% 'tf': all transfer functions
% 'results': estimation error of all features from all signals estimated by all transfer functions
% Steps in this analysis
% load
% butterfilter
% align
% normalize
% tf
% est
% result
% Those step names were used as field names in the structures.
clear, close all; clc
% Load Data
 
% para.dir='C:\Users\win10\Desktop\CABP\RADI data';
para.dir='F:\PW基本处理\数据\中心动脉压与外周动脉压\Radi_aortic&brachial\Radi_aortic&brachial';
para.fs=100; % in Hz
para.load.interpFlag=false; % increase the sampling rate or not
para.load.interpFactor=10; % interpolate factor
[sig,para]=loaddata(para);

% ButterWorth Filtering
 
para.butterfilter.switch=false; % filter the signals or not
para.butterfilter.Wn=0.6; % cutoff of the filter
para.butterfilter.order=4; % order fo the filter
para.butterfilter.outputFlag=true; % output figures or not
para.currentstep='butterfilter';
[sig,feature]=filtering(sig,para);

% Align signals
 
para.align.switch=true; % align the signals or not
para.align.outputFlag=true; % output figures or not
para.align.outputSubjIdx=1; % subject whose results to output
[sig,feature]=alignsig(sig,para,feature);

% Normalize signals 
 
para.currentstep='normalize';
para.scalemethod='scale';  %'scale' or 'noscale'
para.normalizemethod='f2f';  %'p2p' or 'f2f'
sig=normalizesig(sig,para,feature);

% for i=6
% Calculating TFs
 
para.TF.order.align=11; % Model order for unnormalized signals 
para.TF.order.normalize=7; % Model order for normalized signals  
para.TF.deltL=para.fs*10; % segment length
para.TF.outputFlag=true; % output figures or not
para.TF.outputType=1; % 1 for magitude, 2 for phase
para.TF.phase=true; % change phase or not
para.TF.useFullPeriodPulseSig=true; % Segments start from the foot or not 
para.TF.arx.focus='prediction'; % Error to be minimized
para.TF.steps={'align','normalize'}; % Steps from which the signals would be used to calculate TFs
para.TF.regressmethod='multiple'; % 'multiple regression' or 'single regression'
para.TF.mulipletmethod='pls'; % 'pca ' or 'pls' or 'ridge'or 'single'
para.TF.regressnum.g=1; %num of variables for gain
para.TF.regressnum.p=1; %num of variables for phase
[TF,feature]=tfCalc(sig,para,feature);

% Estimation
 
para.est.steps={'align','normalize'};
para.est.TFs={'individual','general','general_cv','adaptive','adaptive_cv'};
sig=est(sig,para,feature,TF);

% Calibration
para.currentstep='calibrate';
sig=calibratesig(sig,para,feature,TF);

% Calculate features

para.currentstep='est';
% para.featureCalc.interp=true;
para.featureCalc.interp=false;
feature=featureCalc(sig,para,feature);

% Results
 
para.result.features={'RMSE','R','SBP','DBP','PP','MBP','ED','AI','pct','FF',}; % features for validation
result=resultCalc(para,feature,sig);

% route1='F:\PW基本处理\调整相角\result\';
% if para.TF.phase
%    route2='相角\'; 
% else
%    route2='幅值\';      
% end
% 
% route='F:\PW基本处理\两级传递函数\result\幅值+相角\';
% save(strcat(route,'feature.mat'),'feature');
% save(strcat(route,'para.mat'),'para');
% save(strcat(route,'result.mat'),'result');
% save(strcat(route,'sig.mat'),'sig');
% save(strcat(route,'TF.mat'),'TF');


% route='F:\PW基本处理\调整相角\result\sensitivity analysis of hyperparameter\num_g\';
% name=num2str(i);
% save(strcat(route,name),'result');
% clear result;clear TF;
% end



