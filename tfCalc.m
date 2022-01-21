function [TF,feature]=tfCalc(sig,para,feature)
% This function calculates TFs from all signals in 'align' (unnormalized
% signals) and 'normalize' (normalized signals) step.
% (1)Calculate TFs from all segments;
% (2)Calculate individual average TFs;
% (3)Calculate group average TFs (general transfer function)
% (4)Calculate adaptive transfer function
% results saved in 'tf';
% tf.segment: TFs from all segments;
% tf.individual: individual average TFs;
% tf.general: general TF;
% tf.general_cv: general TFs for leave-one-out cross-validation;
% tf.adaptive: adaptive TF;
% tf.adaptive_cv: adaptive TFs for leave-one-out cross-validation;
% Note that only the peak frequency is saved for the adaptive TFs. The
% frequency response of the adaptive TF is calculated in
% 'sim_peakfrequency.m'.

steps=para.TF.steps;
regressmethod=para.TF.regressmethod;
mulipletmethod=para.TF.mulipletmethod;
regressnum_g=para.TF.regressnum.g;
regressnum_p=para.TF.regressnum.p;
for stepIdx=1:length(steps)
    currentstep=steps{stepIdx};
    na=para.TF.order.(currentstep);
    nb=para.TF.order.(currentstep);
    nk=0;
    
    % calculate individual and average TFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TF.individual.(currentstep).A=nan(para.Nsubj,na+1);
    TF.individual.(currentstep).B=nan(para.Nsubj,nb);
    for subjIdx=1:para.Nsubj
        aorta=sig.(currentstep){subjIdx}.aorta;
        brachial=sig.(currentstep){subjIdx}.brachial;
        L=length(aorta);
        deltL=para.TF.deltL;
        arxdata={};
        if para.TF.useFullPeriodPulseSig==1
            Nsegment=fix((feature.align{subjIdx}.aorta.bottom(end,1)-feature.align{subjIdx}.aorta.bottom(1,1))/deltL*2-1);
            Nperiod=floor(deltL/2/mean(diff(feature.align{subjIdx}.aorta.bottom(:,1))));
            if Nsegment<=1
                x1=(feature.align{subjIdx}.aorta.bottom(1,1):feature.align{subjIdx}.aorta.bottom(end,1));
                arxdata{1}=iddata(aorta(x1),brachial(x1),1/para.fs);
            else
                for ii=1:Nsegment
                    x1=(feature.align{subjIdx}.aorta.bottom((ii-1)*Nperiod+1,1):feature.align{subjIdx}.aorta.bottom((ii-1)*Nperiod+1,1)+deltL);
                    arxdata{ii}=iddata(aorta(x1),brachial(x1),1/para.fs);
                end
            end
        else
            if L<=deltL
                arxdata{1}=iddata(aorta,brachial,1/para.fs);
            else
                for ii=1:floor(L/(deltL/2))-1
                    x1=deltL/2*(ii-1)+1:deltL/2*(ii+1);
                    arxdata{ii}=iddata(aorta(x1),brachial(x1),1/para.fs);
                end
            end
        end
        A=nan(length(arxdata),na+1);
        B=nan(length(arxdata),nb);
        fpe=nan(length(arxdata),1);
        a1=nan(length(arxdata),1);b1=nan(length(arxdata),1);
        a2=nan(length(arxdata),1);b2=nan(length(arxdata),1);c2=nan(length(arxdata),1);
        a3=nan(length(arxdata),1);b3=nan(length(arxdata),1);c3=nan(length(arxdata),1);d3=nan(length(arxdata),1);        
        for segmentIdx=1:length(arxdata)
            m=arx(arxdata{segmentIdx},[na nb nk],'Focus',para.TF.arx.focus); % ARX model for current segment
            TF.segment.(currentstep){subjIdx}{segmentIdx}=m; 
            
            A(segmentIdx,:)=m.A;
            B(segmentIdx,:)=m.B;
            fpe(segmentIdx,:)=m.Report.Fit.FPE;
            
            p1=fittype('a*x+b','independent','x');
            model=tf(m);
            x=lsim(model,arxdata{segmentIdx}.inputdata);
            y=arxdata{segmentIdx}.outputdata;
            f1=fit(x(50:end),y(50:end),p1);
            a1(segmentIdx,:)=f1.a; b1(segmentIdx,:)=f1.b;
            
            p2=fittype('a*x.^2+b*x+c','independent','x');
            f2=fit(x(50:end),y(50:end),p2);
            a2(segmentIdx,:)=f2.a; b2(segmentIdx,:)=f2.b;c2(segmentIdx,:)=f2.c;
            
            p3=fittype('a*x.^3+b*x.^2+c*x+d','independent','x');
            f3=fit(x(50:end),y(50:end),p3);
            a3(segmentIdx,:)=f3.a; b3(segmentIdx,:)=f3.b;c3(segmentIdx,:)=f3.c;d3(segmentIdx,:)=f3.d;
            
%             myfunc = inline('(beta(1)*exp(beta(2)*x))','beta','x');%三个参数分别为：函数模型(注意需要使用点除和点乘)，待定系数，自变量
%             beta0 = [1,1]';%待定系数的预估值
%             beta = nlinfit(x(50:end),y(50:end),myfunc,beta0);
%             a4(segmentIdx,:)=beta(1); b4(segmentIdx,:)=beta(2);
%             yy= a4*exp(b4*x(50:end));
%             plot(x(50:end),yy)
            
        end
        TF.individual.(currentstep).A(subjIdx,:)=mean(A,1);
        TF.individual.(currentstep).B(subjIdx,:)=mean(B,1);
        TF.individual.(currentstep).w(subjIdx,:)=linspace(0,pi,1000);
        TF.individual.(currentstep).f(subjIdx,:)=linspace(0,pi,1000)/pi*para.fs/2;
        TF.individual.(currentstep).H(subjIdx,:)=freqz(TF.individual.(currentstep).B(subjIdx,:),...
            TF.individual.(currentstep).A(subjIdx,:),TF.individual.(currentstep).w(subjIdx,:));
        % calculate peak frequency for each individual TF
%         [~,x1]=min(TF.individual.(currentstep).H(subjIdx,1:160));
        [~,loc_min]=findpeaks(-1*abs(TF.individual.(currentstep).H(subjIdx,:))); x1=min(loc_min);
        [~,loc_min]=findpeaks(-1*phase(TF.individual.(currentstep).H(subjIdx,:))); x2=loc_min(1);
%         [~,loc_max]=findpeaks(phase(TF.individual.(currentstep).H(subjIdx,:)));  
%         if phase(TF.individual.(currentstep).H(subjIdx,loc_max(1)))<0
%            subjIdx
% %            [~,x]=min(abs(phase(TF.individual.(currentstep).H(subjIdx,loc_max:end))));
%         end
%         [~,x]=min(abs(phase(TF.individual.(currentstep).H(subjIdx,loc_min(1):loc_max(1))))); 
        TF.individual.(currentstep).peak_f(subjIdx,:)=[TF.individual.(currentstep).f(subjIdx,x1),TF.individual.(currentstep).f(subjIdx,x2)];
%         TF.individual.(currentstep).peak_f(subjIdx,:)=TF.individual.(currentstep).f(subjIdx,x1);        
        TF.individual.(currentstep).peak_H(subjIdx,1)=TF.individual.(currentstep).H(subjIdx,x1);
        TF.individual.(currentstep).valley_f(subjIdx,1)=TF.individual.(currentstep).f(subjIdx,loc_min(1));
        TF.individual.(currentstep).valley_P(subjIdx,1)=TF.individual.(currentstep).H(subjIdx,loc_min(1));
        
        TF.individual.(currentstep).coefficient1(subjIdx,:)=[mean(a1),mean(b1)];
        TF.individual.(currentstep).coefficient2(subjIdx,:)=[mean(a2),mean(b2),mean(c2)];
        TF.individual.(currentstep).coefficient3(subjIdx,:)=[mean(a3),mean(b3),mean(c3),mean(d3)];
%         TF.individual.(currentstep).coefficient4(subjIdx,:)=[mean(a4),mean(b4)];

    end
    
    % calculate GTFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subjIdxes=1:para.Nsubj;
    for subjIdx=1:para.Nsubj
        % 'general' GTF over all subjects
        % the same for all rows of A/B/w/f/H
        TF.general.(currentstep).A(subjIdx,:)=mean(TF.individual.(currentstep).A,1);
        TF.general.(currentstep).B(subjIdx,:)=mean(TF.individual.(currentstep).B,1);
        TF.general.(currentstep).w(subjIdx,:)=linspace(0,pi,1000);
        TF.general.(currentstep).f(subjIdx,:)=linspace(0,pi,1000)/pi*para.fs/2;
        TF.general.(currentstep).H(subjIdx,:)=freqz(TF.general.(currentstep).B(subjIdx,:),...
            TF.general.(currentstep).A(subjIdx,:),TF.general.(currentstep).w(subjIdx,:));
        [~,x1]=min(TF.general.(currentstep).H(subjIdx,1:160));
        [~,loc_min]=findpeaks(-1*phase(TF.general.(currentstep).H(subjIdx,:))); 
        [~,loc_max]=findpeaks(phase(TF.general.(currentstep).H(subjIdx,:))); 
        [~,x]=min(abs(phase(TF.general.(currentstep).H(subjIdx,loc_min(1):loc_max(1))))); x2=loc_min(1);
        TF.general.(currentstep).peak_f(subjIdx,:)=[TF.general.(currentstep).f(subjIdx,x1),TF.general.(currentstep).f(subjIdx,x2)];
        TF.general.(currentstep).peak_H(subjIdx,1)=[TF.general.(currentstep).H(subjIdx,x1)];
        
        TF.general.(currentstep).coefficient1(subjIdx,:)=mean(TF.individual.(currentstep).coefficient1,1);
        TF.general.(currentstep).coefficient2(subjIdx,:)=mean(TF.individual.(currentstep).coefficient2,1);
        TF.general.(currentstep).coefficient3(subjIdx,:)=mean(TF.individual.(currentstep).coefficient3,1);
%         TF.general.(currentstep).coefficient4(subjIdx,:)=mean(TF.individual.(currentstep).coefficient4,1);
        
        
        % 'general_cv' GTF for leave-one-out cross-validation
        TF.general_cv.(currentstep).A(subjIdx,:)=mean(TF.individual.(currentstep).A(subjIdxes~=subjIdx,:),1);
        TF.general_cv.(currentstep).B(subjIdx,:)=mean(TF.individual.(currentstep).B(subjIdxes~=subjIdx,:),1);
        TF.general_cv.(currentstep).w(subjIdx,:)=linspace(0,pi,1000);
        TF.general_cv.(currentstep).f(subjIdx,:)=linspace(0,pi,1000)/pi*para.fs/2;
        TF.general_cv.(currentstep).H(subjIdx,:)=freqz(TF.general_cv.(currentstep).B(subjIdx,:),...
            TF.general_cv.(currentstep).A(subjIdx,:),TF.general_cv.(currentstep).w(subjIdx,:));
        [~,x1]=min(TF.general_cv.(currentstep).H(subjIdx,1:160));
        [~,loc_min]=findpeaks(-1*phase(TF.general.(currentstep).H(subjIdx,:))); 
        [~,loc_max]=findpeaks(phase(TF.general.(currentstep).H(subjIdx,:))); 
        [~,x]=min(abs(phase(TF.general.(currentstep).H(subjIdx,loc_min(1):loc_max(1))))); x2=loc_min(1);
        TF.general_cv.(currentstep).peak_f(subjIdx,:)=[TF.general_cv.(currentstep).f(subjIdx,x1),TF.general_cv.(currentstep).f(subjIdx,x2)];
        TF.general_cv.(currentstep).peak_H(subjIdx,:)=[TF.general_cv.(currentstep).H(subjIdx,x1),TF.general_cv.(currentstep).H(subjIdx,x2)];
        
        lag(subjIdx,1)=feature.align{subjIdx}.lag;
        
        TF.general_cv.(currentstep).coefficient1(subjIdx,:)=mean(TF.individual.(currentstep).coefficient1(subjIdxes~=subjIdx,:),1);
        TF.general_cv.(currentstep).coefficient2(subjIdx,:)=mean(TF.individual.(currentstep).coefficient2(subjIdxes~=subjIdx,:),1);
        TF.general_cv.(currentstep).coefficient3(subjIdx,:)=mean(TF.individual.(currentstep).coefficient3(subjIdxes~=subjIdx,:),1);
%         TF.general_cv.(currentstep).coefficient4(subjIdx,:)=mean(TF.individual.(currentstep).coefficient4(subjIdxes~=subjIdx,:),1);
    end
    
    % calculate ATFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add here if you want to apply multiple regression analysis or PCA.
    % calculate ATF over all subjects
    featureNames=fields(feature.subjAverage.align.brachial);
    for ii=1:length(featureNames) 
        y1=TF.individual.(currentstep).peak_f(:,1);
%         y1=abs(TF.individual.(currentstep).peak_H(:,1));
        x1=feature.subjAverage.align.brachial.(featureNames{ii});
        [r1,p1]=corrcoef(x1,y1);
        R1(ii,1)=r1(2,1);
        P1(ii,1)=p1(2,1);
        brachial_mat(:,ii)=feature.subjAverage.align.brachial.(featureNames{ii});
        
%         if any(strcmp(featureNames{ii},'HR'))
%             noutliers=11;
%             plotOp=1;
%             [x_cut,y_cut,rSquares,outliers_idx]=regoutliers(x1,y1,noutliers,plotOp); 
%             x_cut = x_cut( ~ isnan(x_cut));
%             y_cut = y_cut( ~ isnan(y_cut));            
%             [r1_cut,p1_cut]=corrcoef(x_cut,y_cut);
%         end        
        
        y2=TF.individual.(currentstep).peak_f(:,2);
%         y2=phase(TF.individual.(currentstep).valley_P(:,1));
        x2=feature.subjAverage.align.brachial.(featureNames{ii});
        [r2,p2]=corrcoef(x2,y2);
        R2(ii,1)=r2(2,1);
        P2(ii,1)=p2(2,1);
        
        [r3,p3]=corrcoef(y1,y2);
        
        y3=TF.individual.(currentstep).valley_f(:,1);
        [r4,p4]=corrcoef(y1,lag);
        [r5,p5]=corrcoef(y2,lag);
        [r6,p6]=corrcoef(y3,lag);
        y4=TF.individual.(currentstep).valley_P(:,1);
        [r7,p7]=corrcoef(phase(y4),lag);
        [r8,p8]=corrcoef(x1,phase(y4));
        R8(ii,1)=r8(2,1);
        P8(ii,1)=p8(2,1);
        
    end
    feature.subjAverage.align.brachial_cormat=corrcoef(brachial_mat);
    
    switch regressmethod
        case 'single'
%     [~,idx]=max(R);  [~,idy]=min(P);  [~,idx]=sort(R);   [~,idy]=sort(P);
        [~,idx1]=min(P1);
%         idx1=1;
        x1=feature.subjAverage.align.brachial.(featureNames{idx1});
        X1=[ones(size(x1)),x1];
        b1=regress(y1,X1);
        TF.adaptive.(currentstep).peak_f(:,1)=X1*b1;
        TF.adaptive.(currentstep).formula=b1;
        case 'multiple'
        [idx1,~]=find(P1<=0.05);
        if isempty(idx1)
           [~,idx1]=min(P1);
            x1=feature.subjAverage.align.brachial.(featureNames{idx1});
            X1=[ones(size(x1)),x1];
            b1=regress(y1,X1);
            TF.adaptive.(currentstep).peak_f(:,1)=X1*b1;
            TF.adaptive.(currentstep).formula1=b1;
        else
            for idx_i1=1:length(idx1)
                x1(:,idx_i1)=feature.subjAverage.align.brachial.(featureNames{idx_i1});
            end 
            switch mulipletmethod
                case 'pca'
                [x1,ps]=mapminmax(x1',0,1);
                x1=x1';
%                 x=zscore(x);%数据标准化
                [COEFF_x,pca_x,LATENT_x]=pca(x1);  %特征向量矩阵
                LATENT_x=100*LATENT_x/sum(LATENT_x);
                if para.TF.outputFlag
%                     figure;pareto(LATENT_x)
                    [r1,p1]=corrcoef(pca_x(:,1),y1);
                    R1(ii,1)=r1(2,1);
                    P1(ii,1)=p1(2,1);
                end
                X1=[ones(size(pca_x(:,1))),pca_x(:,1:regressnum_g)];
                b1=regress(y1,X1);
                TF.adaptive.(currentstep).peak_f(:,1)=X1*b1;
                case 'pls' 
                [x1_std,mu_x1,sigma1_x1]=zscore(x1);%数据标准化
                [y1_std,mu_y1,sigma1_y1]=zscore(y1);
                [~,~,~,~,~,PCTVAR1] = plsregress(x1,y1,size(x1,2));
                ncomp1=find(cumsum(PCTVAR1(1,:))>0.6,1);
%                 [~,~,~,~,b1,~] = plsregress(x1_std,y1,1);
                ncomp1=regressnum_g;
                X1=[ones(size(x1,1),1),x1];
%                 x1(outliers_idx,:)=[];
%                 y1(outliers_idx,:)=[];
%                 TF.adaptive.(currentstep).peak_f(:,1)=X1*b1;
%                 TF.adaptive.(currentstep).formula=b1;
                [betahat1] = bootstrapplsreg(x1,y1,ncomp1,1000,0.05);
                TF.adaptive.(currentstep).formula1=betahat1;
                TF.adaptive.(currentstep).peak_f(:,1)=X1*betahat1.mean;
            end
                

        end     
        
    end
    
    if para.TF.phase
       switch regressmethod
        case 'single'
%     [~,idx]=max(R);  [~,idy]=min(P);  [~,idx]=sort(R);   [~,idy]=sort(P);
        [~,idx2]=min(P2);
        x2=feature.subjAverage.align.brachial.(featureNames{idx2});
        X2=[ones(size(x2)),x2];
        b2=regress(y2,X2);
        TF.adaptive.(currentstep).peak_f(:,2)=X2*b2;
        TF.adaptive.(currentstep).formula2=b2;
         
        case 'multiple'
        [idx2,~]=find(P2<=0.05);
        if isempty(idx2)
           [~,idx2]=min(P2);
            x2=feature.subjAverage.align.brachial.(featureNames{idx2});
            X2=[ones(size(x2)),x2];
            b2=regress(y,X2);
            TF.adaptive.(currentstep).peak_f(:,2)=X2*b2;
            TF.adaptive.(currentstep).formula2=b2;
        else
            for idx_i2=1:length(idx2)
                x2(:,idx_i2)=feature.subjAverage.align.brachial.(featureNames{idx_i2});
            end 
            switch mulipletmethod
                case 'pca'
%                 [x2,ps]=mapminmax(x2',0,1);
%                 x2=x2';
%                 x=zscore(x);%数据标准化
                [COEFF_x,pca_x,LATENT_x]=pca(x2);  %特征向量矩阵
                LATENT_x=100*LATENT_x/sum(LATENT_x);
                if para.TF.outputFlag
%                     figure;pareto(LATENT_x)
                    [r2,p2]=corrcoef(pca_x(:,1),y2);
                    R2(ii,1)=r2(2,1);
                    P2(ii,1)=p2(2,1);
                end
                X2=[ones(size(pca_x(:,1))),pca_x(:,1:regressnum_p)];
                b2=regress(y2,X2);
                TF.adaptive.(currentstep).peak_f(:,2)=X2*b2;
                case 'pls'
                [x2_std,mu_x2,sigma1_x2]=zscore(x2);%数据标准化
                [y2_std,mu_y2,sigma1_y2]=zscore(y2);                    
%                 [x,ps]=mapminmax(x',0,1);
%                 x=x';
                
                [~,~,~,~,~,PCTVAR2] = plsregress(x2_std,y2_std,size(x2,2));
                ncomp2=find(cumsum(PCTVAR2(1,:))>0.6,1);
%                 [~,~,~,~,b2,~] = plsregress(x2_std,y2_std,5);
                ncomp2=regressnum_p;
                X2=[ones(size(x2,1),1),x2];
%                 TF.adaptive.(currentstep).peak_f(:,2)=X2*b2;
%                 TF.adaptive.(currentstep).formula2=b2;
                [betahat2] = bootstrapplsreg(x2,y2,ncomp2,1000,0.05);
                TF.adaptive.(currentstep).formula2=betahat2;
                TF.adaptive.(currentstep).peak_f(:,2)=X2*betahat2.mean;
                
            end
                

        end     
        
        end 
    end
    
       % calculate ATF for leave-one-out cross-validation
        subjIdxes=1:para.Nsubj;
        for subjIdx=1:para.Nsubj
            if strcmp(mulipletmethod,'pls')&&strcmp(regressmethod,'multiple')
%                 [~,~,~,~,b1,~] = plsregress(x1(subjIdxes~=subjIdx,:),y1(subjIdxes~=subjIdx,:),1);
                [betahat1] = bootstrapplsreg(x1(subjIdxes~=subjIdx,:),y1(subjIdxes~=subjIdx,:),ncomp1,1000,0.05);
                X1=[ones(size(x1,1),1),x1];
                TF.adaptive_cv.(currentstep).peak_f(subjIdx,1)=X1(subjIdx,:)*betahat1.mean;
            elseif strcmp(mulipletmethod,'pca')&&strcmp(regressmethod,'multiple')
                   [x1_nor_tra,ps]=mapminmax(x1(subjIdxes~=subjIdx,:)',0,1);
                   [coeff,scoreTrain,~,~,explained,mu]=pca(x1_nor_tra');
                   idx1=find(cumsum(explained)>90,1);
%                    idx1=2;
                   scoreTrain95=scoreTrain(:,1:idx1);
                   X1Tra=[ones(size(scoreTrain95,1),1),scoreTrain95];
                   b1=regress(y1(subjIdxes~=subjIdx),X1Tra); 
                   x1_nor_tes=mapminmax('apply',x1(subjIdx,:)',ps);
    %             X1=[ones(size(x1_nor_tes',1),1),x1_nor_tes'];
    %             TF.adaptive_cv.(currentstep).peak_f(subjIdx,1)=X1*b1;
                  scoreTest95=(x1_nor_tes'-mu)*coeff(:,1:idx1);
                  X1Tes=[ones(size(scoreTest95,1),1),scoreTest95];
                  TF.adaptive_cv.(currentstep).peak_f(subjIdx,1)=X1Tes*b1;
                
            else
                X1=[ones(size(x1)),x1];
                b1=regress(y1(subjIdxes~=subjIdx),X1(subjIdxes~=subjIdx,:));
                TF.adaptive_cv.(currentstep).peak_f(subjIdx,1)=X1(subjIdx,:)*b1;
            end
            
        end
        
       if para.TF.phase 
         for subjIdx=1:para.Nsubj
            if strcmp(mulipletmethod,'pls')&&strcmp(regressmethod,'multiple')  
%                 [~,~,~,~,b2,~]=plsregress(x2(subjIdxes~=subjIdx,:),y2(subjIdxes~=subjIdx,:),1);
                [betahat2] = bootstrapplsreg(x2(subjIdxes~=subjIdx,:),y2(subjIdxes~=subjIdx,:),ncomp2,1000,0.05);
                X2=[ones(size(x2,1),1),x2];
                TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)=X2(subjIdx,:)*betahat2.mean;
            elseif strcmp(mulipletmethod,'pca')&&strcmp(regressmethod,'multiple')
                   [x2_nor_tra,ps]=mapminmax(x2(subjIdxes~=subjIdx,:)',0,1);
                   [coeff,scoreTrain,~,~,explained,mu]=pca(x2_nor_tra');
                   idx2=find(cumsum(explained)>90,1);
%                    idx2=2;
                   scoreTrain95=scoreTrain(:,1:idx2);
                   X2Tra=[ones(size(scoreTrain95,1),1),scoreTrain95];
                   b2=regress(y2(subjIdxes~=subjIdx),X2Tra);
                   x2_nor_tes=mapminmax('apply',x2(subjIdx,:)',ps);          
%             TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)=X2*b2;
%                 TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)=[1,x2(subjIdx,1:end-1),TF.adaptive_cv.(currentstep).peak_f(subjIdx,1)]*b2;
%                 TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)=[1,TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)]*b2;
                 scoreTest95=(x2_nor_tes'-mu)*coeff(:,1:idx2);
                 X2Tes=[ones(size(scoreTest95,1),1),scoreTest95];
                 TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)=X2Tes*b2; 
            
            
            else
                X2=[ones(size(x2)),x2];
%                 [~,~,~,~,b2,~] = plsregress(x2(subjIdxes~=subjIdx,:),y2(subjIdxes~=subjIdx,:));
                b2=regress(y2(subjIdxes~=subjIdx),X2(subjIdxes~=subjIdx,:));
                TF.adaptive_cv.(currentstep).peak_f(subjIdx,2)=X2(subjIdx,:)*b2;
            end
        end  
           
       end
        
        
   
    % output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if para.TF.outputFlag
       outputType=para.TF.outputType;
        figure, % output TFs for all segments
        w=linspace(0,pi,1000);
        f=w/pi*para.fs/2;
        for subjIdx=1:para.Nsubj
            subplot(6,6,subjIdx)
            hold on
            for segmentIdx=1:length(TF.segment.(currentstep){subjIdx})
                H=freqz(TF.segment.(currentstep){subjIdx}{segmentIdx}.B,...
                    TF.segment.(currentstep){subjIdx}{segmentIdx}.A,w);
                plot(f,abs(H))
                peak_f=TF.individual.(currentstep).peak_f(subjIdx,1);
                peak_H=TF.individual.(currentstep).peak_H(subjIdx,:);
                plot(peak_f,abs(peak_H),'o','linewidth',1)
            end
            ylabel('Magnitude')
            xlabel('Frequency(Hz)')
            xlim([0 10])
            ylim([0 2.5])
            title(['Subject ',num2str(subjIdx)])
            set(gca,'Fontsize',8)
        end
        
        figure, % output GTF along with individual TFs
        hold on
        switch outputType
                case 1
                    for subjIdx=1:para.Nsubj
                        plot(TF.individual.(currentstep).f(subjIdx,:),abs(TF.individual.(currentstep).H(subjIdx,:)),'b','linewidth',2)
%                         plot(diff(abs(TF.individual.(currentstep).H(subjIdx,:))),'b','linewidth',2)
                    end
                    plot(TF.general.(currentstep).f(subjIdx,:),abs(TF.general.(currentstep).H(subjIdx,:)),'r-','linewidth',3)
%                     plot(diff(abs(TF.general.(currentstep).H(subjIdx,:))),'r-','linewidth',3)
                    ylabel('Phase')
                    xlabel('Frequency(Hz)')
            %         xlim([0 10])
                    set(gca,'Fontsize',8)
                    
                case 2
                    for subjIdx=1:para.Nsubj
                        plot(TF.individual.(currentstep).f(subjIdx,:),phase(TF.individual.(currentstep).H(subjIdx,:)),'b','linewidth',2)
                    end
                    plot(TF.general.(currentstep).f(subjIdx,:),phase(TF.general.(currentstep).H(subjIdx,:)),'r-','linewidth',3)
                    ylabel('Phase')
                    xlabel('Frequency(Hz)')
            %         xlim([0 10])
                    set(gca,'Fontsize',8)
        end
        
        
%         figure, % output FPE of TF for all segments
%         hold on
%         for subjIdx=1:para.Nsubj
%             for segmentIdx=1:length(TF.segment.(currentstep){subjIdx})
%                 scatter(subjIdx,TF.segment.(currentstep){subjIdx}{segmentIdx}.Report.Fit.FPE,...
%                     'o','markeredgecolor','b','markerfacecolor','b')
%                 
%             end
%             ylabel('FPE')
%             xlabel('Subject')
%         end
        
%         boxplot(TF.segment.(currentstep){subjIdx}{segmentIdx}.Report.Fit.FPE)
        
        figure, % the relationship between TF peak frequency and parameters of pulse waveforms
        for ii=1:length(featureNames)
            subplot(4,4,ii)
            y=TF.individual.(currentstep).peak_f(:,outputType);
%             y=abs(TF.individual.(currentstep).peak_H);
            x1=feature.subjAverage.align.brachial.(featureNames{ii});
            plot(x1,y,'r.')
            [r1,p1]=corrcoef(x1,y);
            R1=r1(2,1);
            P1=p1(2,1);
            title(['R=',num2str(round(R1,2)),' P=',num2str(round(P1,2))]);
%             text([],[],['R=',num2str(round(R,2)),'P=',num2str(round(P,2))]);
            xlabel(featureNames{ii})
            ylabel('F_P (Hz)')
            set(gca,'Fontsize',8)
        end
        figure, % the relationship between estimated and measured peak frequency
        plot(TF.individual.(currentstep).peak_f(:,outputType),TF.adaptive.(currentstep).peak_f(:,outputType),'o');hold on;
        plot(TF.individual.(currentstep).peak_f(:,outputType),TF.individual.(currentstep).peak_f(:,outputType));hold off;
        xlabel('Measured F_P (Hz)')
        ylabel('Estimated F_P (Hz)')
        set(gca,'Fontsize',8)
    end
end   
        

end
