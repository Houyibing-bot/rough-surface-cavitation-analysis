function Result = cavitation33_analysis(data, opts)
% 完整六列逻辑回归：完整模型AUC + 有限改进实验 + 原值k + 条件阈值 + 高清图。
% 在MATLAB导入六列数值矩阵data后运行：Result = cavitation33_analysis(data);
% 也可直接点运行，无参数时读取基础工作区data，结果保存为Result。
% 列顺序：W(m),H(m),k(m^2/s^2),sigma_local,V(0/1),Y(0/1)。
% 基础MATLAB R2021b可用，所有分析和绘图辅助函数已包含在本文件。
% 可选：opts=struct('ShowFigures',true,'CVRepeats',10,'ThresholdBootstrap',500);
%       Result=cavitation33_analysis(data,opts);
% 不清空工作区、不改写data。新图显示并保留；CSV/PDF/PNG/FIG自动保存。

if nargin<1
    if evalin('base','exist(''data'',''var'')')~=1
        error('cav25:MissingData','Import six numeric columns as data, then call Result=cavitation33_analysis(data).');
    end
    data=evalin('base','data');
end
if nargin<2,opts=struct();end
assert(isstruct(opts)&&isscalar(opts),'opts must be a scalar structure.');
opts=defaults(opts,struct('LengthUnit','m','AssociationBootstrap',2000, ...
    'CVRepeats',10,'OuterFolds',5,'InnerFolds',4,'Lambdas',logspace(-4,2,10), ...
    'ThresholdBootstrap',500,'CorrelationBootstrap',2000,'Seed',20261111,'ShowFigures',true, ...
    'FigureWidth',10,'FontSize',14,'ExportDPI',900, ...
    'OutputDir','','RawPressureFile','','CheckWithGlmfit',true, ...
    'MakeFigures',true,'TargetRecall',0.95,'ReferenceVelocity',45,'ReferenceLength',0.05, ...
    'CostRatios',[1 2 5 10],'ScenarioGroups',[3 9 15 30],'ScenarioICC',[0 .01 .025 .05 .1 .2 .5]));
if ~isnumeric(data)||~isreal(data)||~ismatrix(data)||isempty(data)||size(data,2)~=6
    error('cav25:DataShape','data must be a real numeric matrix with exactly six columns W,H,k,sigma,V,Y.');
end
A=double(full(data));assert(all(isfinite(A(:))),'NaN/Inf found; no records are silently removed.');
if strcmpi(opts.LengthUnit,'mm'),A(:,1:2)=A(:,1:2)/1000;
elseif ~strcmpi(opts.LengthUnit,'m'),error('LengthUnit must be m or mm.');end
assert(size(A,1)>=30,'At least 30 records are required for this default workflow.');
assert(all(A(:,1)>0 & A(:,2)>=0 & A(:,3)>=0),'W>0, H>=0 and k>=0 are required.');
assert(all(ismember(A(:,5:6),[0 1]),'all'),'V and Y must use binary 0/1 labels.');
assert(numel(unique(A(:,5)))==2 && numel(unique(A(:,6)))==2,'Both V groups and Y classes are required.');
assert(opts.ThresholdBootstrap>=0&&fix(opts.ThresholdBootstrap)==opts.ThresholdBootstrap,'ThresholdBootstrap must be a nonnegative integer.');
assert(opts.FigureWidth>=7&&opts.FigureWidth<=20,'FigureWidth must be 7 to 20 inches.');
assert(opts.FontSize>=9&&opts.ExportDPI>=150,'Use FontSize>=9 and ExportDPI>=150.');
assert(isscalar(opts.TargetRecall)&&isfinite(opts.TargetRecall)&&opts.TargetRecall>0&&opts.TargetRecall<=1,'TargetRecall must be in (0,1].');
assert(isscalar(opts.ReferenceVelocity)&&isfinite(opts.ReferenceVelocity)&&opts.ReferenceVelocity>0,'ReferenceVelocity must be positive m/s.');
assert(isscalar(opts.ReferenceLength)&&isfinite(opts.ReferenceLength)&&opts.ReferenceLength>0,'ReferenceLength must be positive m.');
assert(all(isfinite(opts.CostRatios))&&all(opts.CostRatios>0),'CostRatios must be positive.');
assert(all(isfinite(opts.ScenarioGroups))&&all(opts.ScenarioGroups>=1)&all(opts.ScenarioGroups<=size(A,1)),'ScenarioGroups must be between 1 and N.');
assert(all(isfinite(opts.ScenarioICC))&&all(opts.ScenarioICC>=0)&all(opts.ScenarioICC<=1),'ScenarioICC must lie in [0,1].');
base=fileparts(mfilename('fullpath'));
if isempty(opts.OutputDir),opts.OutputDir=fullfile(base,['results_',datestr(now,'yyyymmdd_HHMMSS')]);end
out=char(opts.OutputDir);if ~isfolder(out),mkdir(out);end
fprintf('Input: %d records, six columns. Output: %s\n',size(A,1),out);
fprintf('Stage 1/6: maximum-likelihood association analysis with RAW k.\n');
o=struct('BootstrapN',opts.AssociationBootstrap,'Seed',20261109,'MakeFigures',false, ...
    'OutputDir',fullfile(out,'association'),'CheckWithGlmfit',opts.CheckWithGlmfit, ...
    'RawPressureFile',opts.RawPressureFile);
R=cav25_core(A,o,base);
fprintf('Stage 2/6: repeated nested row-wise cross-validation of ridge models.\n');
V=cav25_validate(A,opts,fullfile(out,'validation'));
fprintf('Stage 3/6: conditional thresholds and retuned bootstrap.\n');
T=cav25_thresholds(A,V,opts,fullfile(out,'thresholds'));
fprintf('Stage 4/6: Pearson, Spearman and correlation bootstrap.\n');
C=cav25_correlations(A,opts,fullfile(out,'correlation'));
fprintf('Stage 5/6: AUC-directed tuning and limited coupling candidate.\n');
U=cav25_refinement(A,V,opts,fullfile(out,'refinement'));
Result=struct('refinement',U,'predictorSpecification','raw_k','correlation',C,'data',A,'options',opts,'association',R,'validation',V,'thresholds',T, ...
    'scope','Same-field state classification; row-wise internal diagnostics, not grouped or external validation.');
fprintf('Stage 6/6: complete-model performance and supplementary figures.\n');
if opts.MakeFigures
    cav25_publication_figures(A,R,V,T,out,opts);
    cav25_correlation_figure(C,out,opts);
end
Result.performance=cav25_auc_overview(A,R,V,U,out,opts);
fprintf('Reviewer 3.3 additions: operating points, dimensionless boundaries and design-effect scenarios.\n');
Result.reviewer33=rv33_analysis(A,V,T,U,opts,fullfile(out,'reviewer33'));
save(fullfile(out,'auc_results.mat'),'Result','-v7');
cav25_write_report(Result,out);
fprintf('COMPLETE. Figures, tables and models: %s\n',out);
if nargout==0,assignin('base','Result',Result);end
end

function R = cav25_core(inputData, opts, base)
% Input is passed explicitly from the script workspace.
if nargin<2, opts=struct(); end
opts=defaults(opts,struct('LengthUnit','m','Sheet',1,'BootstrapN',2000, ...
    'Seed',20261109,'UseArchivedIndices',true,'BootstrapIndexFile','', ...
    'RawPressureFile','','OutputDir','','MakeFigures',true, ...
    'CheckWithGlmfit',true));
if isempty(opts.OutputDir)
    opts.OutputDir=fullfile(base,['results_',datestr(now,'yyyymmdd_HHMMSS')]);
end
out=char(opts.OutputDir); if ~isfolder(out), mkdir(out); end
assert(opts.BootstrapN>=0 && opts.BootstrapN==fix(opts.BootstrapN),'BootstrapN must be a nonnegative integer.');
A=double(full(inputData)); % Copy; the workspace variable data is not changed.
assert(isreal(A) && ismatrix(A),'data must be a real two-dimensional numeric matrix.');
assert(size(A,2)==6,'Exactly six numeric columns are required: W,H,k,sigma,V,Y.');
assert(all(isfinite(A(:))),'Missing/non-numeric data found. No rows are silently removed.');
assert(size(A,1)>=12,'At least 12 records are required.');
if strcmpi(opts.LengthUnit,'mm'), A(:,1:2)=A(:,1:2)/1000;
elseif ~strcmpi(opts.LengthUnit,'m'), error('LengthUnit must be m or mm.'); end
assert(all(A(:,1)>0 & A(:,2)>=0 & A(:,3)>=0),'W must be positive; H and k must be nonnegative.');
assert(all(ismember(A(:,5:6),[0,1]),'all'),'V and Y must be binary 0/1.');
y=A(:,6); v=A(:,5); n=numel(y);
assert(numel(unique(y))==2 && numel(unique(v))==2,'Both outcomes and both vortex groups are required.');
R=struct('data',A,'source','Workspace numeric matrix data','options',opts,'n',n, ...
    'positive',sum(y),'negative',sum(1-y),'duplicates',n-size(unique(A,'rows'),1));
R.scope='Contemporaneous conditional association in selected monitoring records; apparent fit only. No causal, temporal, grouped or external validation is claimed.';
R.intervalNote='Wald and Wilson intervals assume independent records. Bootstrap percentile ranges are conditional on successful finite fits, without optimism correction or cluster adjustment.';
F=[1000*A(:,1:2),A(:,3),A(:,4)];
hasLog=all(A(:,4)>0);
if hasLog, F(:,5)=log(A(:,4)); else, F(:,5)=NaN; end
mu=mean(F,1); sd=std(F,0,1);
assert(all(sd(1:4)>0),'A continuous predictor is constant; the specified model cannot be fitted.');
Z=(F-mu)./sd;
R.scaling=struct('mean',mu,'sd',sd,'names',{{'W_mm','H_mm','k_m2s2','sigma_local','ln_sigma'}});
names={'G','GK','Full_raw','Full_log','Interaction_raw'};
xx={Z(:,1:2),Z(:,1:3),[Z(:,1:4),v],[Z(:,1:3),Z(:,5),v],[Z(:,1:4),v,Z(:,4).*v]};
terms={{'z_W','z_H'},{'z_W','z_H','z_k'}, ...
    {'z_W','z_H','z_k','z_sigma','V'}, ...
    {'z_W','z_H','z_k','z_lnsigma','V'}, ...
    {'z_W','z_H','z_k','z_sigma','V','z_sigma_times_V'}};
if ~hasLog
    warning('cav25:NonpositiveSigma','Nonpositive sigma found: Full_log is skipped; raw sigma is retained without shifting or deleting records.');
    names(4)=[];xx(4)=[];terms(4)=[];
end
for j=1:numel(names)
    mj=model_summary(xx{j},y,names{j},terms{j});
    if j==1, R.models=mj; else, R.models(j)=mj; end
    fprintf('%-18s apparent AUC %.8f   Brier %.8f\n',names{j},mj.auc,mj.brier);
end
R.constantBrier=mean(y)*(1-mean(y));

% Descriptive statistics use sample SD and NumPy-compatible type-7 quantiles.
labels={'W_mm';'H_mm';'k_m2s2';'sigma_local';'ln_sigma'};
vals={F(:,1),F(:,2),A(:,3),A(:,4),F(:,5)};
S=NaN(5,7); corrY=NaN(5,1);
for j=1:5
    x=vals{j}; if all(isfinite(x)), S(j,:)=[mean(x),std(x),q7(x,.5),q7(x,.25),q7(x,.75),min(x),max(x)]; c=corrcoef(x,y);corrY(j)=c(1,2); end
end
R.descriptive=array2table(S,'VariableNames',{'Mean','SampleSD','Median','Q1','Q3','Minimum','Maximum'},'RowNames',labels);
R.correlations=table(labels,corrY,'VariableNames',{'Variable','PointBiserial_r'});
cc=corrcoef(v,y);R.phi=cc(1,2);
cnt=[sum(v==0&y==0),sum(v==0&y==1);sum(v==1&y==0),sum(v==1&y==1)];
R.VYcounts=cnt;R.contingencyOR=NaN;R.contingencyOR95=[NaN NaN];
if all(cnt(:)>0)
    lor=log(cnt(2,2)*cnt(1,1)/(cnt(2,1)*cnt(1,2)));
    R.contingencyOR=exp(lor);R.contingencyOR95=exp(lor+[-1 1]*1.959963984540054*sqrt(sum(1./cnt(:))));
end
R.vortexFractions=zeros(2,5);
for vv=0:1
    ss=sum(y(v==vv));nn=sum(v==vv);
    R.vortexFractions(vv+1,:)=[nn,ss,ss/nn,wilson(ss,nn)];
end
edges=q7(A(:,4),[0 1/3 2/3 1]);R.tertileEdges=edges;
rows=zeros(6,9);r=0;
for vv=0:1
    for b=1:3
        r=r+1; sel=v==vv & A(:,4)>=edges(b) & (A(:,4)<edges(b+1) | (b==3 & A(:,4)<=edges(b+1)));
        nn=sum(sel);ss=sum(y(sel));ci=wilson(ss,nn);
        rows(r,:)=[vv,b,edges(b),edges(b+1),nn,ss,ss/max(nn,1),ci];
        if nn==0,rows(r,7)=NaN;end
    end
end
R.tertiles=array2table(rows,'VariableNames',{'V','Bin','Lower','Upper','N','Positive','Fraction','WilsonLow','WilsonHigh'});
R.commonSupport=[max(min(A(v==0,4)),min(A(v==1,4))),min(max(A(v==0,4)),max(A(v==1,4)))];

% Univariate cut-points are not outputs of the multivariable model.
C=zeros(4,10);
for j=1:4
    x=A(:,j); if j<=2,x=x*1000;end
    thresholds=sort(unique(x),'descend');best=-Inf;
    for t=thresholds'
        if j==4, pred=x<=t;else,pred=x>=t;end
        cm=confusion(y,pred);J=cm(5)+cm(6)-1;
        if J>best+1e-12,best=J;C(j,:)=[t,cm,J,sign(4-j)];end
    end
end
R.cuts=array2table(C(:,1:9),'VariableNames',{'Cutpoint','TP','FN','TN','FP','Sensitivity','Specificity','Accuracy','YoudenJ'},'RowNames',{'W_mm','H_mm','k_m2s2','sigma_local'});

% Paired row resamples. Original scaling is fixed; the model is refitted and
% evaluated on EACH SAME resample. This is NOT cross-validation or correction.
[IX,basis]=bootstrap_indices(A,opts,base);B=size(IX,1);M=numel(names);
if opts.UseArchivedIndices && B>0 && startsWith(basis,'MATLAB')
    warning('cav25:NewBootstrap', ...
        'Archived resamples were not reused (data, row order, replicate count or reference files differ). A reproducible new bootstrap is used; report its actual results.');
end
R.bootstrap=struct('basis',basis,'indices',IX,'AUC',NaN(B,M),'Brier',NaN(B,M),'status',zeros(B,M),'failureReason',{cell(B,M)});
for b=1:B
    ii=IX(b,:);yy=y(ii);
    for j=1:M
        fit=mle([ones(n,1),xx{j}(ii,:)],yy);
        if fit.ok
            R.bootstrap.AUC(b,j)=auc_rank(yy,fit.p);R.bootstrap.Brier(b,j)=mean((yy-fit.p).^2);R.bootstrap.status(b,j)=1;
        else
            R.bootstrap.failureReason{b,j}=fit.reason;
        end
    end
    if mod(b,200)==0 || b==B,fprintf('Bootstrap %d/%d\n',b,B);end
end
bs=NaN(M,7);audit=zeros(M,3);rare=find(v==1 & y==0);
for j=1:M
    ok=R.bootstrap.status(:,j)==1;aa=R.bootstrap.AUC(ok,j);bb=R.bootstrap.Brier(ok,j);
    bs(j,:)=[B,sum(ok),sum(~ok),q7(aa,[.025 .975]),q7(bb,[.025 .975])];
    failed=find(~ok);
    for b=failed'
        ii=IX(b,:);audit(j,1)=audit(j,1)+~any(ismember(ii,rare));
        sig=A(ii,4);vv=v(ii);yy=y(ii);pp=sig(vv==1&yy==1);nn=sig(vv==1&yy==0);
        sep=isempty(pp)||isempty(nn)||max(pp)<min(nn)||max(nn)<min(pp);
        audit(j,2)=audit(j,2)+sep;
    end
    audit(j,3)=sum(~ok)-audit(j,2);
end
R.bootstrapSummary=array2table(bs,'VariableNames',{'Requested','FiniteFits','FailedFits','AUCLow','AUCHigh','BrierLow','BrierHigh'},'RowNames',names);
R.bootstrapFailureAudit=array2table(audit,'VariableNames',{'Missing_V1Y0','SigmaSeparatedWithinV1','OtherFailure'},'RowNames',names);
R.rareCellRecordIDs=rare;
pairs={'GK','G';'Full_raw','GK';'Full_log','GK';'Interaction_raw','Full_raw';'Full_log','Full_raw'};br={};
for q=1:size(pairs,1)
    a=find(strcmp(names,pairs{q,1}));b=find(strcmp(names,pairs{q,2}));if isempty(a)||isempty(b),continue;end
    good=R.bootstrap.status(:,a)==1 & R.bootstrap.status(:,b)==1;
    da=R.bootstrap.AUC(good,a)-R.bootstrap.AUC(good,b);db=R.bootstrap.Brier(good,a)-R.bootstrap.Brier(good,b);
    br(end+1,:)=[{pairs{q,1},pairs{q,2},sum(good),R.models(a).auc-R.models(b).auc,R.models(a).brier-R.models(b).brier},num2cell([q7(da,[.025 .975]),q7(db,[.025 .975])])]; %#ok<AGROW>
end
R.pairedDifferences=cell2table(br,'VariableNames',{'Model','Comparator','PairedFinite','DeltaAUC','DeltaBrier','DeltaAUCLow','DeltaAUCHigh','DeltaBrierLow','DeltaBrierHigh'});

% Leave-one-record-out influence: no model selection or data deletion.
lr={};
for j=find(ismember(names,{'Full_raw','Full_log','Interaction_raw'}))
    m=R.models(j);
    for i=1:n
        keep=(1:n)'~=i;f=mle([ones(n-1,1),xx{j}(keep,:)],y(keep));
        da=NaN;av=NaN;bv=NaN;bsig=NaN;bvort=NaN;
        if f.ok,da=max(abs((f.beta-m.beta)./m.se));av=auc_rank(y(keep),f.p);bv=mean((y(keep)-f.p).^2);bsig=f.beta(5);bvort=f.beta(6);end
        lr(end+1,:)={m.name,i,A(i,3),A(i,4),v(i),y(i),f.ok,da,av,bv,bsig,bvort,m.hat(i)}; %#ok<AGROW>
    end
end
R.influence=cell2table(lr,'VariableNames',{'Model','OmittedRecord','k','Sigma','V','Y','FiniteFit','MaxBetaChangeOverFullSE','AUC','Brier','SigmaBeta','VBeta','FullHat'});
[~,R.maxKrecord]=max(A(:,3));
R.pressure=pressure_check(opts.RawPressureFile,A,Z);

% Optional independent verification using MathWorks glmfit, if available.
R.glmfitCheck=table();
if opts.CheckWithGlmfit && exist('glmfit','file')==2 && license('test','Statistics_Toolbox')
    checks={};gopt=statset('glmfit');gopt.MaxIter=1000;gopt.TolX=1e-10;
    for j=1:M
        [b,~,st]=glmfit(xx{j},y,'binomial','link','logit','options',gopt);
        m=R.models(j);err=max(abs(b-m.beta));serr=max(abs(st.se-m.se));
        checks(j,:)={m.name,err,serr,err<1e-6 && serr<1e-5};
    end
    R.glmfitCheck=cell2table(checks,'VariableNames',{'Model','MaxBetaDifference','MaxSEDifference','Passed'});
end
export_results(R,out);
save(fullfile(out,'analysis_results.mat'),'R','-v7');
fprintf('\nCompleted. Results: %s\n',out);
fprintf('Raw-k MLE results are apparent; the next stage reports nested row-wise diagnostics.\n');
end

function o=defaults(o,d)
f=fieldnames(d);for j=1:numel(f),if ~isfield(o,f{j}),o.(f{j})=d.(f{j});end,end
unknown=setdiff(fieldnames(o),f);assert(isempty(unknown),'Unknown option field. See help/README.');
end

function A=read_numeric(file,sheet)
[~,~,ext]=fileparts(file);
if any(strcmpi(ext,{'.xlsx','.xls','.xlsm'})),A=readmatrix(file,'Sheet',sheet,'UseExcel',false);
else,A=readmatrix(file);end
% A textual header can be auto-detected by readmatrix. An entirely nonnumeric
% first row is a header; interior missing rows are errors, never discarded.
while ~isempty(A) && all(isnan(A(1,:))),A(1,:)=[];end
assert(~isempty(A),'No numeric data found.');
end

function f=mle(X,y)
f=struct('ok',false,'reason','iteration_limit_or_separation','beta',[], 'p',[],'cov',[],'ll',NaN,'iterations',0);
if numel(unique(y))<2,f.reason='single_outcome';return;end
if rank(X)<size(X,2),f.reason='rank_deficient';return;end
b=zeros(size(X,2),1);b(1)=log(mean(y)/(1-mean(y)));
for it=1:100
    p=logistic(X*b);info=X'*(X.*(p.*(1-p)));score=X'*(y-p);
    if ~all(isfinite(info(:))) || rcond(info)<eps,f.reason='singular_information_or_separation';return;end
    step=info\score;old=loglik(X,y,b);fac=1;
    while loglik(X,y,b+fac*step)<old-1e-11 && fac>1e-8,fac=fac/2;end
    b=b+fac*step;
    if max(abs(fac*step))<1e-9 && max(abs(score))<1e-7,f.ok=true;break;end
    if ~all(isfinite(b)) || max(abs(b))>1e5,break;end
end
if ~f.ok,return;end
p=logistic(X*b);info=X'*(X.*(p.*(1-p)));
[~,flag]=chol((info+info')/2);
if flag~=0,f.ok=false;f.reason='nonpositive_information';return;end
f.beta=b;f.p=p;f.cov=info\eye(size(info));f.ll=loglik(X,y,b);f.iterations=it;f.reason='finite_MLE';
end
function p=logistic(e),p=exp(-max(0,-e)-log1p(exp(-abs(e))));end
function l=loglik(X,y,b),e=X*b;l=sum(y.*e-max(e,0)-log1p(exp(-abs(e))));end
function a=auc_rank(y,p)
pos=p(y==1);neg=p(y==0)';a=mean((pos>neg)+.5*(pos==neg),'all');
end
function c=confusion(y,pred)
pred=logical(pred);tp=sum(pred&y==1);fn=sum(~pred&y==1);tn=sum(~pred&y==0);fp=sum(pred&y==0);
c=[tp,fn,tn,fp,tp/(tp+fn),tn/(tn+fp),(tp+tn)/numel(y)];
end
function q=q7(x,p)
x=sort(x(:));q=NaN(size(p));if isempty(x),return;end
h=1+(numel(x)-1)*p;lo=floor(h);hi=ceil(h);q=reshape(x(lo),size(p)).*(hi-h)+reshape(x(hi),size(p)).*(h-lo);
same=lo==hi;q(same)=x(lo(same));
end
function ci=wilson(s,n)
if n==0,ci=[NaN NaN];return;end
z=1.959963984540054;p=s/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;r=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d;ci=[c-r,c+r];
end
function m=model_summary(x,y,name,terms)
X=[ones(numel(y),1),x];f=mle(X,y);
assert(f.ok,'Model %s has no finite stable MLE (%s). Do not report a failed fit.',name,f.reason);
m=f;m.name=name;m.terms=[{'Intercept'},terms];m.se=sqrt(diag(f.cov));m.ci=[m.beta-1.959963984540054*m.se,m.beta+1.959963984540054*m.se];
m.OR=exp(m.beta);m.OR95=exp(max(-700,min(700,m.ci)));m.auc=auc_rank(y,m.p);m.brier=mean((y-m.p).^2);m.cm=confusion(y,m.p>=.5);
m.aic=-2*m.ll+2*numel(m.beta);m.vif=NaN(numel(m.beta),1);
for j=1:size(x,2),a=[ones(size(x,1),1),x(:,setdiff(1:size(x,2),j))];res=x(:,j)-a*(a\x(:,j));m.vif(j+1)=sum((x(:,j)-mean(x(:,j))).^2)/sum(res.^2);end
m.hat=sum((X*m.cov).*X,2).*m.p.*(1-m.p);m.maxScore=max(abs(X'*(y-m.p)));m.infoCondition=cond(X'*(X.*(m.p.*(1-m.p))));
[~,order]=sortrows([m.p,(1:numel(y))'],[1 2]);ng=5;nn=floor(numel(y)/ng)*ones(1,ng);nn(1:mod(numel(y),ng))=nn(1:mod(numel(y),ng))+1;
m.calibration=zeros(ng,4);start=0;
for g=1:ng,ii=order(start+(1:nn(g)));start=start+nn(g);m.calibration(g,:)=[mean(m.p(ii)),mean(y(ii)),numel(ii),sum(y(ii))];end
t=sort(unique(m.p),'descend');m.roc=zeros(numel(t)+1,2);
for j=1:numel(t),c=confusion(y,m.p>=t(j));m.roc(j+1,:)=[1-c(6),c(5)];end
end

function [IX,basis]=bootstrap_indices(A,o,base)
B=o.BootstrapN;n=size(A,1);IX=[];basis='MATLAB twister generated row indices';
if ~isempty(o.BootstrapIndexFile)
    IX=readmatrix(o.BootstrapIndexFile);basis='User-provided row-index matrix';
elseif o.UseArchivedIndices && B==2000 && n==105 && isfile(fullfile(base,'reference','records6_numeric.csv')) && isfile(fullfile(base,'reference','bootstrap_indices.csv'))
    ref=readmatrix(fullfile(base,'reference','records6_numeric.csv'));
    if isequal(size(ref),size(A)) && max(abs(ref-A)./(1+abs(ref)),[],'all')<1e-13
        IX=readmatrix(fullfile(base,'reference','bootstrap_indices.csv'));
        basis='Archived 2000 paired row resamples from the manuscript reanalysis (NumPy PCG64 seed 20261109), reused exactly';
    end
end
if isempty(IX)
    old=rng;clean=onCleanup(@()rng(old));rng(o.Seed,'twister');IX=randi(n,B,n); %#ok<NASGU>
end
assert(isequal(size(IX),[B,n]) && all(IX(:)>=1 & IX(:)<=n & IX(:)==fix(IX(:))),'Bootstrap index matrix must be B-by-n with integer indices in 1:n.');
end

function P=pressure_check(file,A,Z)
P=struct('performed',false,'reason','Six columns do not contain local speed or the source pressure. The pressure-offset sensitivity is skipped unless RawPressureFile is supplied.');
if isempty(file),return;end
D=read_numeric(file,1);assert(size(D,1)==size(A,1)&&size(D,2)>=11,'Pressure workbook must have aligned rows and at least A:K.');
assert(all(isfinite(D(:,1:11)),'all'),'Missing pressure workbook values.');
E=D(:,[1 2 6 9 10 11]);assert(max(abs(E-A)./(1+abs(A)),[],'all')<1e-10,'Pressure workbook six columns do not match the main data row by row.');
assert(all(D(:,5)>0),'Local saved velocity values must be positive.');
assert(max(abs(D(:,4)-D(:,3)-100000))<1e-5,'Column D must equal column C + 100000 Pa for this archived convention.');
recomputed=(D(:,4)-2300)./(500*D(:,5).^2);
assert(max(abs(recomputed-A(:,4))./(1+abs(A(:,4))))<1e-10,'Supplied sigma is inconsistent with the saved pressure/speed calculation.');
s0=(D(:,3)-2300)./(500*D(:,5).^2);mu=mean(s0);sd=std(s0);
P.performed=true;P.reason='Archived +100000 Pa convention retained as primary. Omission is sensitivity only, not a corrected pressure datum.';
P.sigma0=s0;P.mean=mu;P.sd=sd;P.reconstructionMaxError=max(abs(recomputed-A(:,4)));
P.model=model_summary([Z(:,1:3),(s0-mu)/sd,A(:,5)],A(:,6),'Pressure_offset_alternative',{'z_W','z_H','z_k','z_sigma0','V'});
P.rawPressure=D(:,3);P.convertedPressure=D(:,4);P.localVelocity=D(:,5);
end

function export_results(R,out)
writeT(R.descriptive,'descriptive');writeT(R.correlations,'correlations');writeT(R.cuts,'univariate_cutpoints');writeT(R.tertiles,'sigma_vortex_tertiles');
writeT(R.bootstrapSummary,'bootstrap_summary');writeT(R.bootstrapFailureAudit,'bootstrap_failure_audit');writeT(R.pairedDifferences,'paired_bootstrap_differences');writeT(R.influence,'leave_one_record_out');
if ~isempty(R.glmfitCheck),writeT(R.glmfitCheck,'glmfit_check');end
writeT(array2table([R.scaling.mean;R.scaling.sd],'VariableNames',R.scaling.names,'RowNames',{'Mean','SampleSD'}),'scaling');
cf={};perf={};cal={};ro={};sc=R.data;bs=(1:size(R.bootstrap.indices,1))';
scoreNames={'W_m','H_m','k_m2s2','Sigma','V','Y'};bootNames={'Replicate'};
for j=1:numel(R.models)
    m=R.models(j);
    for k=1:numel(m.beta)
        or=m.OR(k);ci=m.OR95(k,:);if k==1,or=NaN;ci=[NaN NaN];end
        cf(end+1,:)={m.name,m.terms{k},m.beta(k),m.se(k),m.ci(k,1),m.ci(k,2),or,ci(1),ci(2),m.vif(k)}; %#ok<AGROW>
    end
    perf(end+1,:)=[{m.name,m.auc,m.brier},num2cell(m.cm),{m.aic,m.maxScore,m.infoCondition}]; %#ok<AGROW>
    for k=1:5,cal(end+1,:)=[{m.name,k},num2cell(m.calibration(k,:))];end %#ok<AGROW>
    for k=1:size(m.roc,1),ro(end+1,:)=[{m.name},num2cell(m.roc(k,:))];end %#ok<AGROW>
    sc(:,end+1)=m.p;scoreNames{end+1}=[m.name,'_p']; %#ok<AGROW>
    bs=[bs,R.bootstrap.AUC(:,j),R.bootstrap.Brier(:,j),R.bootstrap.status(:,j)]; %#ok<AGROW>
    bootNames=[bootNames,{[m.name,'_AUC'],[m.name,'_Brier'],[m.name,'_finite']}]; %#ok<AGROW>
end
writeT(cell2table(cf,'VariableNames',{'Model','Term','Beta','SE','BetaLow95','BetaHigh95','OR','ORLow95','ORHigh95','VIF'}),'coefficients');
writeT(cell2table(perf,'VariableNames',{'Model','ApparentAUC','ApparentBrier','TP','FN','TN','FP','Sensitivity','Specificity','Accuracy','AIC','MaxScore','InformationCondition'}),'model_performance');
writeT(cell2table(cal,'VariableNames',{'Model','Bin','MeanProbability','ObservedFraction','N','Positive'}),'calibration');
writeT(cell2table(ro,'VariableNames',{'Model','FalsePositiveRate','TruePositiveRate'}),'roc');
writeT(array2table(sc,'VariableNames',scoreNames),'record_probabilities');
writeT(array2table(bs,'VariableNames',bootNames),'bootstrap_replicates');
writematrix(R.bootstrap.indices,fullfile(out,'bootstrap_indices.csv'));
if R.pressure.performed
    m=R.pressure.model;
    writeT(table(R.pressure.sigma0,m.p,'VariableNames',{'SigmaWithoutOffset','FittedProbability'}),'pressure_sensitivity_records');
    writeT(table(m.terms',m.beta,m.se,'VariableNames',{'Term','Beta','SE'}),'pressure_sensitivity_coefficients');
end
fid=fopen(fullfile(out,'analysis_report.txt'),'w','n','UTF-8');clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Reviewer 2.5 six-column MATLAB analysis: RAW k\nN=%d; Y=1: %d; Y=0: %d; full-row duplicates: %d (retained)\n',R.n,R.positive,R.negative,R.duplicates);
fprintf(fid,'%s\n%s\n\nBootstrap indices: %s\n',R.scope,R.intervalNote,R.bootstrap.basis);
for m=R.models
    fprintf(fid,'\n%s: AUC %.10f; Brier %.10f; TP/FN/TN/FP %d/%d/%d/%d; sensitivity %.6f; specificity %.6f; accuracy %.6f\n',m.name,m.auc,m.brier,m.cm);
end
m=R.models(strcmp({R.models.name},'Full_raw'));i=5;scale=R.scaling.sd(4);
fprintf(fid,'\nFull_raw sigma OR per one ORIGINAL sigma unit: %.8f [%.8f, %.8f]\n',exp(m.beta(i)/scale),exp(m.ci(i,:)/scale));
fprintf(fid,'V-Y phi %.8f; unadjusted OR %.8f [%.8f, %.8f]\n',R.phi,R.contingencyOR,R.contingencyOR95);
for v=0:1,f=R.vortexFractions(v+1,:);fprintf(fid,'V=%d: %d/%d = %.8f; Wilson95 [%.8f, %.8f]\n',v,f(2),f(1),f(3:5));end
fprintf(fid,'Maximum-k record ID: %d (see leave_one_record_out.csv)\n',R.maxKrecord);
fprintf(fid,'\nPressure sensitivity: %s\n',R.pressure.reason);
if R.pressure.performed,fprintf(fid,'Alternative AUC %.10f; Brier %.10f; sigma reconstruction error %.3g; sigma0 range [%.10f, %.10f]\n',R.pressure.model.auc,R.pressure.model.brier,R.pressure.reconstructionMaxError,min(R.pressure.sigma0),max(R.pressure.sigma0));end
fprintf(fid,'\nNo pressure datum is selected by its AUC. No row split, external test, causal effect, independent risk estimate or universal inception threshold is claimed.\n');
    function writeT(t,name)
        writetable(t,fullfile(out,[name,'.csv']),'WriteRowNames',~isempty(t.Properties.RowNames));
    end
end




function V=cav25_validate(A,opts,out)
% Nested, stratified ROW-level internal validation. No surface grouping exists.
d=struct('CVRepeats',10,'OuterFolds',5,'InnerFolds',4,'Lambdas',logspace(-4,2,10),'Seed',20261111);
nf=fieldnames(d);for j=1:numel(nf),if ~isfield(opts,nf{j}),opts.(nf{j})=d.(nf{j});end,end
y=A(:,6);n=numel(y);K=opts.OuterFolds;J=opts.InnerFolds;R=opts.CVRepeats;
assert(all([K,J,R]==fix([K,J,R]))&&K>=2&&J>=2&&R>=1,'CV fold/repeat counts must be integers.');
assert(all(isfinite(opts.Lambdas))&&all(opts.Lambdas>0),'All ridge penalties must be positive.');
opts.Lambdas=unique(opts.Lambdas(:)','sorted');
assert(min(sum(y==0),sum(y==1))>=K,'Too few records per outcome for the requested outer folds.');
old=rng;restore=onCleanup(@()rng(old));rng(opts.Seed,'twister'); %#ok<NASGU>
V=struct('options',opts,'folds',zeros(n,R),'innerFolds',zeros(n,K,R));
V.scope='Repeated stratified nested ROW-level internal validation; no independent surface, temporal or external validation. Repeats are dependent split-sensitivity analyses, not independent samples.';
V.penalty='Mean negative log likelihood + lambda/2 * sum(non-intercept coefficients squared). Continuous standardization is fitted inside each training fold; intercept is unpenalized.';
V.thresholdRule='Inner-CV Youden probability cutpoint; ties choose closest to 0.5, then the larger cutpoint. Outer test labels never choose lambda or threshold.';
names={'G','GK','Full_raw'};if all(A(:,4)>0),names{end+1}='Full_log';end
F=cell(1,numel(names));S=F;FN=F;
blank=struct('name','','oof',NaN(n,R),'oofEta',NaN(n,R),'predictedLabel',false(n,R),'thresholdByRecord',NaN(n,R),'repeatMetrics',table(),'recallPredicted',false(n,R), ...
    'recallThresholdByRecord',NaN(n,R),'recallInnerSensitivityByRecord',NaN(n,R),'final',struct());
V.models=repmat(blank,1,numel(names));
for m=1:numel(names)
    [F{m},S{m},FN{m}]=cv25_features(A,names{m});V.models(m).name=names{m};
end
V.baseline=struct('name','Training_prevalence','oof',NaN(n,R),'oofEta',NaN(n,R),'repeatMetrics',table());
fs={};
for r=1:R
    fold=cv25_stratified_folds(y,K);V.folds(:,r)=fold;
    for k=1:K
        test=fold==k;train=~test;yt=y(train);nt=sum(train);
        assert(min(sum(yt==0),sum(yt==1))>=J,'Too few training records per outcome for inner folds.');
        inner=cv25_stratified_folds(yt,J);V.innerFolds(train,k,r)=inner;
        piTrain=mean(yt);V.baseline.oof(test,r)=piTrain;V.baseline.oofEta(test,r)=log(piTrain/(1-piTrain));
        for m=1:numel(names)
            T=cv25_tune(F{m}(train,:),yt,inner,opts.Lambdas,S{m});
            f=cv25_fit(F{m}(train,:),yt,T.lambda,S{m});
            [pp,ee]=cv25_predict(f,F{m}(test,:));
            V.models(m).oof(test,r)=pp;V.models(m).oofEta(test,r)=ee;
            V.models(m).predictedLabel(test,r)=pp>=T.threshold;
            V.models(m).thresholdByRecord(test,r)=T.threshold;
            [rt,rs]=rv33_recall_cut(yt,T.p,opts.TargetRecall);
            V.models(m).recallPredicted(test,r)=pp>=rt;
            V.models(m).recallThresholdByRecord(test,r)=rt;
            V.models(m).recallInnerSensitivityByRecord(test,r)=rs;
            fs(end+1,:)={names{m},r,k,nt,sum(test),T.lambda,T.threshold,T.loss(T.index),T.youden,sum(f.constantPredictors)}; %#ok<AGROW>
        end
    end
    for m=1:numel(names)
        t=cv25_metrics(y,V.models(m).oof(:,r),V.models(m).predictedLabel(:,r),V.models(m).oofEta(:,r),r);
        V.models(m).repeatMetrics=[V.models(m).repeatMetrics;t];
    end
    t=cv25_metrics(y,V.baseline.oof(:,r),V.baseline.oof(:,r)>=0.5,V.baseline.oofEta(:,r),r);
    V.baseline.repeatMetrics=[V.baseline.repeatMetrics;t];
    fprintf('Nested row-CV repeat %d/%d completed (outer %d, inner %d).\n',r,R,K,J);
end
V.foldSummary=cell2table(fs,'VariableNames',{'Model','Repeat','OuterFold','TrainN','TestN','Lambda','ProbabilityThreshold','InnerOOFLogLoss','InnerOOFYouden','ConstantTrainingPredictors'});
% Fix the final tuning split independently of the number of outer repeats.
rng(opts.Seed+701,'twister');
V.finalInnerFolds=cv25_stratified_folds(y,J);
for m=1:numel(names)
    T=cv25_tune(F{m},y,V.finalInnerFolds,opts.Lambdas,S{m});
    f=cv25_fit(F{m},y,T.lambda,S{m});[f.p,f.eta]=cv25_predict(f,F{m});
    f.threshold=T.threshold;f.featureNames=FN{m};f.innerOOF=T.p;f.innerEta=T.eta;
    [f.recallThreshold,f.innerRecallAtTargetCut]=rv33_recall_cut(y,T.p,opts.TargetRecall);
    f.tuningTable=table(opts.Lambdas(:),T.loss(:),'VariableNames',{'Lambda','InnerOOFLogLoss'});
    f.selectedLambdaIndex=T.index;f.innerYouden=T.youden;
    V.models(m).final=f;
end
ss={};
for m=1:numel(names)+1
    if m<=numel(names),t=V.models(m).repeatMetrics;name=names{m};else,t=V.baseline.repeatMetrics;name=V.baseline.name;end
    ss(end+1,:)={name,mean(t.AUC),std(t.AUC),min(t.AUC),max(t.AUC),mean(t.Brier),std(t.Brier),min(t.Brier),max(t.Brier),mean(t.LogLoss),mean(t.Sensitivity),mean(t.Specificity),mean(t.Accuracy),mean(t.BalancedAccuracy)}; %#ok<AGROW>
end
V.summary=cell2table(ss,'VariableNames',{'Model','MeanAUC','SDAUC','MinAUC','MaxAUC','MeanBrier','SDBrier','MinBrier','MaxBrier','MeanLogLoss','MeanSensitivity','MeanSpecificity','MeanAccuracy','MeanBalancedAccuracy'});
if ~isempty(out),cv25_export(V,A,out);end
end

function [F,standardize,names]=cv25_features(A,name)
F=[1000*A(:,1:2),A(:,3),A(:,4),A(:,5)];
names={'W_mm','H_mm','k_m2s2','sigma_local','V'};
switch name
    case 'G',F=F(:,1:2);names=names(1:2);
    case 'GK',F=F(:,1:3);names=names(1:3);
    case 'Full_raw'
    case 'Full_log'
        assert(all(A(:,4)>0),'Log-sigma sensitivity model requires positive sigma.');
        F(:,4)=log(A(:,4));names{4}='ln_sigma';
    otherwise,error('cv25:UnknownModel','Unknown model %s.',name);
end
standardize=true(1,size(F,2));if size(F,2)==5,standardize(5)=false;end
end

function f=cv25_fit(F,y,lambda,standardize)
mu=zeros(1,size(F,2));sd=ones(1,size(F,2));
mu(standardize)=mean(F(:,standardize),1);rawSD=std(F,0,1);sd(standardize)=rawSD(standardize);
constant=rawSD==0;sd(sd==0)=1;
X=[ones(numel(y),1),(F-mu)./sd];
fit=cv25_solve(X,y,lambda);
if ~fit.ok,error('cv25:FitFailure','Ridge fit failed (%s), lambda=%.6g, n=%d.',fit.status,lambda,numel(y));end
f=struct('beta',fit.beta,'mu',mu,'sd',sd,'lambda',lambda,'standardize',standardize, ...
    'constantPredictors',constant,'iterations',fit.iterations,'gradientNorm',fit.gradientNorm,'status',fit.status);
end

function [p,eta]=cv25_predict(f,F)
eta=[ones(size(F,1),1),(F-f.mu)./f.sd]*f.beta;p=cv25_sigmoid(eta);
end

function T=cv25_tune(F,y,innerFold,lambdas,standardize)
% innerFold can assign all duplicated source records to the same fold.
folds=unique(innerFold(:))';assert(all(folds>=1)&&numel(folds)>=2,'At least two nonempty inner folds are required.');
n=numel(y);L=numel(lambdas);ee=NaN(n,L);pp=NaN(n,L);
for k=folds
    test=innerFold==k;train=~test;
    assert(numel(unique(y(train)))==2,'Each inner training split requires both outcomes.');
    for a=1:L
        f=cv25_fit(F(train,:),y(train),lambdas(a),standardize);
        [pp(test,a),ee(test,a)]=cv25_predict(f,F(test,:));
    end
end
assert(all(isfinite(ee(:))),'Incomplete or nonfinite inner predictions.');
loss=mean(max(ee,0)+log1p(exp(-abs(ee)))-y.*ee,1);
% With numerical objective ties prefer the larger penalty for stability.
lo=min(loss);ids=find(loss<=lo+1e-12);[~,j]=max(lambdas(ids));ix=ids(j);
[threshold,youden]=cv25_cut(y,pp(:,ix));
T=struct('lambda',lambdas(ix),'index',ix,'threshold',threshold,'youden',youden,'loss',loss,'p',pp(:,ix),'eta',ee(:,ix));
end

function fold=cv25_stratified_folds(y,K)
assert(min(sum(y==0),sum(y==1))>=K,'Too few records for stratification.');
fold=zeros(size(y));
for c=0:1
    ii=find(y==c);ii=ii(randperm(numel(ii)));fold(ii)=mod((0:numel(ii)-1)',K)+1;
end
end

function [threshold,J]=cv25_cut(y,p)
candidates=unique([0;0.5;1;p(:)]);v=zeros(size(candidates));
for j=1:numel(candidates)
    pred=p>=candidates(j);v(j)=sum(pred&y==1)/sum(y==1)-sum(pred&y==0)/sum(y==0);
end
J=max(v);ids=find(abs(v-J)<1e-12);dist=abs(candidates(ids)-0.5);ids=ids(abs(dist-min(dist))<1e-12);threshold=max(candidates(ids));
end

function s=cv25_solve(X,y,lambda)
% Convex average-log-loss objective; intercept is never penalized.
n=numel(y);p=size(X,2);pen=lambda*[0;ones(p-1,1)];
s=struct('ok',false,'status','iteration_limit','beta',NaN(p,1),'iterations',0,'gradientNorm',Inf);
if numel(unique(y))<2,s.status='single_outcome';return;end
b=zeros(p,1);b(1)=log(mean(y)/(1-mean(y)));
for it=1:200
    eta=X*b;pr=cv25_sigmoid(eta);g=X'*(pr-y)/n+pen.*b;
    H=X'*(X.*(pr.*(1-pr)))/n+diag(pen);
    if any(~isfinite(H(:)))||rcond(H)<1e-14,s.status='singular_information';return;end
    if norm(g,Inf)<1e-9,s.ok=true;s.status='converged';break;end
    step=H\g;alpha=1;old=cv25_objective(X,y,b,pen);direction=g'*step;
    while alpha>=2^-30
        trial=b-alpha*step;value=cv25_objective(X,y,trial,pen);
        if isfinite(value)&&value<=old-1e-4*alpha*direction+1e-14,break;end
        alpha=alpha/2;
    end
    if alpha<2^-30,s.status='line_search_failure';return;end
    b=trial;
    if any(~isfinite(b))||norm(b,Inf)>1e8,s.status='unbounded_or_nonfinite';return;end
end
s.beta=b;s.iterations=it;s.gradientNorm=norm(g,Inf);
end

function value=cv25_objective(X,y,b,pen)
eta=X*b;value=mean(max(eta,0)+log1p(exp(-abs(eta)))-y.*eta)+0.5*sum(pen.*b.^2);
end

function p=cv25_sigmoid(eta)
p=exp(-max(0,-eta)-log1p(exp(-abs(eta))));
end

function t=cv25_metrics(y,p,pred,eta,repeat)
pos=p(y==1);neg=p(y==0)';auc=mean((pos>neg)+0.5*(pos==neg),'all');
brier=mean((y-p).^2);ll=mean(max(eta,0)+log1p(exp(-abs(eta)))-y.*eta);
c=cv25_confusion(y,pred);c05=cv25_confusion(y,p>=0.5);[ci,cs,cj,status]=cv25_calibration(y,eta);
values=[repeat,auc,brier,ll,ci,cs,cj,c,c05];
columns={'Repeat','AUC','Brier','LogLoss','CalibrationIntercept','CalibrationSlope','CalibrationJointIntercept', ...
    'TP','FN','TN','FP','Sensitivity','Specificity','Accuracy','Precision','NPV','F1','BalancedAccuracy', ...
    'TP_05','FN_05','TN_05','FP_05','Sensitivity_05','Specificity_05','Accuracy_05','Precision_05','NPV_05','F1_05','BalancedAccuracy_05'};
t=array2table(values,'VariableNames',columns);t.CalibrationStatus={status};
end

function c=cv25_confusion(y,pred)
pred=logical(pred);tp=sum(pred&y==1);fn=sum(~pred&y==1);tn=sum(~pred&y==0);fp=sum(pred&y==0);
se=tp/(tp+fn);sp=tn/(tn+fp);precision=NaN;npv=NaN;f1=NaN;
if tp+fp>0,precision=tp/(tp+fp);end
if tn+fn>0,npv=tn/(tn+fn);end
if 2*tp+fp+fn>0,f1=2*tp/(2*tp+fp+fn);end
c=[tp,fn,tn,fp,se,sp,(tp+tn)/numel(y),precision,npv,f1,(se+sp)/2];
end

function [intercept,slope,jointIntercept,status]=cv25_calibration(y,eta)
% Intercept: calibration-in-the-large (slope fixed at 1).
% Slope and joint intercept: unpenalized two-parameter logistic calibration.
intercept=NaN;slope=NaN;jointIntercept=NaN;status='ok';
a=0;ok=false;
for it=1:100
    p=cv25_sigmoid(eta+a);score=sum(y-p);info=sum(p.*(1-p));
    if info<=eps,break;end
    step=score/info;step=max(-10,min(10,step));a=a+step;
    if abs(step)<1e-9,ok=true;break;end
end
if ok,intercept=a;else,status='intercept_not_finite';end
sd=std(eta);mu=mean(eta);
if sd<1e-12,status=[status,';constant_score'];return;end
e0=eta(y==0);e1=eta(y==1);
if max(e0)<=min(e1)||max(e1)<=min(e0),status=[status,';slope_separation'];return;end
f=cv25_solve([ones(numel(y),1),(eta-mu)/sd],y,0);
if f.ok,slope=f.beta(2)/sd;jointIntercept=f.beta(1)-slope*mu;
else,status=[status,';slope_',f.status];end
end

function cv25_export(V,A,out)
if ~isfolder(out),mkdir(out);end
writetable(V.summary,fullfile(out,'nested_cv_summary.csv'));
writetable(V.foldSummary,fullfile(out,'nested_cv_fold_tuning.csv'));
n=size(A,1);R=size(V.folds,2);M=numel(V.models);rows={};inner={};
for r=1:R
    for i=1:n
        rows(end+1,:)={i,r,V.folds(i,r),A(i,6),V.baseline.oof(i,r)}; %#ok<AGROW>
    end
    for k=1:size(V.innerFolds,2)
        ii=find(V.innerFolds(:,k,r)>0);
        for i=ii',inner(end+1,:)={i,r,k,V.innerFolds(i,k,r)};end %#ok<AGROW>
    end
end
t=cell2table(rows,'VariableNames',{'RecordID','Repeat','OuterFold','Y','TrainingPrevalenceProbability'});
writetable(t,fullfile(out,'nested_cv_outer_assignments.csv'));
t=cell2table(inner,'VariableNames',{'RecordID','Repeat','OuterFold','InnerFold'});
writetable(t,fullfile(out,'nested_cv_inner_assignments.csv'));
writematrix([(1:n)',V.finalInnerFolds],fullfile(out,'nested_cv_final_inner_assignments.csv'));
oof={};rep=table();coef={};tuning={};
for m=1:M
    model=V.models(m);tm=model.repeatMetrics;tm.Model=repmat({model.name},height(tm),1);rep=[rep;tm];
    for r=1:R
        for i=1:n
            oof(end+1,:)={model.name,i,r,V.folds(i,r),A(i,6),model.oof(i,r),model.oofEta(i,r),model.thresholdByRecord(i,r),double(model.predictedLabel(i,r))}; %#ok<AGROW>
        end
    end
    f=model.final;
    coef(end+1,:)={model.name,'Intercept',f.beta(1),0,1,f.lambda,f.threshold}; %#ok<AGROW>
    for j=1:numel(f.mu),coef(end+1,:)={model.name,f.featureNames{j},f.beta(j+1),f.mu(j),f.sd(j),f.lambda,f.threshold};end %#ok<AGROW>
    for j=1:height(f.tuningTable),tuning(end+1,:)={model.name,j,f.tuningTable.Lambda(j),f.tuningTable.InnerOOFLogLoss(j),j==f.selectedLambdaIndex};end %#ok<AGROW>
end
tm=V.baseline.repeatMetrics;tm.Model=repmat({V.baseline.name},height(tm),1);rep=[rep;tm];
writetable(rep,fullfile(out,'nested_cv_repeat_metrics.csv'));
t=cell2table(oof,'VariableNames',{'Model','RecordID','Repeat','OuterFold','Y','OOFProbability','OOFLogit','TrainingSelectedThreshold','PredictedLabel'});
writetable(t,fullfile(out,'nested_cv_oof_predictions.csv'));
t=cell2table(coef,'VariableNames',{'Model','Term','Beta','TrainingMean','TrainingSD','Lambda','FinalInnerCVThreshold'});
writetable(t,fullfile(out,'ridge_final_coefficients.csv'));
t=cell2table(tuning,'VariableNames',{'Model','Candidate','Lambda','InnerOOFLogLoss','Selected'});
writetable(t,fullfile(out,'ridge_final_tuning.csv'));
fid=fopen(fullfile(out,'nested_cv_scope.txt'),'w');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n%s\n%s\n',V.scope,V.penalty,V.thresholdRule);
fprintf(fid,'ROC and calibration display all repeats; the confusion figure uses prespecified repeat 1. Repeated-CV ranges are split sensitivity, not 95%% confidence intervals.\n');
fprintf(fid,'Final full-data coefficients and threshold are for subsequent predictions, not held-out test-performance estimation.\n');
fprintf(fid,'CalibrationIntercept fixes slope=1; CalibrationSlope and CalibrationJointIntercept come from a joint unpenalized calibration fit.\n');
end

function T=cav25_thresholds(A,V,opts,out)
% Conditional crossings of the final full ridge model; not physical inception.
if ~isfolder(out),mkdir(out);end
j=find(strcmp({V.models.name},'Full_raw'),1);f=V.models(j).final;
[F,standardize]=cv25_features(A,'Full_raw');n=size(A,1);
T.probabilityThreshold=f.threshold;
T.referenceFeatures=[mean(F(:,1:3),1),mean(F(:,4)),0];
T.referenceRaw=[T.referenceFeatures(1:2)/1000,T.referenceFeatures(3)];
T.sigmaSupport=[min(A(A(:,5)==0,4)),max(A(A(:,5)==0,4));min(A(A(:,5)==1,4)),max(A(A(:,5)==1,4))];
T.commonSigma=[max(T.sigmaSupport(:,1)),min(T.sigmaSupport(:,2))];
names={'W_mm','H_mm','k_m2s2','sigma_local'};units={'mm','mm','m^2/s^2','1'};
rows=cell(n*4,9);counter=0;
for i=1:n
    for q=1:4
        [x,direction]=th25_crossing(f,F(i,:),q,f.threshold);
        group=A(:,5)==A(i,5);
        raw=F(group,q);
        within=isfinite(x)&&x>=min(raw)&&x<=max(raw);
        domain=isfinite(x)&&((q==1&&x>0)||(ismember(q,[2 3])&&x>=0)||q==4);
        counter=counter+1;rows(counter,:)={i,names{q},units{q},f.threshold,x,direction,within,domain,A(i,5)};
    end
end
T.recordThresholds=cell2table(rows,'VariableNames',{'Record','Variable','Unit','ProbabilityCutoff','ConditionalCrossing','PositiveWhen','WithinVGroupMarginalRange','PhysicalDomainValid','V'});
writetable(T.recordThresholds,fullfile(out,'conditional_thresholds_each_record.csv'));
ref=cell(4,7);row=0;
for v=0:1
    point=T.referenceFeatures;point(5)=v;
    for cutoff=[.5,f.threshold]
        [x,direction]=th25_crossing(f,point,4,cutoff);row=row+1;
        ref(row,:)={v,cutoff,x,direction,isfinite(x)&&x>=T.sigmaSupport(v+1,1)&&x<=T.sigmaSupport(v+1,2),T.sigmaSupport(v+1,1),T.sigmaSupport(v+1,2)};
    end
end
T.referenceCrossings=cell2table(ref,'VariableNames',{'V','ProbabilityCutoff','SigmaCrossing','PositiveWhen','WithinVGroupMarginalRange','ObservedMin','ObservedMax'});
writetable(T.referenceCrossings,fullfile(out,'reference_sigma_crossings.csv'));

% Bootstrap the complete tuning procedure, keeping replicate copies of an
% original record in the SAME inner fold. This is not physical grouping.
B=opts.ThresholdBootstrap;old=rng;cleanup=onCleanup(@()rng(old));rng(opts.Seed+904,'twister'); %#ok<NASGU>
boot=NaN(B,7);coef=NaN(B,6);idx=NaN(B,n);folds=NaN(B,n);status=false(B,1);reasons=cell(B,1);
for b=1:B
    ii=randi(n,n,1);idx(b,:)=ii';[ids,~,map]=unique(ii);
    K=min(opts.InnerFolds,min(sum(A(ids,6)==0),sum(A(ids,6)==1)));
    if K<2,reasons{b}='Too few unique records in an outcome class';continue;end
    g=cv25_stratified_folds(A(ids,6),K);inner=g(map);folds(b,:)=inner';
    % Assert resampled copies never leak across the inner folds.
    for id=ids',assert(numel(unique(inner(ii==id)))==1);end
    try
        tune=cv25_tune(F(ii,:),A(ii,6),inner,opts.Lambdas,standardize);
        fit=cv25_fit(F(ii,:),A(ii,6),tune.lambda,standardize);
        if isfield(fit,'ok') && ~fit.ok,error('cav25:BootstrapFit','Ridge fit did not converge');end
        p0=T.referenceFeatures;p0(5)=0;p1=p0;p1(5)=1;
        boot(b,:)=[tune.lambda,tune.threshold,th25_crossing(fit,p0,4,tune.threshold), ...
            th25_crossing(fit,p1,4,tune.threshold),th25_crossing(fit,p0,4,.5),th25_crossing(fit,p1,4,.5),fit.beta(5)];
        % Coefficients on the original full-data scaling, including intercept.
        native=fit.beta(2:end)'./fit.sd;
        coef(b,:)=[fit.beta(1)+sum(native.*(f.mu-fit.mu)),native.*f.sd];
        boot(b,7)=coef(b,5); % Use a fixed scale across bootstrap fits.
        status(b)=true;
    catch err
        reasons{b}=[err.identifier,': ',err.message];
    end
    if mod(b,50)==0||b==B,fprintf('Threshold bootstrap %d/%d\n',b,B);end
end
T.bootstrap=struct('values',boot,'coefficients',coef,'indices',idx,'innerFolds',folds,'success',status,'failureReason',{reasons});
vars={'Lambda','ProbabilityCutoff','SigmaV0Tuned','SigmaV1Tuned','SigmaV0Half','SigmaV1Half','SigmaBetaOnOriginalScaling'};
tab=array2table([(1:B)',double(status),boot],'VariableNames',[{'Replicate','Success'},vars]);
writetable(tab,fullfile(out,'threshold_bootstrap_replicates.csv'));
writematrix(idx,fullfile(out,'threshold_bootstrap_source_indices.csv'));
writematrix(folds,fullfile(out,'threshold_bootstrap_inner_folds.csv'));
sumrows=cell(7,6);
for q=1:7
    x=boot(status & isfinite(boot(:,q)),q);bounds=q7(x,[.025,.5,.975]);
    sumrows(q,:)={vars{q},numel(x),B-numel(x),bounds(1),bounds(2),bounds(3)};
end
T.bootstrapSummary=cell2table(sumrows,'VariableNames',{'Quantity','FiniteValues','MissingOrFailed','Percentile2_5','Median','Percentile97_5'});
writetable(T.bootstrapSummary,fullfile(out,'threshold_bootstrap_summary.csv'));
ci=NaN(6,2);for q=1:6,ci(q,:)=q7(coef(status,q),[.025 .975]);end
T.ridgeCoefficientRanges=table([{'Intercept'};f.featureNames(:)],f.beta,ci(:,1),ci(:,2), ...
    'VariableNames',{'Term','Beta','BootstrapPercentile2_5','BootstrapPercentile97_5'});
writetable(T.ridgeCoefficientRanges,fullfile(out,'ridge_coefficients_bootstrap.csv'));
T.note='Thresholds are conditional statistical crossings. Within a marginal range does not establish joint support. Bootstrap ranges use row independence, include retuning, and are conditional on finite crossings; no group validation or physical-inception uncertainty is implied.';
if sum(status)<B,warning('cav25:ThresholdBootstrapFailures','%d of %d threshold bootstrap fits failed; see saved failureReason.',B-sum(status),B);end
end

function [x,direction]=th25_crossing(f,point,q,cutoff)
direction='none';x=NaN;if cutoff<=0||cutoff>=1||abs(f.beta(q+1))<1e-10,return;end
z=(point-f.mu)./f.sd;other=f.beta(1)+z*f.beta(2:end)-z(q)*f.beta(q+1);
t=f.mu(q)+f.sd(q)*(log(cutoff/(1-cutoff))-other)/f.beta(q+1);
x=t;
if f.beta(q+1)>0,direction='greater_or_equal';else,direction='less_or_equal';end
end

function C=cav25_correlations(A,opts,out)
% Raw-variable correlation analysis. No logarithmic transformation of k.
% Paired row bootstrap ranges adopt working independence of the records.
if ~isfield(opts,'CorrelationBootstrap'),opts.CorrelationBootstrap=2000;end
if ~isfield(opts,'Seed'),opts.Seed=20261111;end
B=opts.CorrelationBootstrap;
assert(isscalar(B)&&isfinite(B)&&B>=0&&B==fix(B),'CorrelationBootstrap must be a nonnegative integer.');
assert(isnumeric(A)&&ismatrix(A)&&size(A,2)==6&&all(isfinite(A(:))),'Six finite numeric columns are required.');
if ~isfolder(out),mkdir(out);end
X=[1000*A(:,1:2),A(:,3:6)];n=size(X,1);
C=struct();C.names={'W_mm','H_mm','k_m2s2','sigma_c','V','Y'};
C.labels={'W','H','k','\sigma_c','V','Y'};
C.Pearson=cr25_corr(X);C.Spearman=cr25_corr(cr25_ranks(X));
C.sampleSize=n;C.bootstrapRequested=B;
C.intervalNote=['Paired row-bootstrap percentile ranges assume independent records and are conditional on finite correlations. ' ...
    'Continuous-variable Pearson correlations with binary Y are point-biserial; the V-Y Pearson correlation is phi. ' ...
    'Spearman correlations use average ranks for ties. Correlation is neither a causal effect nor validation of prediction.'];
old=rng;restore=onCleanup(@()rng(old));rng(opts.Seed+1901,'twister'); %#ok<NASGU>
bp=NaN(B,5);bs=NaN(B,5);indices=zeros(B,n);single=false(B,1);constant=false(B,6);
for b=1:B
    ii=randi(n,n,1);indices(b,:)=ii';XX=X(ii,:);
    constant(b,:)=max(XX,[],1)==min(XX,[],1);single(b)=constant(b,6);
    cp=cr25_corr(XX);cs=cr25_corr(cr25_ranks(XX));
    bp(b,:)=cp(1:5,6)';bs(b,:)=cs(1:5,6)';
    if mod(b,500)==0||b==B,fprintf('Correlation bootstrap %d/%d\n',b,B);end
end
C.bootstrap=struct('Pearson',bp,'Spearman',bs,'indices',indices,'singleOutcome',single,'constantColumns',constant);
pc=NaN(5,2);sc=NaN(5,2);np=zeros(5,1);ns=zeros(5,1);
for j=1:5
    pp=bp(isfinite(bp(:,j)),j);ss=bs(isfinite(bs(:,j)),j);
    pc(j,:)=cr25_quantile(pp,[.025 .975]);sc(j,:)=cr25_quantile(ss,[.025 .975]);
    np(j)=numel(pp);ns(j)=numel(ss);
end
C.summary=table(C.names(1:5)',C.Pearson(1:5,6),pc(:,1),pc(:,2),np,B-np, ...
    C.Spearman(1:5,6),sc(:,1),sc(:,2),ns,B-ns, ...
    'VariableNames',{'Predictor','Pearson','PearsonLower','PearsonUpper','PearsonFinite','PearsonMissing', ...
    'Spearman','SpearmanLower','SpearmanUpper','SpearmanFinite','SpearmanMissing'});
pt=array2table(C.Pearson,'VariableNames',C.names,'RowNames',C.names);
st=array2table(C.Spearman,'VariableNames',C.names,'RowNames',C.names);
writetable(pt,fullfile(out,'pearson_matrix.csv'),'WriteRowNames',true);
writetable(st,fullfile(out,'spearman_matrix.csv'),'WriteRowNames',true);
writetable(C.summary,fullfile(out,'correlations_with_cavitation.csv'));
bn=[strcat('Pearson_',C.names(1:5)),strcat('Spearman_',C.names(1:5))];
bt=array2table([(1:B)',double(single),bp,bs],'VariableNames',[{'Replicate','SingleOutcome'},bn]);
writetable(bt,fullfile(out,'correlation_bootstrap_replicates.csv'));
writematrix(indices,fullfile(out,'correlation_bootstrap_indices.csv'));
writetable(array2table(double(constant),'VariableNames',C.names),fullfile(out,'correlation_bootstrap_constant_columns.csv'));
fid=fopen(fullfile(out,'correlation_scope.txt'),'w','n','UTF-8');assert(fid>0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',C.intervalNote);
fprintf(fid,'Variable order: W(mm), H(mm), raw k(m^2/s^2), sigma_c, visual V, visual Y.\n');
fprintf(fid,'Requested bootstrap samples: %d; single-outcome samples: %d.\n',B,sum(single));
fprintf(fid,'No missing rows are deleted. Constant resampled variables produce NaN correlations and are counted explicitly.\n');
fprintf(fid,'Ranges do not account for unknown surface/time/protrusion clusters. No p-values or causal significance claims are generated.\n');
end

function cav25_correlation_figure(C,out,opts)
% Two publication figures with vector rectangles and no explanatory footnotes.
if ~isfield(opts,'FigureWidth'),opts.FigureWidth=10;end
if ~isfield(opts,'FontSize'),opts.FontSize=14;end
if ~isfield(opts,'ExportDPI'),opts.ExportDPI=900;end
if ~isfield(opts,'ShowFigures'),opts.ShowFigures=true;end
if ~isfolder(out),mkdir(out);end
f=cr25_figure('Correlation matrices',5.4,opts);
matrices={C.Pearson,C.Spearman};positions=[.075 .235 .39 .69;.565 .235 .39 .69];
titles={'(a) Pearson correlation','(b) Spearman correlation'};
for k=1:2
    ax=axes('Parent',f,'Position',positions(k,:));hold(ax,'on');M=matrices{k};
    for r=1:6
        for c=1:6
            v=M(r,c);color=cr25_color(v);
            rectangle(ax,'Position',[c-.5,r-.5,1,1],'FaceColor',color,'EdgeColor',[1 1 1],'LineWidth',.8);
            tc=[.12 .14 .17];if isfinite(v)&&abs(v)>.55,tc=[1 1 1];end
            label='--';if isfinite(v),label=sprintf('%.2f',v);end
            text(ax,c,r,label,'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontName','Times New Roman','FontSize',opts.FontSize-1,'Color',tc);
        end
    end
    set(ax,'XLim',[.5 6.5],'YLim',[.5 6.5],'YDir','reverse','XTick',1:6,'YTick',1:6, ...
        'XTickLabel',C.labels,'YTickLabel',C.labels,'TickLabelInterpreter','tex', ...
        'TickLength',[0 0],'FontName','Times New Roman','FontSize',opts.FontSize, ...
        'LineWidth',.8,'Box','on','DataAspectRatio',[1 1 1]);
    title(ax,titles{k},'FontSize',opts.FontSize+1,'FontWeight','bold');
end
% Vector color scale; MATLAB colorbar/imagesc can embed a raster in the PDF.
ax=axes('Parent',f,'Position',[.30 .095 .40 .035]);hold(ax,'on');edges=linspace(-1,1,129);
for j=1:numel(edges)-1
    rectangle(ax,'Position',[edges(j),0,edges(j+1)-edges(j),1], ...
        'FaceColor',cr25_color(mean(edges(j:j+1))),'EdgeColor','none');
end
set(ax,'XLim',[-1 1],'YLim',[0 1],'XTick',[-1 -.5 0 .5 1],'YTick',[], ...
    'FontName','Times New Roman','FontSize',opts.FontSize-1,'TickDir','out','Box','on','LineWidth',.7);
cr25_save(f,out,'Fig07_correlation_matrices',opts);

f=cr25_figure('Correlation with cavitation',4.6,opts);
ax=axes('Parent',f,'Position',[.20 .18 .73 .65]);hold(ax,'on');
colors=[0 .447 .698;.835 .369 0];labels={'W','H','k','\sigma_c','V'};
line(ax,[0 0],[.5 5.5],'Color',[.52 .55 .58],'LineStyle','--','LineWidth',1);
h=gobjects(2,1);
for k=1:2
    if k==1,p=C.summary.Pearson;lo=C.summary.PearsonLower;hi=C.summary.PearsonUpper;offset=.12;marker='o';
    else,p=C.summary.Spearman;lo=C.summary.SpearmanLower;hi=C.summary.SpearmanUpper;offset=-.12;marker='s';end
    for j=1:5
        yy=6-j+offset;
        if isfinite(lo(j))&&isfinite(hi(j))
            line(ax,[lo(j) hi(j)],[yy yy],'Color',colors(k,:),'LineWidth',1.8);
            line(ax,[lo(j) lo(j)],yy+[-.045 .045],'Color',colors(k,:),'LineWidth',1.2);
            line(ax,[hi(j) hi(j)],yy+[-.045 .045],'Color',colors(k,:),'LineWidth',1.2);
        end
        mark=plot(ax,p(j),yy,marker,'Color',colors(k,:),'MarkerFaceColor',colors(k,:), ...
            'MarkerEdgeColor','w','MarkerSize',8,'LineWidth',1.1);
        if j==1,h(k)=mark;end
    end
end
set(ax,'XLim',[-1 1],'XTick',-1:.25:1,'YLim',[.5 5.5],'YTick',1:5, ...
    'YTickLabel',fliplr(labels),'TickLabelInterpreter','tex','FontName','Times New Roman', ...
    'FontSize',opts.FontSize,'TickDir','out','Box','off','LineWidth',1);
grid(ax,'on');ax.YGrid='off';ax.GridAlpha=.3;
xlabel(ax,'Correlation coefficient','FontSize',opts.FontSize+1);
title(ax,'Correlation with cavitation','FontSize',opts.FontSize+1,'FontWeight','bold');
lg=legend(ax,h,{'Pearson','Spearman'},'Orientation','horizontal','Box','off','FontSize',opts.FontSize-1);
lg.Units='normalized';lg.Position=[.30 .90 .50 .055];
cr25_save(f,out,'Fig08_correlations_with_cavitation',opts);
end

function M=cr25_corr(X)
Z=X-mean(X,1);ss=sum(Z.^2,1);den=sqrt(ss'*ss);M=(Z'*Z)./den;
M(den==0)=NaN;M(isfinite(M))=max(-1,min(1,M(isfinite(M))));
end

function R=cr25_ranks(X)
R=zeros(size(X));n=size(X,1);
for j=1:size(X,2)
    [values,order]=sort(X(:,j));r=zeros(n,1);a=1;
    while a<=n
        b=a;while b<n&&values(b+1)==values(a),b=b+1;end
        r(a:b)=(a+b)/2;a=b+1;
    end
    R(order,j)=r;
end
end

function q=cr25_quantile(x,p)
x=sort(x(:));q=NaN(size(p));if isempty(x),return;end
h=1+(numel(x)-1)*p;lo=floor(h);hi=ceil(h);
q=reshape(x(lo),size(p)).*(hi-h)+reshape(x(hi),size(p)).*(h-lo);
same=lo==hi;q(same)=x(lo(same));
end

function color=cr25_color(v)
if ~isfinite(v),color=[.88 .89 .90];return;end
white=[.985 .985 .985];negative=[.10 .34 .65];positive=[.76 .16 .20];
t=min(1,abs(v));if v<0,edge=negative;else,edge=positive;end
color=(1-t)*white+t*edge;
end

function f=cr25_figure(name,h,opts)
visibility='off';if opts.ShowFigures,visibility='on';end
ss=get(groot,'ScreenSize');if ss(3)<100,ss=[1 1 1600 1000];end
w=min(1450,.92*ss(3));hh=min(860,.85*ss(4));
f=figure('Name',name,'NumberTitle','off','Color','w','Units','pixels', ...
    'Position',[max(10,(ss(3)-w)/2),max(20,(ss(4)-hh)/2),w,hh], ...
    'Visible',visibility,'Renderer','painters','DefaultAxesFontName','Times New Roman', ...
    'DefaultTextFontName','Times New Roman','DefaultAxesFontSize',opts.FontSize,'DefaultTextFontSize',opts.FontSize);
set(f,'PaperUnits','inches','PaperSize',[opts.FigureWidth h], ...
    'PaperPosition',[0 0 opts.FigureWidth h]);
end

function cr25_save(f,out,name,opts)
drawnow;print(f,fullfile(out,[name '.pdf']),'-dpdf','-painters');
print(f,fullfile(out,[name '.png']),'-dpng',sprintf('-r%d',opts.ExportDPI));
savefig(f,fullfile(out,[name '.fig']));
if ~opts.ShowFigures,close(f);end
end

function cav25_publication_figures(A,R,V,T,out,opts) %#ok<INUSD>
% Publication figures; keep windows open unless ShowFigures=false.
if ~isfield(opts,'ShowFigures'),opts.ShowFigures=true;end
if ~isfield(opts,'FigureWidth'),opts.FigureWidth=10;end
if ~isfield(opts,'ExportDPI'),opts.ExportDPI=900;end
if ~isfield(opts,'FontSize'),opts.FontSize=14;end
assert(opts.FigureWidth>0 && opts.ExportDPI>=72 && opts.FontSize>0,'Invalid figure size, DPI or font size.');
if ~isfolder(out),mkdir(out);end
C=[0 .447 .698;.835 .369 0;0 .620 .451;.494 .184 .556];
pp25_data(A,out,opts,C);
pp25_forest(R,out,opts,C);
% Complete-model diagnostics are rendered by cav25_auc_overview.
pp25_confusions(A,V,out,opts,C);
pp25_conditional(A,V,out,opts,C);
pp25_repeats(V,out,opts,C);
fprintf('Five publication figures exported to %s\n',out);
if opts.ShowFigures
    fprintf('Figure windows are kept open. Vector PDFs, %d dpi PNGs and editable FIG files are saved.\n',opts.ExportDPI);
end
end

function f=pp25_figure(name,h,opts)
visibility='off';if opts.ShowFigures,visibility='on';end
w=opts.FigureWidth;paperH=h*w/7.2;
% Enlarge the screen window independently of export resolution; keep its ratio.
rootUnits=get(groot,'Units');set(groot,'Units','pixels');
screen=get(groot,'ScreenSize');set(groot,'Units',rootUnits);
pw=min(1400,max(480,screen(3)-100));ph=pw*paperH/w;
if ph>screen(4)-110,ph=max(360,screen(4)-110);pw=ph*w/paperH;end
left=max(20,(screen(3)-pw)/2);bottom=max(45,(screen(4)-ph)/2);
f=figure('Name',name,'NumberTitle','off','Color','w','Units','pixels', ...
    'Position',[left bottom pw ph],'Visible',visibility,'Renderer','painters', ...
    'DefaultAxesFontName','Times New Roman','DefaultTextFontName','Times New Roman', ...
    'DefaultAxesFontSize',opts.FontSize,'DefaultTextFontSize',opts.FontSize);
setappdata(f,'PublicationFontSize',opts.FontSize);
set(f,'PaperUnits','inches','PaperSize',[w paperH],'PaperPosition',[0 0 w paperH],'PaperPositionMode','manual');
end

function ax=pp25_axes(f,pos)
fs=getappdata(f,'PublicationFontSize');
ax=axes('Parent',f,'Position',pos,'FontName','Times New Roman','FontSize',fs, ...
    'LineWidth',1.0,'Box','off','TickDir','out','Layer','top');
hold(ax,'on');grid(ax,'on');ax.GridColor=[.82 .84 .86];ax.GridAlpha=.45;
end

function pp25_title(ax,s)
title(ax,s,'FontSize',ax.FontSize+1,'FontWeight','bold','HorizontalAlignment','left');
ax.Title.Units='normalized';ax.Title.Position(1)=0;
end

function pp25_save(f,out,name,opts)
% Log-axis minor grids can overwhelm the data; retain only major grids.
axs=findall(f,'Type','axes');
for j=1:numel(axs),set(axs(j),'XMinorGrid','off','YMinorGrid','off');end
drawnow;
print(f,fullfile(out,[name '.pdf']),'-dpdf','-painters');
print(f,fullfile(out,[name '.png']),'-dpng',sprintf('-r%d',round(opts.ExportDPI)));
savefig(f,fullfile(out,[name '.fig']));
if ~opts.ShowFigures,close(f);end
end

function m=pp25_model(V,name)
j=find(strcmp({V.models.name},name),1);
assert(~isempty(j),'Requested plot model %s is absent.',name);
m=V.models(j);
end

function s=pp25_name(name)
switch name
    case 'G',s='Geometry';
    case 'GK',s='Geometry + TKE';
    case 'Full_raw',s='Full (\sigma_c)';
    case 'Full_log',s='Full (ln \sigma_c)';
    otherwise,s=strrep(name,'_',' ');
end
end

function pp25_data(A,out,opts,C)
f=pp25_figure('Data and observed labels',5.5,opts);
y=A(:,6);v=A(:,5);xx={1000*A(:,1),1000*A(:,2),A(:,3),A(:,4)};
yl={'Width W (mm)','Height H (mm)','k (m^2 s^{-2})','Cavitation number \sigma_c'};
titles={'(a) Roughness width','(b) Roughness height','(c) Turbulent energy','(d) Cavitation number'};
positions=[.085 .59 .235 .33;.415 .59 .235 .33;.745 .59 .235 .33; ...
           .085 .10 .235 .33;.415 .10 .235 .33;.745 .10 .235 .33];
for j=1:4
    ax=pp25_axes(f,positions(j,:));x=xx{j};
    for g=0:1
        values=x(y==g);n=numel(values);jit=.14*sin((1:n)'*2.39996323);
        plot(ax,g+1+jit,values,'.','Color',.65*C(g+1,:)+.35,'MarkerSize',8);
        qu=pp25_quantile(values,[.25 .5 .75]);iq=qu(3)-qu(1);
        low=min(values(values>=qu(1)-1.5*iq));high=max(values(values<=qu(3)+1.5*iq));
        line(ax,[g+1 g+1],[low high],'Color',C(g+1,:),'LineWidth',1.1);
        line(ax,g+1+[-.1 .1],[low low],'Color',C(g+1,:),'LineWidth',1.1);
        line(ax,g+1+[-.1 .1],[high high],'Color',C(g+1,:),'LineWidth',1.1);
        patch(ax,g+1+[-.19 .19 .19 -.19],[qu(1) qu(1) qu(3) qu(3)], ...
            .86+.14*C(g+1,:),'EdgeColor',C(g+1,:),'LineWidth',1.1);
        line(ax,g+1+[-.19 .19],[qu(2) qu(2)],'Color',C(g+1,:),'LineWidth',1.6);
    end
    if j==4 && all(x>0),set(ax,'YScale','log');end
    xlim(ax,[.55 2.45]);set(ax,'XTick',[1 2],'XTickLabel',{'Y = 0','Y = 1'});
    ylabel(ax,yl{j},'FontSize',opts.FontSize);pp25_title(ax,titles{j});
    if j==3 && max(x)>8
        % The full linear axis retains every point, including the maximum.
        % An explicitly bounded inset resolves the dense low-k region.
        inset=axes('Parent',f,'Position',[.835 .72 .13 .14], ...
            'Color','w','Box','on','FontName','Times New Roman', ...
            'FontSize',opts.FontSize-3,'LineWidth',.8,'TickDir','out');
        hold(inset,'on');
        for g=0:1
            values=x(y==g);jit=.11*sin((1:numel(values))'*2.39996323);
            plot(inset,g+1+jit,values,'.','Color',C(g+1,:),'MarkerSize',4);
            med=pp25_quantile(values,.5);
            line(inset,g+1+[-.2 .2],[med med],'Color',C(g+1,:),'LineWidth',1.5);
        end
        set(inset,'XTick',[1 2],'XTickLabel',{'0','1'},'YTick',[0 2 4]);
        xlim(inset,[.6 2.4]);ylim(inset,[0 4]);
        title(inset,'k \leq 4','FontSize',opts.FontSize-2,'FontWeight','normal');
    end
end
ax=pp25_axes(f,positions(5,:));counts=[sum(y==0),sum(y==1)];
for j=1:2
    rectangle(ax,'Position',[j-.29 0 .58 counts(j)],'FaceColor',C(j,:),'EdgeColor','none');
end
set(ax,'XTick',[1 2],'XTickLabel',{'Y = 0','Y = 1'});xlim(ax,[.45 2.55]);
ylim(ax,[0 max(counts)*1.23]);ylabel(ax,'Number of records');pp25_title(ax,'(e) Cavitation labels');
for j=1:2,text(ax,j,counts(j)+max(counts)*.035,num2str(counts(j)), ...
        'HorizontalAlignment','center','FontSize',opts.FontSize);end
ax=pp25_axes(f,positions(6,:));
for g=0:1
    n=sum(v==g);s=sum(y(v==g));pr=s/n;ci=pp25_wilson(s,n);
    errorbar(ax,g+1,pr,pr-ci(1),ci(2)-pr,'o','Color',C(g+1,:), ...
        'MarkerFaceColor',C(g+1,:),'MarkerEdgeColor','w','MarkerSize',9,'LineWidth',1.3,'CapSize',8);
    text(ax,g+1,min(1.09,ci(2)+.06),sprintf('%d/%d',s,n),'HorizontalAlignment','center','Color',C(g+1,:));
end
set(ax,'XTick',[1 2],'XTickLabel',{'V = 0','V = 1'});xlim(ax,[.45 2.55]);ylim(ax,[0 1.16]);
set(ax,'YTick',0:.25:1);ylabel(ax,'Observed cavitation fraction');pp25_title(ax,'(f) Vortex groups');
pp25_save(f,out,'Fig01_data_and_labels',opts);
end

function pp25_forest(R,out,opts,C)
f=pp25_figure('Adjusted associations: unpenalized full model',3.55,opts);
ax=pp25_axes(f,[.26 .20 .68 .67]);
if ~isfield(R,'models') || ~any(strcmp({R.models.name},'Full_raw'))
    axis(ax,'off');text(ax,.5,.5,'Finite unpenalized full-model fit unavailable.', ...
        'Units','normalized','HorizontalAlignment','center');
else
    m=R.models(strcmp({R.models.name},'Full_raw'));ix=2:6;yy=5:-1:1;
    set(ax,'XScale','log');line(ax,[1 1],[.4 5.6],'Color',[.55 .55 .55],'LineStyle','--');
    for j=1:5
        line(ax,m.OR95(ix(j),:),[yy(j) yy(j)],'LineWidth',1.5,'Color',C(1,:));
        plot(ax,m.OR(ix(j)),yy(j),'o','MarkerFaceColor',C(1,:), ...
            'MarkerEdgeColor','w','MarkerSize',9);
    end
    set(ax,'YTick',1:5,'YTickLabel',{'V: 0 to 1','\sigma_c','k (m^2 s^{-2})','Height H','Width W'});
    xlim(ax,[10^floor(log10(min(m.OR95(ix,1)))) 10^ceil(log10(max(m.OR95(ix,2))))]);
    ylim(ax,[.4 5.6]);xlabel(ax,'Adjusted odds ratio (log axis)','FontSize',opts.FontSize+1);
    pp25_title(ax,'Adjusted odds ratios');
end
pp25_save(f,out,'Fig02_adjusted_odds_ratios',opts);
end

function pp25_diagnostics(A,V,out,opts,C)
f=pp25_figure('ROC and calibration',4.2,opts);
ax1=pp25_axes(f,[.095 .18 .36 .63]);ax2=pp25_axes(f,[.60 .18 .36 .63]);
line(ax1,[0 1],[0 1],'Color',[.6 .63 .65],'LineStyle','--','HandleVisibility','off');
line(ax2,[0 1],[0 1],'Color',[.6 .63 .65],'LineStyle','--','HandleVisibility','off');
want={'GK','Full_raw','Full_log'};want=want(ismember(want,{V.models.name}));
h=gobjects(numel(want),1);labs=cell(numel(want),1);markers={'o','^','d'};
for j=1:numel(want)
    m=pp25_model(V,want{j});p=m.oof(:,1);roc=pp25_roc(A(:,6),p);
    h(j)=stairs(ax1,roc(:,1),roc(:,2),'Color',C(j,:),'LineWidth',2.0);
    labs{j}=sprintf('%s   AUC %.3f',pp25_name(want{j}),pp25_auc(A(:,6),p));
    cal=pp25_calibration(A(:,6),p,5);
    plot(ax2,cal(:,1),cal(:,2),markers{j},'Color',C(j,:),'MarkerFaceColor',C(j,:), ...
        'MarkerEdgeColor','w','MarkerSize',10,'LineStyle','none');
end
for ax=[ax1 ax2]
    xlim(ax,[0 1]);ylim(ax,[0 1]);set(ax,'XTick',0:.2:1,'YTick',0:.2:1);
end
xlabel(ax1,'False-positive fraction');ylabel(ax1,'True-positive fraction');pp25_title(ax1,'(a) ROC curves');
xlabel(ax2,'Predicted probability');ylabel(ax2,'Observed cavitation fraction');pp25_title(ax2,'(b) Calibration');
lg=legend(ax1,h,labs,'Orientation','horizontal','Box','off','FontSize',opts.FontSize-1);
lg.Units='normalized';lg.Position=[.075 .91 .90 .055];
pp25_save(f,out,'Fig03_nestedCV_ROC_calibration',opts);
end

function pp25_confusions(A,V,out,opts,C)
f=pp25_figure('Classification counts',3.85,opts);
m=pp25_model(V,'Full_raw');preds={m.oof(:,1)>=.5,m.predictedLabel(:,1)};
titles={'(a) Fixed threshold: 0.5','(b) Inner-fold threshold'};
pos=[.12 .21 .33 .63;.62 .21 .33 .63];
for j=1:2
    ax=axes('Parent',f,'Position',pos(j,:));y=A(:,6);p=preds{j};
    cm=[sum(y==0&~p),sum(y==0&p);sum(y==1&~p),sum(y==1&p)];
    % Draw vector rectangles so the PDF contains no raster heatmap.
    hold(ax,'on');scale=max(sum(y==0),sum(y==1));
    for rr=1:2
        for cc=1:2
            fraction=cm(rr,cc)/scale;fillColor=(1-fraction)*[1 1 1]+fraction*C(1,:);
            rectangle(ax,'Position',[cc-.5 rr-.5 1 1],'FaceColor',fillColor,'EdgeColor','none');
        end
    end
    xlim(ax,[.5 2.5]);ylim(ax,[.5 2.5]);
    set(ax,'YDir','reverse','XTick',[1 2],'YTick',[1 2], ...
        'XTickLabel',{'Y = 0','Y = 1'},'YTickLabel',{'Y = 0','Y = 1'}, ...
        'FontName','Times New Roman','FontSize',opts.FontSize-1,'TickLength',[0 0],'Box','on','LineWidth',.8);
    for r=1:2
        for c=1:2
            color=[.1 .12 .15];if cm(r,c)>max(sum(y==0),sum(y==1))*.58,color=[1 1 1];end
            text(ax,c,r-.08,sprintf('%d',cm(r,c)),'HorizontalAlignment','center', ...
                'FontSize',opts.FontSize+12,'FontWeight','bold','Color',color);
            text(ax,c,r+.19,sprintf('%.1f%%',100*cm(r,c)/sum(cm(r,:))), ...
                'HorizontalAlignment','center','FontSize',opts.FontSize,'Color',color);
        end
    end
    xlabel(ax,'Predicted label','FontSize',opts.FontSize);ylabel(ax,'Observed label','FontSize',opts.FontSize);
    pp25_title(ax,titles{j});ax.Title.FontSize=opts.FontSize+1;
end
pp25_save(f,out,'Fig04_nestedCV_confusion_matrices',opts);
end

function pp25_conditional(A,V,out,opts,C)
f=pp25_figure('Conditional cavitation probability',4.15,opts);
ax=pp25_axes(f,[.105 .20 .63 .62]);m=pp25_model(V,'Full_raw');fit=m.final;
lower=max(min(A(A(:,5)==0,4)),min(A(A(:,5)==1,4)));
upper=min(max(A(A(:,5)==0,4)),max(A(A(:,5)==1,4)));
if lower>=upper
    axis(ax,'off');text(ax,.5,.5,'The two vortex groups have no common observed sigma range.', ...
        'Units','normalized','HorizontalAlignment','center');
else
    if lower>0
        sg=logspace(log10(lower),log10(upper),301)';set(ax,'XScale','log');
    else
        sg=linspace(lower,upper,301)';
    end
    hh=gobjects(2,1);crossings=NaN(2,1);curves=zeros(numel(sg),2);
    for vv=0:1
        eta=fit.beta(1)+fit.beta(5)*(sg-fit.mu(4))/fit.sd(4)+fit.beta(6)*vv;
        curves(:,vv+1)=pp25_sigmoid(eta);
        hh(vv+1)=plot(ax,sg,curves(:,vv+1),'Color',C(vv+1,:),'LineWidth',2.1);
        if abs(fit.beta(5))>1e-12 && fit.threshold>0 && fit.threshold<1
            crossing=fit.mu(4)+fit.sd(4)*(log(fit.threshold/(1-fit.threshold))-fit.beta(1)-fit.beta(6)*vv)/fit.beta(5);
            if crossing>=lower && crossing<=upper
                crossings(vv+1)=crossing;
                line(ax,[crossing crossing],[0 fit.threshold],'Color',C(vv+1,:), ...
                    'LineStyle',':','LineWidth',1.1,'HandleVisibility','off');
                plot(ax,crossing,fit.threshold,'o','Color',C(vv+1,:), ...
                    'MarkerFaceColor',C(vv+1,:),'MarkerEdgeColor','w','MarkerSize',8,'HandleVisibility','off');
            end
        end
    end
    line(ax,[lower upper],[fit.threshold fit.threshold],'Color',[.25 .25 .25], ...
        'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    text(ax,1.035,.93,sprintf('p_* = %.3f',fit.threshold),'Units','normalized','FontSize',opts.FontSize+1);
    for vv=0:1
        if isfinite(crossings(vv+1)),s=sprintf('V = %d\n\\sigma_c^* = %.3g',vv,crossings(vv+1));
        else,s=sprintf('V = %d\nNo crossing\nin common range',vv);end
        text(ax,1.035,.72-.27*vv,s,'Units','normalized','FontSize',opts.FontSize,'Color',C(vv+1,:), ...
            'VerticalAlignment','top','Interpreter','tex');
    end
    xlim(ax,[lower upper]);ylim(ax,[0 1]);set(ax,'YTick',0:.2:1);
    xlabel(ax,'Cavitation number \sigma_c');ylabel(ax,'Fitted conditional probability');
    pp25_title(ax,'Conditional cavitation probability');
    lg=legend(ax,hh,{'V = 0','V = 1'},'Orientation','horizontal','Box','off');
    lg.Units='normalized';lg.Position=[.22 .91 .43 .055];
    tab=array2table([sg,curves],'VariableNames',{'Sigma','Probability_V0','Probability_V1'});
    writetable(tab,fullfile(out,'figure_conditional_probability_curves.csv'));
end
pp25_save(f,out,'Fig05_conditional_probability_boundary',opts);
end

function pp25_repeats(V,out,opts,C)
f=pp25_figure('Cross-validation performance',4.15,opts);
pos=[.095 .25 .36 .61;.60 .25 .36 .61];names={V.models.name};n=numel(names);
labels=cellfun(@pp25_name,names,'UniformOutput',false);metrics={'AUC','Brier'};
for k=1:2
    ax=pp25_axes(f,pos(k,:));vals=[];
    for j=1:n,vals(:,j)=V.models(j).repeatMetrics.(metrics{k});end %#ok<AGROW>
    nr=size(vals,1);
    for r=1:nr,plot(ax,1:n,vals(r,:),'-','Color',[.84 .86 .88],'LineWidth',.6,'HandleVisibility','off');end
    for j=1:n
        jit=.11*sin((1:nr)'*2.39996323);cc=C(1+mod(j-1,size(C,1)),:);
        plot(ax,j+jit,vals(:,j),'o','Color',cc,'MarkerFaceColor',cc,'MarkerEdgeColor','w','MarkerSize',7.5);
        med=pp25_quantile(vals(:,j),.5);line(ax,j+[-.20 .20],[med med],'Color',[.15 .15 .15],'LineWidth',2.1);
    end
    set(ax,'XTick',1:n,'XTickLabel',labels,'XTickLabelRotation',22,'FontSize',opts.FontSize-1);
    xlim(ax,[.6 n+.4]);ylabel(ax,metrics{k});
    if k==1
        pp25_title(ax,'(a) AUC');
        yl=ylim(ax);ylim(ax,[max(0,min(.45,yl(1))) min(1,max(.95,yl(2)))]);
        line(ax,[.6 n+.4],[.5 .5],'Color',[.55 .55 .55],'LineStyle','--');
    else
        pp25_title(ax,'(b) Brier score');
        if isfield(V,'baseline') && isfield(V.baseline,'repeatMetrics')
            base=median(V.baseline.repeatMetrics.Brier);
            hbase=line(ax,[.6 n+.4],[base base],'Color',[.45 .45 .45],'LineStyle','--');
            legend(ax,hbase,{'Prevalence baseline'},'Location','northeast','Box','off','FontSize',opts.FontSize-2);
        end
    end
end
pp25_save(f,out,'Fig06_nestedCV_partition_sensitivity',opts);
end

function p=pp25_sigmoid(e)
p=exp(-max(0,-e)-log1p(exp(-abs(e))));
end

function a=pp25_auc(y,p)
pos=p(y==1);neg=p(y==0)';a=mean((pos>neg)+.5*(pos==neg),'all');
end

function roc=pp25_roc(y,p)
thresholds=sort(unique(p),'descend');roc=zeros(numel(thresholds)+1,2);
for j=1:numel(thresholds)
    pred=p>=thresholds(j);roc(j+1,:)=[sum(pred&y==0)/sum(y==0),sum(pred&y==1)/sum(y==1)];
end
end

function C=pp25_calibration(y,p,g)
[~,order]=sortrows([p,(1:numel(y))'],[1 2]);sizes=floor(numel(y)/g)*ones(g,1);
sizes(1:mod(numel(y),g))=sizes(1:mod(numel(y),g))+1;start=0;C=zeros(g,4);
for j=1:g
    ix=order(start+(1:sizes(j)));start=start+sizes(j);
    C(j,:)=[mean(p(ix)),mean(y(ix)),numel(ix),sum(y(ix))];
end
end

function q=pp25_quantile(x,p)
x=sort(x(:));q=NaN(size(p));if isempty(x),return;end
h=1+(numel(x)-1)*p;lo=floor(h);hi=ceil(h);
q=reshape(x(lo),size(p)).*(hi-h)+reshape(x(hi),size(p)).*(h-lo);
same=lo==hi;q(same)=x(lo(same));
end

function ci=pp25_wilson(s,n)
z=1.959963984540054;p=s/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;
r=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d;ci=[c-r,c+r];
end

function cav25_write_report(Result,out)
V=Result.validation;T=Result.thresholds;R=Result.association;o=Result.options;
fid=fopen(fullfile(out,'READ_RESULTS.txt'),'w','n','UTF-8');assert(fid>0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Cavitation-label logistic analysis: actual computed results\n');
fprintf(fid,'N=%d; positives=%d; negatives=%d\n\n',R.n,R.positive,R.negative);
fprintf(fid,'All regressions use k directly (linear standardization only), with no log transform of k.\n\n');
fprintf(fid,'1. Raw-k unpenalized maximum-likelihood analysis\n');
for m=1:numel(R.models)
    a=R.models(m);fprintf(fid,'%s: apparent AUC %.8f; Brier %.8f; TP/FN/TN/FP %d/%d/%d/%d at p=0.5\n',a.name,a.auc,a.brier,a.cm(1:4));
end
fprintf(fid,'The raw-k full-model coefficients and Wald intervals are in association/coefficients.csv.\n');
fprintf(fid,'Those coefficients are not the penalized coefficients of the models below.\n\n');
fprintf(fid,'2. Ridge models: %d repeats of stratified %d-fold outer / %d-fold inner row CV\n',o.CVRepeats,o.OuterFolds,o.InnerFolds);
fprintf(fid,'Training-fold scaling; inner pooled log loss selects lambda; inner OOF Youden selects probability cutoff.\n');
fprintf(fid,'Objective = mean logistic negative log likelihood + lambda/2 * sum(nonintercept beta^2).\n');
fprintf(fid,'Outer test labels never select lambda or cutoff. Continuous predictors are standardized; V stays 0/1.\n');
fprintf(fid,'Each repeat contains one outer-fold prediction per original record. Repeat variability is not a confidence interval.\n\n');
for m=1:numel(V.models)
    f=V.models(m).final;
    fprintf(fid,'%s final lambda=%.12g; final internally tuned probability cutoff=%.12g\n',V.models(m).name,f.lambda,f.threshold);
    fprintf(fid,'Feature order: %s\n',strjoin(f.featureNames,', '));
    fprintf(fid,'Final beta including intercept: %s\n',mat2str(f.beta',10));
    fprintf(fid,'Scaling means: %s\nScaling SDs: %s\n',mat2str(f.mu,10),mat2str(f.sd,10));
    fprintf(fid,'%s\n',evalc('disp(V.models(m).repeatMetrics)'));
    if f.lambda==min(o.Lambdas)||f.lambda==max(o.Lambdas)
        fprintf(fid,'NOTE: selected lambda is at the prescribed grid boundary; inspect sensitivity to extending the grid before making stability claims.\n');
    end
end
fprintf(fid,'3. Multivariable conditional boundaries\n');
fprintf(fid,'For a target probability t, solve beta0 + sum(beta_j*z_j) = log(t/(1-t)).\n');
fprintf(fid,'For raw sigma: sigma_t = mu_sigma + sd_sigma*(logit(t)-beta0-betaW*zW-betaH*zH-betaK*zk-betaV*V)/betaSigma.\n');
fprintf(fid,'If betaSigma<0, classify Y=1 at sigma<=sigma_t; reverse if positive.\n');
fprintf(fid,'Conditional reference W=%.9g m, H=%.9g m, k=%.9g m^2/s^2 (arithmetic sample means).\n',T.referenceRaw);
fprintf(fid,'%s\n',evalc('disp(T.referenceCrossings)'));
fprintf(fid,'%s\n',evalc('disp(T.bootstrapSummary)'));
fprintf(fid,'The %d threshold bootstraps repeat scaling, lambda selection and probability-threshold selection.\n',o.ThresholdBootstrap);
fprintf(fid,'Copies of the same resampled record stay in one inner fold. This copy grouping is NOT grouping by surface/time/peak.\n');
fprintf(fid,'Percentile ranges condition on successful finite crossings and assume row independence. Outside-range crossings are flagged, not clipped.\n');
fprintf(fid,'A crossing inside marginal ranges does not establish joint covariate support or physical inception.\n');
fprintf(fid,'Univariate Youden cutpoints in association/univariate_cutpoints.csv are separate apparent descriptive cutpoints.\n\n');
fprintf(fid,'4. Interpretation\n');
fprintf(fid,'The model classifies contemporaneous manually labelled monitoring records. Local k, sigma and V are not guaranteed pre-cavitation information.\n');
fprintf(fid,'Surface/time/protrusion IDs are unavailable. These row-wise diagnostics do not satisfy grouped validation and cannot establish transport to new surfaces.\n');
fprintf(fid,'No conclusion of reliability is forced from AUC. Inspect Brier/log loss against the outer-training prevalence baseline, calibration, false negatives and split/threshold stability together.\n');
fprintf(fid,'5. Correlation analysis\nPearson and Spearman matrices and correlations with Y are in correlation/. Continuous-Y Pearson values are point-biserial; binary V-Y is phi. Bootstrap ranges use working row independence.\n');
fprintf(fid,'6. Figures\nAll Fig*.pdf, Fig*.png and Fig*.fig files are at this report location. ShowFigures=true retains all ten windows.\n');
fprintf(fid,'Reopen without refitting: cavitation25_open_figures(Result).\n');
end

function U=cav25_refinement(A,V,opts,out)
% Fixed, exploratory AUC-prioritized refinements, evaluated on existing splits.
% No outer outcome is used to choose the family, penalty, or threshold.
if ~isfield(opts,'Lambdas'),opts.Lambdas=V.options.Lambdas;end
lambdas=unique(opts.Lambdas(:)','sorted');
assert(all(isfinite(lambdas))&&all(lambdas>0),'Positive finite penalties required.');
y=A(:,6);n=numel(y);R=size(V.folds,2);K=size(V.innerFolds,2);
assert(size(A,2)==6&&isequal(size(V.folds),[n R]));
U=struct('scope','Exploratory nested ROW-level internal assessment; no grouped or external validation.', ...
    'selectionRule','Inner OOF AUC descending, log loss ascending, additive before interaction, lambda descending; ties within 1e-12.', ...
    'features','W_mm, H_mm, raw k, sigma_c, V; candidate interaction sigma_c*V. Training-only scaling; V unchanged.', ...
    'folds',V.folds,'innerFolds',V.innerFolds,'lambdas',lambdas);
names={'Full_AUC','Coupling_AUC'};
blank=struct('name','','oof',NaN(n,R),'oofEta',NaN(n,R), ...
    'predictedLabel',false(n,R),'thresholdByRecord',NaN(n,R), ...
    'repeatMetrics',table(),'final',struct());
U.models=repmat(blank,1,2);
for m=1:2,U.models(m).name=names{m};end
foldRows={};candidateRows={};
for r=1:R
    for k=1:K
        test=V.folds(:,r)==k;train=~test;
        assert(all(V.innerFolds(test,k,r)==0),'Outer records found in inner tuning folds.');
        inner=V.innerFolds(train,k,r);assert(all(inner>0));
        for m=1:2
            T=rf25_tune(A(train,:),inner,lambdas,m==2);
            [F,s,fn]=rf25_features(A,T.family);
            f=cv25_fit(F(train,:),y(train),T.lambda,s);
            [p,e]=cv25_predict(f,F(test,:));
            U.models(m).oof(test,r)=p;U.models(m).oofEta(test,r)=e;
            U.models(m).thresholdByRecord(test,r)=T.threshold;
            U.models(m).predictedLabel(test,r)=p>=T.threshold;
            foldRows(end+1,:)={names{m},r,k,sum(train),sum(test),T.family,T.lambda,T.threshold,T.auc,T.logloss,T.failedFits,f.gradientNorm}; %#ok<AGROW>
            candidateRows=[candidateRows;rf25_candidate_rows(T,names{m},r,k)]; %#ok<AGROW>
        end
    end
    for m=1:2
        z=U.models(m);
        U.models(m).repeatMetrics=[z.repeatMetrics;cv25_metrics(y,z.oof(:,r),z.predictedLabel(:,r),z.oofEta(:,r),r)];
    end
    fprintf('Exploratory refinement: repeat %d/%d completed.\n',r,R);
end
U.finalInnerFolds=V.finalInnerFolds;
for m=1:2
    T=rf25_tune(A,V.finalInnerFolds,lambdas,m==2);
    [F,s,fn]=rf25_features(A,T.family);f=cv25_fit(F,y,T.lambda,s);
    [f.p,f.eta]=cv25_predict(f,F);f.family=T.family;
    f.threshold=T.threshold;f.featureNames=fn;f.innerOOF=T.p;f.innerEta=T.eta;
    f.tuningTable=T.table;f.selectedCandidate=T.index;f.innerAUC=T.auc;f.innerLogLoss=T.logloss;
    f.betaRaw=[f.beta(1)-sum((f.mu./f.sd).*f.beta(2:end)');f.beta(2:end)./f.sd'];
    f.betaInputUnits=f.betaRaw;f.betaInputUnits(2:3)=1000*f.betaRaw(2:3);
    U.models(m).final=f;
    candidateRows=[candidateRows;rf25_candidate_rows(T,names{m},0,0)]; %#ok<AGROW>
end
U.foldSummary=cell2table(foldRows,'VariableNames',{'Model','Repeat','OuterFold','TrainN','TestN','Family','Lambda','ProbabilityThreshold','InnerAUC','InnerLogLoss','FailedInnerFits','GradientNorm'});
U.candidateMetrics=cell2table(candidateRows,'VariableNames',{'Model','Repeat','OuterFold','Candidate','Family','Lambda','InnerAUC','InnerLogLoss','FailedFits','FailureStatus','Selected'});
base=V.models(strcmp({V.models.name},'Full_raw'));
ss=rf25_summary_row(base.repeatMetrics,'Full_logloss');
for m=1:2,ss=[ss;rf25_summary_row(U.models(m).repeatMetrics,names{m})];end %#ok<AGROW>
U.summary=cell2table(ss,'VariableNames',{'Model','MeanAUC','SDAUC','MinAUC','MaxAUC','MeanBrier','SDBrier','MeanLogLoss','MeanSensitivity','MeanSpecificity','MeanAccuracy','MeanCalibrationSlope','MeanCalibrationIntercept'});
U.totalFailedInnerFits=sum(U.candidateMetrics.FailedFits);
U.outerFitFailures=0; % Any outer fitting failure aborts rather than deleting a repeat.
if ~isempty(out),rf25_export(U,A,out);end
end

function [F,s,names]=rf25_features(A,family)
F=[1000*A(:,1:2),A(:,3:5)];
names={'W_mm','H_mm','k_m2s2','sigma_c','V'};s=[true true true true false];
if strcmp(family,'interaction')
    F=[F,A(:,4).*A(:,5)];s=[s,true];names{end+1}='sigma_c_x_V';
else
    assert(strcmp(family,'additive'),'Unknown family.');
end
end

function T=rf25_tune(A,inner,lambdas,allowInteraction)
y=A(:,6);families={'additive'};if allowInteraction,families{2}='interaction';end
folds=unique(inner(:))';assert(numel(folds)>=2&&all(folds>0));
n=numel(y);L=numel(lambdas);nc=L*numel(families);
P=NaN(n,nc);E=P;row=cell(nc,6);failed=zeros(nc,1);auc=-Inf(nc,1);ll=Inf(nc,1);
for g=1:numel(families)
    [F,s]=rf25_features(A,families{g});
    for a=1:L
        ix=(g-1)*L+a;status={};
        for k=folds
            te=inner==k;tr=~te;
            assert(numel(unique(y(tr)))==2&&numel(unique(y(te)))==2,'Each stratified inner fold requires both classes.');
            try
                f=cv25_fit(F(tr,:),y(tr),lambdas(a),s);
                [P(te,ix),E(te,ix)]=cv25_predict(f,F(te,:));
            catch err
                failed(ix)=failed(ix)+1;status{end+1}=err.identifier; %#ok<AGROW>
            end
        end
        if failed(ix)==0&&all(isfinite(E(:,ix)))
            e=E(:,ix);p=P(:,ix);pos=p(y==1);neg=p(y==0)';
            auc(ix)=mean((pos>neg)+.5*(pos==neg),'all');
            ll(ix)=mean(max(e,0)+log1p(exp(-abs(e)))-y.*e);
        elseif failed(ix)==0
            failed(ix)=1;status={'nonfinite_prediction'};
        end
        row(ix,:)={families{g},lambdas(a),auc(ix),ll(ix),failed(ix),strjoin(status,';')};
    end
end
assert(any(isfinite(ll)),'All inner candidates failed; no result will be silently substituted.');
ids=find(auc>=max(auc)-1e-12);lo=min(ll(ids));ids=ids(ll(ids)<=lo+1e-12);
add=ids(ids<=L);if ~isempty(add),ids=add;end
lam=cell2mat(row(ids,2));[~,j]=max(lam);ix=ids(j);
[threshold,J]=cv25_cut(y,P(:,ix));
T=struct('family',row{ix,1},'lambda',row{ix,2},'auc',auc(ix),'logloss',ll(ix), ...
    'threshold',threshold,'youden',J,'p',P(:,ix),'eta',E(:,ix),'index',ix,'failedFits',sum(failed));
T.table=cell2table(row,'VariableNames',{'Family','Lambda','InnerAUC','InnerLogLoss','FailedFits','FailureStatus'});
T.table.Selected=false(nc,1);T.table.Selected(ix)=true;
end

function rows=rf25_candidate_rows(T,name,r,k)
rows=cell(height(T.table),11);
for a=1:height(T.table)
    tt=T.table;
    rows(a,:)={name,r,k,a,tt.Family{a},tt.Lambda(a),tt.InnerAUC(a),tt.InnerLogLoss(a),tt.FailedFits(a),tt.FailureStatus{a},tt.Selected(a)};
end
end

function row=rf25_summary_row(t,name)
row={name,mean(t.AUC),std(t.AUC),min(t.AUC),max(t.AUC),mean(t.Brier),std(t.Brier),mean(t.LogLoss),mean(t.Sensitivity),mean(t.Specificity),mean(t.Accuracy),mean(t.CalibrationSlope,'omitnan'),mean(t.CalibrationIntercept,'omitnan')};
end

function rf25_export(U,A,out)
if ~isfolder(out),mkdir(out);end
writetable(U.summary,fullfile(out,'refinement_summary.csv'));
writetable(U.foldSummary,fullfile(out,'refinement_outer_tuning.csv'));
writetable(U.candidateMetrics,fullfile(out,'refinement_inner_candidate_metrics.csv'));
oo={};rep=table();coeff={};
for m=1:numel(U.models)
    z=U.models(m);t=z.repeatMetrics;t.Model=repmat({z.name},height(t),1);rep=[rep;t]; %#ok<AGROW>
    for r=1:size(z.oof,2)
        for i=1:size(A,1)
            oo(end+1,:)={z.name,i,r,U.folds(i,r),A(i,6),z.oof(i,r),z.oofEta(i,r),z.thresholdByRecord(i,r),double(z.predictedLabel(i,r))}; %#ok<AGROW>
        end
    end
    f=z.final;terms=[{'Intercept'},f.featureNames];
    for j=1:numel(terms)
        if j==1,mu=0;sd=1;else,mu=f.mu(j-1);sd=f.sd(j-1);end
        coeff(end+1,:)={z.name,f.family,terms{j},f.beta(j),mu,sd,f.betaRaw(j),f.betaInputUnits(j),f.lambda,f.threshold}; %#ok<AGROW>
    end
end
writetable(rep,fullfile(out,'refinement_repeat_metrics.csv'));
writetable(cell2table(oo,'VariableNames',{'Model','RecordID','Repeat','OuterFold','Y','OOFProbability','OOFLogit','TrainingSelectedThreshold','PredictedLabel'}),fullfile(out,'refinement_oof_predictions.csv'));
writetable(cell2table(coeff,'VariableNames',{'Model','Family','Term','BetaStandardized','TrainingMean','TrainingSD','BetaRaw_mm','BetaInput_m','Lambda','ProbabilityThreshold'}),fullfile(out,'refinement_final_coefficients.csv'));
fid=fopen(fullfile(out,'refinement_scope.txt'),'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n%s\n%s\n',U.scope,U.selectionRule,U.features);
fprintf(fid,'Two fixed exploratory strategies; all fits preserve the original five predictors. No selected outer repeat or model is promoted on the basis of its outer AUC.\n');
fprintf(fid,'Full_AUC optimizes only the penalty for inner pooled OOF AUC. Coupling_AUC selects additive vs one sigma_c*V interaction and the penalty in the same inner loop.\n');
fprintf(fid,'AUC ties use inner log loss, then the additive family, then stronger penalty. Each numerical tie tolerance is 1e-12.\n');
fprintf(fid,'The final model uses the original fixed full-data inner folds. All threshold and final fitted metrics are distinct from outer performance.\n');
fprintf(fid,'Failed inner fits: %d. Candidates with a failed fold receive AUC=-Inf and log loss=Inf; failures are exported. Outer/final fit failures abort.\n',U.totalFailedInnerFits);
fprintf(fid,'Coefficient units: BetaRaw_mm uses W/H in mm, BetaInput_m uses W/H in m. The intercept has already been corrected for centering.\n');
end

function H=cav25_auc_overview(A,R,V,U,out,opts)
% Explicit full-model performance; exact mean ROC over every fixed repeat.
y=A(:,6);m=pp25_model(V,'Full_raw');a=R.models(strcmp({R.models.name},'Full_raw'));
rows={};
rows(end+1,:)={'Full_MLE','Apparent',a.auc,NaN,a.brier,1};
rows(end+1,:)={'Full_ridge','Apparent',pp25_auc(y,m.final.p),NaN,mean((y-m.final.p).^2),1};
rows(end+1,:)={'Full_ridge','Nested_row_CV',mean(m.repeatMetrics.AUC),std(m.repeatMetrics.AUC),mean(m.repeatMetrics.Brier),size(m.oof,2)};
for j=1:numel(U.models)
    q=U.models(j);t=q.repeatMetrics;
    rows(end+1,:)={q.name,'Nested_row_CV',mean(t.AUC),std(t.AUC),mean(t.Brier),size(q.oof,2)}; %#ok<AGROW>
end
H.summary=cell2table(rows,'VariableNames',{'Model','Evaluation','AUC','SplitSD','Brier','Repeats'});
H.scope='SplitSD measures partition sensitivity, not confidence limits. Apparent MLE and CV ridge are distinct fitted estimators. No target AUC is imposed.';
writetable(H.summary,fullfile(out,'FULL_MODEL_PERFORMANCE.csv'));
H.fullMeanROC=au25_meanroc(y,m.oof);
writematrix(H.fullMeanROC,fullfile(out,'full_model_mean_ROC_FPR_TPR.csv'));
assert(abs(trapz(H.fullMeanROC(:,1),H.fullMeanROC(:,2))-mean(m.repeatMetrics.AUC))<1e-12,'Mean ROC area does not match mean repeat AUC.');
q=U.models(strcmp({U.models.name},'Full_AUC'));
H.aucTunedMeanROC=au25_meanroc(y,q.oof);
writematrix(H.aucTunedMeanROC,fullfile(out,'auc_tuned_mean_ROC_FPR_TPR.csv'));
assert(abs(trapz(H.aucTunedMeanROC(:,1),H.aucTunedMeanROC(:,2))-mean(q.repeatMetrics.AUC))<1e-12);
C=[0 .447 .698;.835 .369 0;0 .620 .451;.494 .184 .556];
if opts.MakeFigures
    au25_fullfigure(A,a,m,q,H.fullMeanROC,H.aucTunedMeanROC,out,opts,C);
    au25_comparators(A,V,out,opts,C);
    au25_refinementfigure(A,V,U,out,opts,C);
end
fid=fopen(fullfile(out,'AUC_INTERPRETATION.txt'),'w','n','UTF-8');assert(fid>0);cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'FULL MODEL = W + H + raw k + sigma_c + V.\n\n%s\n',evalc('disp(H.summary)'));
fprintf(fid,'The MLE full-data apparent AUC is %.8f. This already lies near 0.88 and is not a cross-validated AUC.\n',a.auc);
fprintf(fid,'The original ridge model has first-repeat AUC %.8f and mean repeated nested-row-CV AUC %.8f (split SD %.8f).\n',m.repeatMetrics.AUC(1),mean(m.repeatMetrics.AUC),std(m.repeatMetrics.AUC));
fprintf(fid,'The revised ROC curve averages every repeated ROC as a function of FPR; its area equals the arithmetic mean of repeat AUCs, including score ties.\n');
fprintf(fid,'Fig03 presents complete-model performance. Fig09 retains the GK/full/log-sigma comparisons as a supplementary ablation and functional-form analysis. Fig10 reports the limited AUC-directed refinement experiment.\n');
fprintf(fid,'Model/penalty selection and Youden cutoffs use inner predictions only; changing a probability cutoff cannot increase an ROC AUC.\n');
fprintf(fid,'Refinement rules are fixed before these experiments. No best repeat or seed is selected; no observations are removed. Candidate results are exploratory because this dataset has already informed model development.\n');
fprintf(fid,'AUC-directed tuning may improve ranking while worsening Brier/log loss. It is not automatically substituted for the primary probability model or its bootstrapped thresholds.\n');
fprintf(fid,'No grouping identifiers are available. Row-wise CV cannot establish new-surface or temporal generalization or remove dependence among nearby monitoring points.\n');
fprintf('FULL MODEL PERFORMANCE (apparent and internal validation are separate):\n');disp(H.summary);
end

function M=au25_meanroc(y,P)
% Preserve vertical jumps AND linear segments from tied scores. The area of
% this mean polyline is exactly the mean empirical trapezoidal ROC area.
nr=size(P,2);rr=cell(nr,1);xx=[];
for r=1:nr,rr{r}=pp25_roc(y,P(:,r));xx=[xx;rr{r}(:,1)];end %#ok<AGROW>
x=unique(xx);lo=zeros(numel(x),nr);hi=lo;
for r=1:nr
    R=rr{r};u=unique(R(:,1));yl=zeros(size(u));yh=yl;
    for k=1:numel(u),ii=R(:,1)==u(k);yl(k)=min(R(ii,2));yh(k)=max(R(ii,2));end
    for k=1:numel(x)
        same=find(u==x(k),1);
        if ~isempty(same),lo(k,r)=yl(same);hi(k,r)=yh(same);
        else
            left=find(u<x(k),1,'last');right=left+1;
            value=yh(left)+(yl(right)-yh(left))*(x(k)-u(left))/(u(right)-u(left));
            lo(k,r)=value;hi(k,r)=value;
        end
    end
end
M=zeros(2*numel(x),2);M(1:2:end,:)=[x,mean(lo,2)];M(2:2:end,:)=[x,mean(hi,2)];
end

function au25_unitaxes(ax)
xlim(ax,[0 1]);ylim(ax,[0 1]);set(ax,'XTick',0:.2:1,'YTick',0:.2:1);
line(ax,[0 1],[0 1],'Color',[.6 .63 .65],'LineStyle','--','HandleVisibility','off');
end

function au25_fullfigure(A,a,m,q,meanroc,aucroc,out,opts,C)
f=pp25_figure('Full-model performance',6.8,opts);
positions=[.09 .60 .36 .31;.60 .60 .36 .31;.09 .12 .36 .31;.60 .12 .36 .31];
y=A(:,6);ax=pp25_axes(f,positions(1,:));au25_unitaxes(ax);
roc=pp25_roc(y,a.p);h=plot(ax,roc(:,1),roc(:,2),'Color',C(1,:),'LineWidth',2.3);
xlabel(ax,'False-positive fraction');ylabel(ax,'True-positive fraction');pp25_title(ax,'(a) Full model: apparent ROC');
legend(ax,h,sprintf('MLE   AUC %.3f',a.auc),'Location','southeast','Box','off','FontSize',opts.FontSize-1);
ax=pp25_axes(f,positions(2,:));au25_unitaxes(ax);
for r=1:size(m.oof,2)
    roc=pp25_roc(y,m.oof(:,r));plot(ax,roc(:,1),roc(:,2),'Color',.76+.24*C(2,:),'LineWidth',.6,'HandleVisibility','off');
end
h=plot(ax,meanroc(:,1),meanroc(:,2),'Color',C(2,:),'LineWidth',2.4);
for r=1:size(q.oof,2)
    roc=pp25_roc(y,q.oof(:,r));plot(ax,roc(:,1),roc(:,2),'Color',.85+.15*C(3,:),'LineWidth',.5,'HandleVisibility','off');
end
h(2)=plot(ax,aucroc(:,1),aucroc(:,2),'Color',C(3,:),'LineWidth',2.4);
xlabel(ax,'False-positive fraction');ylabel(ax,'True-positive fraction');pp25_title(ax,'(b) Full model: validation ROC');
legend(ax,h,{sprintf('Log-loss tuned: %.3f',mean(m.repeatMetrics.AUC)),sprintf('AUC-tuned: %.3f',mean(q.repeatMetrics.AUC))},'Location','southeast','Box','off','FontSize',opts.FontSize-2);
ax=pp25_axes(f,positions(3,:));au25_unitaxes(ax);
cal=pp25_calibration(y,a.p,5);plot(ax,cal(:,1),cal(:,2),'o-','Color',C(1,:),'MarkerFaceColor',C(1,:),'MarkerEdgeColor','w','MarkerSize',9,'LineWidth',1.8);
xlabel(ax,'Predicted probability');ylabel(ax,'Observed cavitation fraction');pp25_title(ax,'(c) Apparent calibration');
text(ax,.06,.91,sprintf('Brier = %.3f',a.brier),'Units','normalized','FontSize',opts.FontSize);
ax=pp25_axes(f,positions(4,:));au25_unitaxes(ax);
for r=1:size(m.oof,2)
    cc=pp25_calibration(y,m.oof(:,r),5);
    plot(ax,cc(:,1),cc(:,2),'-','Color',.78+.22*C(2,:),'LineWidth',.7,'HandleVisibility','off');
end
% Pooled display only: every record has equal weight across the fixed repeats.
% It does not create n*R independent records and no binomial CI is attached.
pooledY=repmat(y,size(m.oof,2),1);cal=pp25_calibration(pooledY,m.oof(:),5);
h=plot(ax,cal(:,1),cal(:,2),'s-','Color',C(2,:),'MarkerFaceColor',C(2,:),'MarkerEdgeColor','w','MarkerSize',9,'LineWidth',1.8);
for r=1:size(q.oof,2)
    cc=pp25_calibration(y,q.oof(:,r),5);
    plot(ax,cc(:,1),cc(:,2),'-','Color',.86+.14*C(3,:),'LineWidth',.5,'HandleVisibility','off');
end
cal=pp25_calibration(repmat(y,size(q.oof,2),1),q.oof(:),5);
h(2)=plot(ax,cal(:,1),cal(:,2),'o-','Color',C(3,:),'MarkerFaceColor',C(3,:),'MarkerEdgeColor','w','MarkerSize',9,'LineWidth',1.8);
xlabel(ax,'Predicted probability');ylabel(ax,'Observed cavitation fraction');pp25_title(ax,'(d) Validation calibration');
legend(ax,h,{sprintf('Log loss: Brier %.3f',mean(m.repeatMetrics.Brier)),sprintf('AUC: Brier %.3f',mean(q.repeatMetrics.Brier))},'Location','northwest','Box','off','FontSize',opts.FontSize-2);
pp25_save(f,out,'Fig03_full_model_performance',opts);
end

function au25_comparators(A,V,out,opts,C)
f=pp25_figure('Model comparison: supplementary analysis',4.3,opts);
ax1=pp25_axes(f,[.09 .19 .37 .66]);ax2=pp25_axes(f,[.60 .19 .36 .66]);au25_unitaxes(ax1);
want={'GK','Full_raw','Full_log'};want=want(ismember(want,{V.models.name}));h=gobjects(numel(want),1);labels=cell(size(want));
for j=1:numel(want)
    m=pp25_model(V,want{j});curve=au25_meanroc(A(:,6),m.oof);
    h(j)=plot(ax1,curve(:,1),curve(:,2),'Color',C(j,:),'LineWidth',2.1);
    labels{j}=sprintf('%s: %.3f',au25_label(want{j}),mean(m.repeatMetrics.AUC));
    t=m.repeatMetrics;dx=.10*sin((1:height(t))'*2.39996323);
    plot(ax2,j+dx,t.AUC,'o','Color',C(j,:),'MarkerFaceColor',.6+.4*C(j,:),'MarkerSize',5,'HandleVisibility','off');
    line(ax2,j+[-.23 .23],[mean(t.AUC),mean(t.AUC)],'LineWidth',2.3,'Color',C(j,:));
end
xlabel(ax1,'False-positive fraction');ylabel(ax1,'True-positive fraction');pp25_title(ax1,'(a) Mean validation ROC');
legend(ax1,h,labels,'Location','southeast','Box','off','FontSize',opts.FontSize-2);
set(ax2,'XTick',1:numel(want),'XTickLabel',cellfun(@au25_label,want,'UniformOutput',false),'XTickLabelRotation',20);
xlim(ax2,[.5 numel(want)+.5]);ylim(ax2,[0 1]);ylabel(ax2,'AUC');pp25_title(ax2,'(b) Repeated validation');
pp25_save(f,out,'Fig09_model_comparison',opts);
end

function s=au25_label(name)
switch name
    case 'GK',s='W, H, k';
    case 'Full_raw',s='Full model';
    case 'Full_log',s='Full (ln \sigma_c)';
    case 'Full_AUC',s='AUC-tuned';
    case 'Coupling_AUC',s='Coupling candidate';
    otherwise,s=strrep(name,'_',' ');
end
end

function au25_refinementfigure(A,V,U,out,opts,C)
f=pp25_figure('AUC-directed model refinement',4.4,opts);
ax1=pp25_axes(f,[.09 .19 .37 .65]);ax2=pp25_axes(f,[.60 .19 .36 .65]);au25_unitaxes(ax1);
models=[{pp25_model(V,'Full_raw')},arrayfun(@(q){q},U.models)];h=gobjects(numel(models),1);labels=cell(size(models));short=labels;
for j=1:numel(models)
    m=models{j};curve=au25_meanroc(A(:,6),m.oof);
    h(j)=plot(ax1,curve(:,1),curve(:,2),'Color',C(j,:),'LineWidth',2.1);
    short{j}=au25_label(m.name);if j==1,short{j}='Log-loss tuned';end
    labels{j}=sprintf('%s: %.3f',short{j},mean(m.repeatMetrics.AUC));
    b=m.repeatMetrics.Brier;dx=.1*sin((1:numel(b))'*2.39996323);
    plot(ax2,j+dx,b,'o','Color',C(j,:),'MarkerFaceColor',.6+.4*C(j,:),'MarkerSize',5);
    line(ax2,j+[-.23 .23],[mean(b),mean(b)],'Color',C(j,:),'LineWidth',2.3);
end
base=mean(V.baseline.repeatMetrics.Brier);line(ax2,[.5 numel(models)+.5],[base base],'LineStyle','--','Color',[.45 .45 .45]);
xlabel(ax1,'False-positive fraction');ylabel(ax1,'True-positive fraction');pp25_title(ax1,'(a) Mean validation ROC');
legend(ax1,h,labels,'Location','southeast','Box','off','FontSize',opts.FontSize-2);
set(ax2,'XTick',1:numel(models),'XTickLabel',short,'XTickLabelRotation',20);
xlim(ax2,[.5 numel(models)+.5]);ylabel(ax2,'Brier score');pp25_title(ax2,'(b) Probability accuracy');
pp25_save(f,out,'Fig10_auc_refinement',opts);
end

function [threshold,achieved]=rv33_recall_cut(y,p,target)
% Fixed rule: largest cutoff achieving the target in TRAINING inner OOF only.
assert(numel(y)==numel(p)&&all(isfinite(p))&&any(y==1));
positive=sort(p(y==1),'descend');needed=ceil(target*numel(positive)-1e-12);
threshold=positive(needed);achieved=mean(p(y==1)>=threshold);
assert(achieved+1e-12>=target,'Training recall target was not attained.');
end

function D=rv33_analysis(A,V,T,U,opts,out)
% No added columns or physical cluster assignments are inferred from six columns.
if ~isfolder(out),mkdir(out);end
n=size(A,1);y=A(:,6);m=pp25_model(V,'Full_raw');f=m.final;R=size(m.oof,2);
D=struct('primaryModel','Full_raw_logloss_ridge','N',n,'Positive',sum(y),'Negative',sum(1-y));
D.scope='Contemporaneous state association; nested ROW-wise internal assessment. No pre-inception, spatial, grouped, external or safety validation.';
D.thresholdRule='Largest inner-OOF probability cutoff with training sensitivity at least TargetRecall; outer labels are used only for evaluation.';
D.referenceVelocity=opts.ReferenceVelocity;D.referenceLength=opts.ReferenceLength;D.targetRecall=opts.TargetRecall;
D.finalRecallThreshold=f.recallThreshold;D.finalInnerRecall=f.innerRecallAtTargetCut;

% Positive constant rescaling commutes with training-only standardization.
Z=[A(:,1)/opts.ReferenceLength,A(:,2)/opts.ReferenceLength,A(:,3)/opts.ReferenceVelocity^2,A(:,2)./A(:,1),A(:,4:6)];
D.dimensionlessData=array2table(Z,'VariableNames',{'W_over_L','H_over_L','k_over_U2','H_over_W','sigma_c','V','Y'});
writetable(D.dimensionlessData,fullfile(out,'dimensionless_records.csv'));
F=[1000*A(:,1:2),A(:,3:5)];scale=[1000*opts.ReferenceLength,1000*opts.ReferenceLength,opts.ReferenceVelocity^2,1,1];
Fnd=F./scale;mu=f.mu./scale;sd=f.sd./scale;
pnd=cv25_sigmoid([ones(n,1),(Fnd-mu)./sd]*f.beta);
D.dimensionlessMaxProbabilityDifference=max(abs(pnd-f.p));
assert(D.dimensionlessMaxProbabilityDifference<1e-10,'Dimensionless reparameterization changed probabilities.');
D.ratioNote='H/W is a local aspect ratio, not the effective slope of the full surface. It is exported descriptively, not added alongside H and W to the primary model.';

% Threshold and coefficient bootstrap use the SAME retuned ridge estimator.
D.ridgeCoefficients=T.ridgeCoefficientRanges;
D.ridgeCoefficients.OR=exp(D.ridgeCoefficients.Beta);
D.ridgeCoefficients.ORPercentile2_5=exp(D.ridgeCoefficients.BootstrapPercentile2_5);
D.ridgeCoefficients.ORPercentile97_5=exp(D.ridgeCoefficients.BootstrapPercentile97_5);
% exp(intercept) is baseline odds, not an odds ratio for a predictor.
D.ridgeCoefficients.OR(1)=NaN;D.ridgeCoefficients.ORPercentile2_5(1)=NaN;D.ridgeCoefficients.ORPercentile97_5(1)=NaN;
writetable(D.ridgeCoefficients,fullfile(out,'PRIMARY_ridge_coefficients_OR.csv'));
B=size(T.bootstrap.coefficients,1);bounds=NaN(B,10);point=NaN(1,10);names={'W_over_L','H_over_L','k_over_U2','sigma_c','H_over_W_fixed_W'};
physical= false(B,10);support=false(B,10);pointPhys=false(1,10);pointSupport=false(1,10);
ref=T.referenceFeatures;refHOverW=ref(2)/ref(1); %#ok<NASGU>
for vv=0:1
    ref(5)=vv;
    for q=1:5
        j=vv*5+q;qq=min(q,4);if q==5,qq=2;end
        [raw,~]=th25_crossing(f,ref,qq,f.threshold);
        [point(j),pointPhys(j),pointSupport(j)]=rv33_transform_crossing(raw,q,vv,A,opts,ref);
        for b=1:B
            if ~T.bootstrap.success(b),continue;end
            fb=f;fb.beta=T.bootstrap.coefficients(b,:)';cut=T.bootstrap.values(b,2);
            raw=th25_crossing(fb,ref,qq,cut);
            [bounds(b,j),physical(b,j),support(b,j)]=rv33_transform_crossing(raw,q,vv,A,opts,ref);
        end
    end
end
rows={};
for vv=0:1
    for q=1:5
        j=vv*5+q;x=bounds(:,j);finite=isfinite(x);ci=q7(x(finite),[.025 .5 .975]);
        rows(end+1,:)={vv,names{q},point(j),pointPhys(j),pointSupport(j),sum(finite),B-sum(finite),sum(physical(:,j)),sum(support(:,j)),ci(1),ci(2),ci(3)}; %#ok<AGROW>
    end
end
D.dimensionlessBoundaries=cell2table(rows,'VariableNames',{'V','Quantity','Crossing','PhysicalDomainValid','WithinVGroupMarginalRange','FiniteBootstrap','MissingOrFailed','PhysicalBootstrap','WithinMarginalBootstrap','Percentile2_5','Median','Percentile97_5'});
D.dimensionlessBoundaryNote='Reference W,H,k are the full-sample means, fixed across row bootstraps. Other predictors are held at reference values. H/W is obtained from a conditional H crossing at fixed W. Negative/out-of-domain crossings are reported, not clipped. Percentile ranges condition on finite crossings and do not establish physical inception or joint support.';
writetable(D.dimensionlessBoundaries,fullfile(out,'dimensionless_conditional_boundaries.csv'));
writematrix(bounds,fullfile(out,'dimensionless_boundary_bootstrap_values.csv'));

% All operating-point rules use identical outer predicted probabilities.
rules={'Fixed_0p5','Youden_inner','Recall_inner'};rows={};costrows={};oo={};foldrows={};
for r=1:R
    pp=m.oof(:,r);preds={pp>=.5,m.predictedLabel(:,r),m.recallPredicted(:,r)};
    for a=1:3
        c=cv25_confusion(y,preds{a});
        rows(end+1,:)=[{rules{a},r},num2cell(c)]; %#ok<AGROW>
        for cr=opts.CostRatios(:)'
            costrows(end+1,:)={rules{a},r,cr,(cr*c(2)+c(4))/n}; %#ok<AGROW>
        end
    end
    for i=1:n
        oo(end+1,:)={i,r,V.folds(i,r),y(i),pp(i),m.recallThresholdByRecord(i,r),m.recallInnerSensitivityByRecord(i,r),double(m.recallPredicted(i,r))}; %#ok<AGROW>
    end
    for k=1:opts.OuterFolds
        idx=V.folds(:,r)==k;th=unique(m.recallThresholdByRecord(idx,r));ach=unique(m.recallInnerSensitivityByRecord(idx,r));assert(isscalar(th)&&isscalar(ach));
        foldrows(end+1,:)={r,k,sum(idx),th,ach}; %#ok<AGROW>
    end
end
D.operatingPoints=cell2table(rows,'VariableNames',{'Rule','Repeat','TP','FN','TN','FP','Sensitivity','Specificity','Accuracy','PPV','NPV','F1','BalancedAccuracy'});
D.costs=cell2table(costrows,'VariableNames',{'Rule','Repeat','CostFNOverFP','MeanCostPerRecord'});
% Deterministic constant-label comparators on the observed case mix.
costRatios=opts.CostRatios(:);
D.costConstantBaselines=table(costRatios,repmat(sum(y==0)/n,numel(costRatios),1),costRatios*sum(y==1)/n, ...
    'VariableNames',{'CostFNOverFP','AlwaysPositiveMeanCost','AlwaysNegativeMeanCost'});
D.costBaselineNote='Constant-label comparators use the observed class proportions. They are diagnostic cost baselines, not fitted models or established engineering decision rules.';
D.recallOOF=cell2table(oo,'VariableNames',{'Record','Repeat','OuterFold','Y','Probability','TrainingSelectedCut','InnerOOFSensitivity','PredictedLabel'});
D.recallFolds=cell2table(foldrows,'VariableNames',{'Repeat','OuterFold','TestN','TrainingSelectedCut','InnerOOFSensitivity'});
writetable(D.operatingPoints,fullfile(out,'operating_points_by_repeat.csv'));
writetable(D.costs,fullfile(out,'cost_sensitivity_by_repeat.csv'));
writetable(D.costConstantBaselines,fullfile(out,'cost_constant_baselines.csv'));
writetable(D.recallOOF,fullfile(out,'recall_outer_predictions.csv'));
writetable(D.recallFolds,fullfile(out,'recall_inner_cutoffs_by_fold.csv'));
summary={};metrics={'Sensitivity','Specificity','PPV','NPV','BalancedAccuracy','TP','FN','TN','FP'};
for a=1:3
    t=D.operatingPoints(strcmp(D.operatingPoints.Rule,rules{a}),:);
    for j=1:numel(metrics)
        x=t.(metrics{j});x=x(isfinite(x));
        if isempty(x),v=[NaN NaN NaN NaN];else,v=[mean(x),std(x),min(x),max(x)];end
        summary(end+1,:)=[{rules{a},metrics{j},numel(x)},num2cell(v)]; %#ok<AGROW>
    end
end
D.operatingSummary=cell2table(summary,'VariableNames',{'Rule','Metric','FiniteRepeats','Mean','SplitSD','Minimum','Maximum'});
writetable(D.operatingSummary,fullfile(out,'operating_point_summary.csv'));
first=D.operatingPoints(strcmp(D.operatingPoints.Rule,'Recall_inner')&D.operatingPoints.Repeat==1,:);
D.firstRepeatRecallWilson=wilson(first.TP,first.TP+first.FN);
D.bestPossibleRecallWilson=wilson(sum(y),sum(y));

% Paired differences describe sensitivity to the SAME splits, not independent tests.
diffrows={};
for j=1:numel(U.models)
    t=U.models(j).repeatMetrics;primary=m.repeatMetrics;
    for r=1:R,diffrows(end+1,:)={U.models(j).name,r,t.AUC(r)-primary.AUC(r),t.Brier(r)-primary.Brier(r),t.LogLoss(r)-primary.LogLoss(r)};end %#ok<AGROW>
end
D.pairedDifferences=cell2table(diffrows,'VariableNames',{'Alternative','Repeat','DeltaAUC','DeltaBrier','DeltaLogLoss'});
writetable(D.pairedDifferences,fullfile(out,'paired_model_differences.csv'));

% Effective sample size is unidentified; these are explicit assumptions only.
rows={};
for G=opts.ScenarioGroups(:)'
    mbar=n/G;
    for rho=opts.ScenarioICC(:)'
        de=1+(mbar-1)*rho;ne=n/de;
        rows(end+1,:)={G,mbar,rho,de,ne,ne*mean(y),ne*(1-mean(y)),sqrt(de)}; %#ok<AGROW>
    end
end
D.designEffect=cell2table(rows,'VariableNames',{'AssumedEqualGroups','AssumedMeanClusterSize','AssumedICC','DesignEffect','ScenarioEffectiveN','ProportionallyScaledPositive','ProportionallyScaledNegative','HeuristicSEMultiplier'});
D.designEffectNote='Equal-sized exchangeable-cluster scenarios for a generic mean; neither ICC nor group membership is estimated. This scalar design effect is not an exact covariance correction for each logistic coefficient. Effective event counts are proportional illustrations, not observed independent events. No post-hoc observed power is claimed.';
writetable(D.designEffect,fullfile(out,'design_effect_scenarios.csv'));
D.positivePerSlope=sum(y)/5;D.negativePerSlope=sum(1-y)/5;D.nominalSlopes=5;
fid=fopen(fullfile(out,'REVIEWER33_REPORT.txt'),'w','n','UTF-8');assert(fid>0);cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Primary: five-slope Full_raw ridge, lambda tuned by inner log loss. RAW k is retained.\n%s\n\n',D.scope);
fprintf(fid,'N=%d, positives=%d, negatives=%d; nominal counts per five slopes: %.4g / %.4g (not a power calculation).\n',n,sum(y),sum(1-y),D.positivePerSlope,D.negativePerSlope);
fprintf(fid,'Reference U=%.8g m/s, L=%.8g m. Constant rescaling prediction difference=%.3g.\n%s\n\n',opts.ReferenceVelocity,opts.ReferenceLength,D.dimensionlessMaxProbabilityDifference,D.ratioNote);
fprintf(fid,'PRIMARY ridge coefficients and odds ratios (row-bootstrap retuning; exploratory ranges):\n%s\n',evalc('disp(D.ridgeCoefficients)'));
fprintf(fid,'%s\n%s\n',D.dimensionlessBoundaryNote,evalc('disp(D.dimensionlessBoundaries)'));
fprintf(fid,'Recall target=%.3f; selected final cutoff=%.8g; full-data inner recall=%.8g (NOT outer test recall).\n',opts.TargetRecall,D.finalRecallThreshold,D.finalInnerRecall);
fprintf(fid,'%s\n%s\n',D.thresholdRule,evalc('disp(D.operatingSummary)'));
fprintf(fid,'First-repeat recall Wilson interval under working row independence: [%.6g,%.6g].\n',D.firstRepeatRecallWilson);
fprintf(fid,'Even all %d positives correct would give a two-sided Wilson lower bound %.6g (working independence).\n',sum(y),D.bestPossibleRecallWilson(1));
fprintf(fid,'Repeat SD measures partition sensitivity; do not treat R*N predictions as independent observations. Cost ratios are sensitivity assumptions, not established engineering losses.\n');
fprintf(fid,'%s\n%s\n',D.costBaselineNote,evalc('disp(D.costConstantBaselines)'));
fprintf(fid,'No DeLong test is applied to correlated repeated-CV predictions or selected/reused data. Paired deltas are descriptive only.\n\n%s\n%s\n',D.designEffectNote,evalc('disp(D.designEffect)'));
if opts.MakeFigures,rv33_figures(D,opts,out);end
fprintf('Reviewer 3.3 output: %s\n',out);
end

function [value,valid,within]=rv33_transform_crossing(raw,q,vv,A,opts,ref)
value=NaN;valid=false;within=false;if ~isfinite(raw),return;end
switch q
    case 1,value=raw/(1000*opts.ReferenceLength);observed=A(A(:,5)==vv,1)/opts.ReferenceLength;valid=value>0;
    case 2,value=raw/(1000*opts.ReferenceLength);observed=A(A(:,5)==vv,2)/opts.ReferenceLength;valid=value>=0;
    case 3,value=raw/opts.ReferenceVelocity^2;observed=A(A(:,5)==vv,3)/opts.ReferenceVelocity^2;valid=value>=0;
    case 4,value=raw;observed=A(A(:,5)==vv,4);valid=true;
    case 5,value=raw/ref(1);observed=A(A(:,5)==vv,2)./A(A(:,5)==vv,1);valid=value>=0;
end
within=valid&&value>=min(observed)&&value<=max(observed);
end

function rv33_figures(D,opts,out)
C=[0 .447 .698;.835 .369 0;0 .620 .451;.494 .184 .556];rules={'Fixed_0p5','Youden_inner','Recall_inner'};labels={'0.5','Youden','Recall target'};
f=pp25_figure('Operating points',4.9,opts);ax=pp25_axes(f,[.10 .19 .38 .67]);
for j=1:3
    t=D.operatingPoints(strcmp(D.operatingPoints.Rule,rules{j}),:);
    x=j+[-.14 .14];vals=[mean(t.Sensitivity),mean(t.Specificity)];sd=[std(t.Sensitivity),std(t.Specificity)];
    for a=1:2,errorbar(ax,x(a),vals(a),sd(a),'o','Color',C(a,:),'MarkerFaceColor',C(a,:),'LineWidth',1.5,'MarkerSize',7,'CapSize',7);end
end
plot(ax,NaN,NaN,'o','Color',C(1,:),'MarkerFaceColor',C(1,:),'DisplayName','Sensitivity');plot(ax,NaN,NaN,'o','Color',C(2,:),'MarkerFaceColor',C(2,:),'DisplayName','Specificity');
allh=get(ax,'Children');legend(ax,flipud(allh(1:2)),{'Sensitivity','Specificity'},'Location','southwest','FontSize',opts.FontSize-1,'Box','off');
set(ax,'XTick',1:3,'XTickLabel',labels,'YLim',[0 1.03],'XLim',[.5 3.5]);ylabel(ax,'Classification rate');pp25_title(ax,'(a) Operating points');
ax2=pp25_axes(f,[.59 .19 .36 .67]);
for j=1:3
    means=zeros(size(opts.CostRatios));sds=means;
    for k=1:numel(opts.CostRatios),t=D.costs(strcmp(D.costs.Rule,rules{j})&D.costs.CostFNOverFP==opts.CostRatios(k),:);means(k)=mean(t.MeanCostPerRecord);sds(k)=std(t.MeanCostPerRecord);end
    errorbar(ax2,opts.CostRatios,means,sds,'-o','Color',C(j,:),'MarkerFaceColor',C(j,:),'LineWidth',1.6,'MarkerSize',6,'DisplayName',labels{j});
end
xlabel(ax2,'Cost ratio C_{FN}/C_{FP}');ylabel(ax2,'Cost per record');pp25_title(ax2,'(b) Cost sensitivity');legend(ax2,'Location','northwest','Box','off','FontSize',opts.FontSize-1);
pp25_save(f,out,'Fig11_operating_points',opts);
f=pp25_figure('Design effect scenarios',4.7,opts);ax=pp25_axes(f,[.14 .18 .78 .68]);rho=linspace(0,max(opts.ScenarioICC),101);
for j=1:numel(opts.ScenarioGroups)
    G=opts.ScenarioGroups(j);m=D.N/G;ne=D.N./(1+(m-1)*rho);
    plot(ax,rho,ne,'LineWidth',2,'Color',C(mod(j-1,size(C,1))+1,:),'DisplayName',sprintf('G = %g, m = %.2f',G,m));
end
xlabel(ax,'Assumed intracluster correlation');ylabel(ax,'Scenario effective sample size');pp25_title(ax,'Dependence sensitivity');
ylim(ax,[0 D.N*1.05]);legend(ax,'Location','northeast','Box','off','FontSize',opts.FontSize-1);
pp25_save(f,out,'Fig12_design_effect_scenarios',opts);
fid=fopen(fullfile(out,'FIGURE_CAPTIONS.txt'),'w','n','UTF-8');assert(fid>0);cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Fig11. Operating-point performance of the same primary model probabilities. Error bars show split SD across dependent repeated nested row-CV partitions, not confidence intervals. Recall cutoffs are selected using inner OOF data only. Costs use assumed relative FN/FP losses.\n');
fprintf(fid,'Fig12. Assumed equal-group exchangeable-correlation scenarios: N_eff=N/[1+(N/G-1)rho]. These curves are not estimates of actual grouping or effective sample size.\n');
end
