function [cf_,gof,out] = spgfit(daten, para0, para_point, dat_point)
%
% [cf_,gof,out] = spgfit(data, para, para_point, dat_point)
%
% cf_: parameter output
% gof: Goodness of Fit, rmse
% out: no. of iterations, etc
%
% data: data (in cell form)
% Input: [time(s) force HA emg1 emg2 etc] 
%
% para0: initial parameters
% para_point: no. of parameters
% dat_point: no of channels

global para_global para_point_global dat_point_global;

%% 
if nargin < 2
    % Fitting initial values
else
    st_ = para0;
end;
%% 
if nargin < 3 || isempty(para_point)
    para_point = 1:length(st_);
end;

dat_split = size(daten{1},2);
if nargin < 4
    dat_point = 1:dat_split-3;
end;

% st_ = para0;
para_point_global = para_point;
para_global = st_;
st_ = st_(para_point_global);

dat_point_global = dat_point;

%% 
if ~iscell(daten)
    daten = {daten};
end;

wh = length(daten);
chan = length(dat_point);
n_dt = 0; 
for i=1:wh
    n_dt=n_dt+(length(daten{i})*chan);
end;

dt = daten{1}(5,1) - daten{1}(4,1); 
t_end = (n_dt-1) * dt;
fittime = (0:dt:t_end)';

expdat=[];
for j=1:wh
    dat_j = daten{j}(:,4:dat_split);
    dat_j = dat_j(:,dat_point_global);
    expdat(end+1:end+(length(daten{j})*length(dat_point_global))) = dat_j(:);
end;
emgs=expdat';

%% Fitting

lb = [      -Inf -Inf ... % W
            -Inf -Inf ... % Wv
            -Inf -Inf ... % Wa
            -Inf -Inf ... % Wva
                 -Inf ... % r1
                 -Inf ... % r2
                 -Inf ... % A11
                 -Inf ... % A12
           0 0 0 0]; ...  % Ta, b
lb = lb(para_point);

ub = [      +Inf +Inf ... % W
            +Inf +Inf ... % Wv
            +Inf +Inf ... % Wa
            +Inf +Inf ... % Wva
                 +Inf ... % r1
                 +Inf ... % r2
                 +Inf ... % A11
                 +Inf ... % A12
      +Inf +Inf +Inf +Inf];  % Tr, Ta, b, th
  


ub = ub(para_point);

fo_ = fitoptions('method','NonlinearLeastSquares', ...
                 'MaxFunEvals', 1000, 'Display', 'iter', ...
                 'Lower',lb,'Upper',ub,'TolFun',0.001,...
                 'TolX',0.001,'MaxFunEvals',2000,'MaxIter',50);
set(fo_,'Startpoint',st_);

ok_ = isfinite(fittime) & isfinite(emgs); 
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end;

global expdata;
expdata = daten;

coeff_cell = {'W1','W2','Wv1','Wv2','Wa3','Wa4','Wva3','Wva4','r1','r2','a11','a12', ...
             'Ta1','Ta2','b1','b2'};

coeff_fit  = coeff_cell(para_point);

fct_str = 'spgmodel(fittime,[';          
for i=1:length(st_)
    fct_str = [fct_str coeff_fit{i} ',']; 
end;
fct_str = [fct_str(1:end-1) '])'];       

% Fittype
ft_ = fittype(fct_str, ...
     'dependent',{'simemg'},'independent',{'fittime'},...
     'coefficients',coeff_fit);

%% 

[cf_,gof,out] = fit(fittime(ok_),emgs(ok_),ft_,fo_); 
