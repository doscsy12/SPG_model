function simout = spgmodel(fittime, para0)
% 
global expdata;  
global para_global para_point_global dat_point_global;

parameter = para_global;
parameter(para_point_global) = para0;

%% Parameters definition
W   = parameter(1:2);
Wv  = parameter(3:4);
Wa = parameter(5:6);
Wva = parameter(7:8);
r1   = parameter(9);
r2   = parameter(10);
a   = [parameter(11) parameter(12)
       parameter(12) parameter(11)]; 
% Tr (1x2) so x0 = data(1,3:end) (1x2). 
Ta   = parameter(13:14);
b   = parameter(15:16);

%% Simulation 
opt = simset('solver','ode4','SrcWorkspace','Current');      %'SrcWorkspace','Current' important!!! 
simout = [];
for wh = 1:length(expdata)
    data = expdata{wh};
    simtime = [data(1,1) data(end,1)];
    s0 = data(1,2); s1 = data(1,3);  x0 = data(1,4:end); v0 = 0; % 
    sim('spg_hip.mdl',simtime,opt);
    simemg = squeeze(simemg)';
    simemg = simemg(:,dat_point_global);
    simemg = simemg(:);
    simout(end+1:end+length(simemg),1) = simemg;
end;
