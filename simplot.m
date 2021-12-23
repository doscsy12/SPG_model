function simplot(daten, parameter, show_input, datenname)

if ~iscell(daten),                   
    for wh = 1:size(daten,3)        
        datcell{wh} = daten(:,:,wh); 
    end;
    daten = datcell;                 
end;

n_set  = length(daten);
n_chan = size(daten{1},2) - 3; % ( -3 for first 3 columns of data)

if nargin < 4
    for i=1:n_chan, datenname{i} = ['Channel ' num2str(i)]; end;
    if nargin < 3, show_input = 0; end;
end;

%% Parameters
W   = parameter(1:2);
Wv  = parameter(3:4);
Wa = parameter(5:6);
Wva = parameter(7:8);
r1   = parameter(9);
r2   = parameter(10);
a   = [parameter(11) parameter(12)
       parameter(12) parameter(11)];
 % Tr(1x2) so x0 = data(1,3:end) (1x2). 
Ta   = parameter(13:14);
b   = parameter(15:16);

%%
for set=1:n_set
    data = daten{set};
    simtime = [data(1,1) data(end,1)];
    s0 = data(1,2); s1 = data(1,3);  x0 = data(1,4:end); v0 = 0; 
    opt = simset('solver','ode4','SrcWorkspace','Current');
    sim('spg_hip.mdl',simtime,opt);
    simemg = squeeze(simemg);

    for channel=4:size(data,2);
        subplot(n_chan,n_set,(channel-4)*n_set+set);
        if show_input
            plot(data(:,1), data(:,2), 'k', ...
                 data(:,1), data(:,channel), 'k', ...
                 data(:,1), simemg(channel-3,:), 'b'),
        else
            plot(data(:,1), data(:,channel), 'k', ...
                 data(:,1), simemg(channel-3,:), 'b'),
        end;
        xlim(simtime);
        if set==1, ylabel(datenname{channel-3}), end;
        if channel==n_chan+3, xlabel(['Data ' num2str(set)]), end;
    end;
end;

