% For processing of data done with (insole) pressure sensors, Qualisys 3D motion system and EMG
% sensors
% ext trigger for Qualisys and EMG, where EMG starts 20ms earlier
% int trigger for pressure sensors and EMG
% inputs:   psensor x(no)R, x(no)L
%           EMG data
%           hip angles LHA RHA

% clear all
delete(figure(1));
delete(figure(2));
clear('LForce','LForce1','LForce2','LStep');
clear('RForce','RForce1','RForce2','RStep');

%% label EMG during walking

A=data; % EMG data

% only from trigger onwards
trig=A(:,1); % ext trigger
N=length(data); 

m1=time(N);
freq_emg=N/m1;

ptrig=A(:,2); % int trigger
i=1;

while i<=N; % start of pressure
    if ptrig(i,:)>0;
        break;
    end;
    i=i+1;
end;

j=i;
while j<N; % end of pressure
    if ptrig(j,:)<0;
        break;
    end;
    j=j+1;
end;

%% for subj no 1,2,3 - Go to LINE 318
% check label for EMG sensors
% trig=[]; % first set of lab
% trig=A(i:j,1); %channel 0
% ptrig=A(i:j,2);
% RRF=A(i:j,3);
% RGM=A(i:j,4);
% RBF=A(i:j,5);
% RVL=A(i:j,6);
% RVM=A(i:j,7);
% RTA=A(i:j,8);
% RGA=A(i:j,9);
% RSO=A(i:j,10);
% LRF=A(i:j,12);
% LGM=A(i:j,13); 
% LBF=A(i:j,14);
% LVL=A(i:j,15);
% LVM=A(i:j,16);
% LTA=A(i:j,17);
% LGA=A(i:j,18);
% LSO=A(i:j,19);

%% for subject no 4-10 - Go to LINE 341
% check label for EMG sensors
trig=[]; % first set of lab
trig=A(i:j,1); %channel 0
ptrig=A(i:j,2);
RRF=A(i:j,4);
RGM=A(i:j,5);
RBF=A(i:j,6);
RVL=A(i:j,7);
RVM=A(i:j,8);
RTA=A(i:j,9);
RGA=A(i:j,10);
RSO=A(i:j,12);
LRF=A(i:j,13);
LGM=A(i:j,14); 
LBF=A(i:j,15);
LVL=A(i:j,16);
LVM=A(i:j,17);
LTA=A(i:j,18);
LGA=A(i:j,19);
LSO=A(i:j,20);

%% Convert Force (from presssure sensor) to N*steps frames
% Output: C(interpolated data) and x(time)
% input:    RForce
% output:   Force_R (used for input)

%% Right foot force
Y1 = input('psensor trial =');
Y=eval(['x',num2str(Y1),'R']); %input
N=length(Y);
f=N*10.0503; %change to 2000 Hz? (2000/199)

x1 = linspace(1,N-1,f);
x=1:N;

x2 = 1:1:f;

C = interp1(x,Y,x1);

plot(x, Y,'o',x1, C);
RForce=C;

%% Left foot force
% Y1 = input('psensor trial =');
Y=eval(['x',num2str(Y1),'L']); %input
% Y=x11L(:,2);
N=length(Y);
f=N*10.0503; %change to 2000 frames?

x1 = linspace(1,N-1,f);
x=1:N;

x2 = 1:1:f;

C = interp1(x,Y,x1);

plot(x, Y,'o',x1, C);
LForce=C;

% Force in 2000Hz, has not determined HS and TO yet. Next cell...

%% processing psensor
% input:    LHA
%           RHA
%           loading eg x_L, x_R

pfreq=199; % freq of psensor

%% Determine HS and TO from Force (and plot)

NL=length(LForce(:,2));
L=1;
Li=1; Lj=1;
LHS=[]; LTO=[]; 
while L<=(NL-1); 
    if L<(NL-51) && LForce(L,2)<=0.5 && LForce(L+1,2)>0.5 && (LForce(L+50,2)-LForce(L+0,2))>=1;
        LHS(Li,:)=L; % start of Lheel (blue) contact
        if Li>1 && (LHS(Li,:)-LHS(Li-1,:))<1200
            LHS(Li-1,:)=LHS(Li,:);
            LHS(Li,:)=[];
            Li=Li;
        else
        Li=Li+1;
        end;
    end;
    if L>31 && L<(NL-41) && LForce(L,2)>0.5 && LForce(L+1,2)<=0.5 && (LForce(L-30,2)-LForce(L+40,2))>=1;
        LTO(Lj,:)=L; % end of Ltoe (pink) contact
        if Lj>1 && (LTO(Lj,:)-LTO(Lj-1,:))<1200
            LTO(Lj,:)=[];
            Lj=Lj;
        else
        Lj=Lj+1;
        end;
    end;
    L=L+1;
end;

NR=length(RForce(:,2));
R=1;
Ri=1; Rj=1;
RHS=[]; RTO=[]; 
while R<=(NR-1); 
    if R<(NR-51) && RForce(R,2)<=4 && RForce(R+1,2)>4 && (RForce(R+50,2)-RForce(R+1,2))>=1;
        RHS(Ri,:)=R; % start of Rheel (blue) contact
        if Ri>1 && (RHS(Ri,:)-RHS(Ri-1,:))<1200
            RHS(Ri-1,:)=RHS(Ri,:);
            RHS(Ri,:)=[];
            Ri=Ri;
        else
        Ri=Ri+1;
        end;
    end;
    
    if R>31 && R<(NR-21) && RForce(R,2)>4 && RForce(R+1,2)<=4 && (RForce(R-30,2)-RForce(R+20,2))>=1;
        RTO(Rj,:)=R; % end of Rtoe (pink) contact
        if Rj>1 && (RTO(Rj,:)-RTO(Rj-1,:))<1200
            RTO(Rj,:)=[];
            Rj=Rj;
        else
        Rj=Rj+1;
        end;
    end;
    R=R+1;
end;

%plot Force
figure(1);
LForce1(1:length(LForce(:,2)),1)=NaN;
LForce1(LHS,1)=1; %LHS
LForce2(1:length(LForce(:,2)),1)=NaN;
LForce2(LTO,1)=1; %LTO
plot(1:length(LForce),LForce(:,2),'-k',1:length(LForce1),LForce1,'o',1:length(LForce2),LForce2,'mo');ylabel('LForce'); 
%plot LHS blue, plot LTO pink
 
figure(2);
RForce1(1:length(RForce(:,2)),1)=NaN;
RForce1(RHS,1)=1; %RHS
RForce2(1:length(RForce(:,2)),1)=NaN;
RForce2(RTO,1)=1; %RTO
plot(1:length(RForce),RForce(:,2),'-k',1:length(RForce1),RForce1,'o',1:length(RForce2),RForce2,'mo');ylabel('RForce'); 
%plot RHS blue, plot RTO pink

% Manual correction? 
if (length(LHS)-length(LTO))>2;
    LeftWrong=1;
    disp('Left wrong!');
else
    LeftWrong=0;
end;

if (length(RHS)-length(RTO))>2;
    RightWrong=1;
    disp('Right wrong!');
else
    RightWrong=0;
end;

if LeftWrong==1 || RightWrong==1;
    break;
end;

%% check if RHS and RTO look about correct

    a=1;
    while a<=(length(LHS)-1); % LHS(1)<LTO(1)<LHS(2)
        if LTO(a)>LHS(a) && LTO(a)>LHS(a+1);
            disp('Check LTO(pink)-Figure1');   
            break;
        end;    
        a=a+1;
    end;

    a=1;
    while a<=(length(LTO)-1);   % LTO(1)<LHS(2)<LTO(2)
        if LHS(a)>LTO(a) && LHS(a)>LTO(a+1);
            disp('Check LHS(blue)-Figure1');
            break;
        end;
        a=a+1;
    end;

    b=1;
    while b<=(length(RHS)-1);
        if RTO(b)>RHS(b) && RTO(b)>RHS(b+1);
            disp('Check RTO(pink)-Figure2');   
            break;
        end;  
        b=b+1;
    end;

    b=1;
    while b<=(length(RTO)-1); 
        if RHS(b)<RTO(b) && RHS(b)>RTO(b+1);
            disp('Check RHS(blue)-Figure2');
            break;
        end;
        b=b+1;
    end;

%% Determine steps
% input:    LHA
%           RHA
%           loading eg x6L, x6R

%% determine steps

if (length(LTO)-length(LHS)<2);
    if LHS(1)<LTO(1) && isequal(length(LHS),length(LTO));
        LStep=[LHS LTO];
    elseif LHS(1)<LTO(1) && (length(LHS)-length(LTO)==1); 
        LStep=[LHS(1:(end-1),:) LTO(1:end,:)]; 
    elseif LHS(1)>LTO(1) && (length(LHS)-length(LTO)==1); 
        LStep=[LHS(1:(end),:) LTO(2:end-1,:)];
    elseif LHS(1)>LTO(1) && (length(LTO)-length(LHS)==1); 
        LStep=[LHS(1:(end),:) LTO(2:end,:)];    
    elseif LHS(1)>LTO(1) && isequal(length(LHS),length(LTO));
        LStep=[LHS(1:(end-1),:) LTO(2:end,:)];
    else
        LStep=[LHS(1:(end-2),:) LTO(2:end,:)];
    end;
else
    disp('Determine LStep yourself');
end;

if (length(RHS)-length(RTO)<2);
    if RHS(1)<RTO(1) && isequal(length(RHS),length(RTO));
        RStep=[RHS RTO];    
    elseif RHS(1)<RTO(1) && (length(RHS)-length(RTO)==1); 
        RStep=[RHS(1:(end-1),:) RTO(1:end,:)];
    elseif RHS(1)>RTO(1) && isequal(length(RHS),length(RTO));
        RStep=[RHS(1:(end-1),:) RTO(2:end,:)];
    else
        RStep=[RHS(1:(end),:) RTO(2:end,:)];
    end;
else 
    disp('Determine RStep yourself');
end;

RNoSteps=length(RStep); disp('No of Right strides ='); disp(RNoSteps);
LNoSteps=length(LStep); disp('No of Left strides ='); disp(LNoSteps);

%% Calculate output force (for input to SPG model)

RForce1=RForce(RStep(1,1):RStep(end,1),:); %complete steps
LForce1=LForce(LStep(1,1):LStep(end,1),:);

% only need 6 steps
f_step=input('first step='); %pick 1st step recorded
l_step=input('last step=')+1; %till 6th step recorded, RHS to next RHS - 1 stride/gait cycle
total_steps=l_step-f_step

Force_R=RForce(RStep(f_step,1):RStep(l_step,1),2:end); 
Force_L=LForce(LStep(f_step,1):LStep(l_step,1),2:end);

figure(3); plot(Force_R,'k'); hold on; plot(Force_L); ylabel('Force'); hold off;

%% determine start of EMG corresponding with Force (from insole)
%% subject no 1,2,3

% % from EMG, 
% trig=[]; 
% trig =A(i+RStep(f_step,1):i+RStep(l_step,1),1); %channel 0 in emg
% ptrig=A(i+RStep(f_step,1):i+RStep(l_step,1),2);
% RRF=A(i+RStep(f_step,1):i+RStep(l_step,1),3);
% RGM=A(i+RStep(f_step,1):i+RStep(l_step,1),4);
% RBF=A(i+RStep(f_step,1):i+RStep(l_step,1),5);
% RVL=A(i+RStep(f_step,1):i+RStep(l_step,1),6);
% RVM=A(i+RStep(f_step,1):i+RStep(l_step,1),7);
% RTA=A(i+RStep(f_step,1):i+RStep(l_step,1),8);
% RGA=A(i+RStep(f_step,1):i+RStep(l_step,1),9);
% RSO=A(i+RStep(f_step,1):i+RStep(l_step,1),10);
% LRF=A(i+LStep(f_step,1):i+LStep(l_step,1),12);
% LGM=A(i+LStep(f_step,1):i+LStep(l_step,1),13); 
% LBF=A(i+LStep(f_step,1):i+LStep(l_step,1),14);
% LVL=A(i+LStep(f_step,1):i+LStep(l_step,1),15);
% LVM=A(i+LStep(f_step,1):i+LStep(l_step,1),16);
% LTA=A(i+LStep(f_step,1):i+LStep(l_step,1),17);
% LGA=A(i+LStep(f_step,1):i+LStep(l_step,1),18);
% LSO=A(i+LStep(f_step,1):i+LStep(l_step,1),19);

%% subject no 4-10
% % from EMG, 
trig=[]; 
trig =A(i+RStep(f_step,1):i+RStep(l_step,1),1); %channel 0 in emg
ptrig=A(i+RStep(f_step,1):i+RStep(l_step,1),2);
RRF=A(i+RStep(f_step,1):i+RStep(l_step,1),4);
RGM=A(i+RStep(f_step,1):i+RStep(l_step,1),5);
RBF=A(i+RStep(f_step,1):i+RStep(l_step,1),6);
RVL=A(i+RStep(f_step,1):i+RStep(l_step,1),7);
RVM=A(i+RStep(f_step,1):i+RStep(l_step,1),8);
RTA=A(i+RStep(f_step,1):i+RStep(l_step,1),9);
RGA=A(i+RStep(f_step,1):i+RStep(l_step,1),10);
RSO=A(i+RStep(f_step,1):i+RStep(l_step,1),12);
LRF=A(i+LStep(f_step,1):i+LStep(l_step,1),13);
LGM=A(i+LStep(f_step,1):i+LStep(l_step,1),14); 
LBF=A(i+LStep(f_step,1):i+LStep(l_step,1),15);
LVL=A(i+LStep(f_step,1):i+LStep(l_step,1),16);
LVM=A(i+LStep(f_step,1):i+LStep(l_step,1),17);
LTA=A(i+LStep(f_step,1):i+LStep(l_step,1),18);
LGA=A(i+LStep(f_step,1):i+LStep(l_step,1),19);
LSO=A(i+LStep(f_step,1):i+LStep(l_step,1),20);

%% analyse EMG 

% figure(1);
% subplot(9,1,1), plot(trig,'k'); ylabel('trigger'); 
% subplot(9,1,2), plot(RRF); ylabel('Right RF');
% subplot(9,1,3), plot(RGM); ylabel('Right GMax');
% subplot(9,1,4), plot(RBF); ylabel('Right BF');
% subplot(9,1,5), plot(RVL); ylabel('Right VL');
% subplot(9,1,6), plot(RVM); ylabel('Right VM');
% subplot(9,1,7), plot(RTA); ylabel('Right TA'); 
% subplot(9,1,8), plot(RGA); ylabel('Right GA'); 
% subplot(9,1,9), plot(RSO); ylabel('Right Sol'); 
% figure(2);
% subplot(9,1,1), plot(ptrig); ylabel('p trigger');
% subplot(9,1,2), plot(LRF); ylabel('Left RF');
% subplot(9,1,3), plot(LGM); ylabel('Left GMax');
% subplot(9,1,4), plot(LBF); ylabel('Left BF');
% subplot(9,1,5), plot(LVL); ylabel('Left VL');
% subplot(9,1,6), plot(LVM); ylabel('Left VM');
% subplot(9,1,7), plot(LTA); ylabel('Left TA');
% subplot(9,1,8), plot(LGA); ylabel('Left GA');
% subplot(9,1,9), plot(LSO); ylabel('Left Sol');

%%
% Remove any DC offset of the signal,  y2 is the signal without DC offset.
% centers (to zero) data

y2_RRF=detrend(RRF);
y2_RGM=detrend(RGM);
y2_RBF=detrend(RBF);
y2_RVL=detrend(RVL);
y2_RVM=detrend(RVM);
y2_RTA=detrend(RTA);
y2_RGA=detrend(RGA);
y2_RSO=detrend(RSO);
y2_LRF=detrend(LRF);
y2_LGM=detrend(LGM);
y2_LBF=detrend(LBF);
y2_LVL=detrend(LVL);
y2_LVM=detrend(LVM);
y2_LTA=detrend(LTA);
y2_LGA=detrend(LGA);
y2_LSO=detrend(LSO);

%%
% Rectification of the EMG signal, rec_y is the rectified signal. 

rec_y_RRF=abs(y2_RRF);
rec_y_RGM=abs(y2_RGM);
rec_y_RBF=abs(y2_RBF);
rec_y_RVL=abs(y2_RVL);
rec_y_RVM=abs(y2_RVM);
rec_y_RTA=abs(y2_RTA);
rec_y_RGA=abs(y2_RGA);
rec_y_RSO=abs(y2_RSO);
rec_y_LRF=abs(y2_LRF);
rec_y_LGM=abs(y2_LGM);
rec_y_LBF=abs(y2_LBF);
rec_y_LVL=abs(y2_LVL);
rec_y_LVM=abs(y2_LVM);
rec_y_LTA=abs(y2_LTA);
rec_y_LGA=abs(y2_LGA);
rec_y_LSO=abs(y2_LSO);

% figure(1);
% subplot(9,1,1), plot(trig,'k'); ylabel('trigger'); hold off;
% subplot(9,1,2), plot(rec_y_RRF); ylabel('Right RF'); hold off;
% subplot(9,1,3), plot(rec_y_RGM); ylabel('Right GMax'); 
% subplot(9,1,4), plot(rec_y_RBF); ylabel('Right BF'); 
% subplot(9,1,5), plot(rec_y_RVL); ylabel('Right VL'); 
% subplot(9,1,6), plot(rec_y_RVM); ylabel('Right VM'); 
% subplot(9,1,7), plot(rec_y_RTA); ylabel('Right TA'); 
% subplot(9,1,8), plot(rec_y_RGA); ylabel('Right GA'); 
% subplot(9,1,9), plot(rec_y_RSO); ylabel('Right Sol'); 
% figure(2);
% subplot(9,1,1), plot(ptrig); ylabel('p trigger');  
% subplot(9,1,2), plot(rec_y_LRF,'y'); ylabel('Left RF'); hold off;
% subplot(9,1,3), plot(rec_y_LGM,'y'); ylabel('Left GMax'); hold off;
% subplot(9,1,4), plot(rec_y_LBF,'y'); ylabel('Left BF'); 
% subplot(9,1,5), plot(rec_y_LVL,'y'); ylabel('Left VL'); 
% subplot(9,1,6), plot(rec_y_LVM,'y'); ylabel('Left VM'); 
% subplot(9,1,7), plot(rec_y_LTA,'y'); ylabel('Left TA'); 
% subplot(9,1,8), plot(rec_y_LGA,'y'); ylabel('Left GA'); 
% subplot(9,1,9), plot(rec_y_LSO,'y'); ylabel('Left Sol');

% Linear Envelope of the EMG signal, Construct a low pass filter 
% of a cut off frequency, 40Hz. 
% Here the sampling frequency is 2000Hz, 5th order filter.
%%
[b,a]=butter(5,30/2000,'low');

% The next step is to filter the signals to obtain the linear envelope.
% The command filtfilt performs filtering in both directions to eliminate any phase shift of the signal.
% filt_y is the filtered signal.

filter_y_RRF=filtfilt(b,a,rec_y_RRF);
filter_y_RGM=filtfilt(b,a,rec_y_RGM);
filter_y_RBF=filtfilt(b,a,rec_y_RBF);
filter_y_RVL=filtfilt(b,a,rec_y_RVL);
filter_y_RVM=filtfilt(b,a,rec_y_RVM);
filter_y_RTA=filtfilt(b,a,rec_y_RTA);
filter_y_RGA=filtfilt(b,a,rec_y_RGA);
filter_y_RSO=filtfilt(b,a,rec_y_RSO);
filter_y_LRF=filtfilt(b,a,rec_y_LRF);
filter_y_LGM=filtfilt(b,a,rec_y_LGM);
filter_y_LBF=filtfilt(b,a,rec_y_LBF);
filter_y_LVL=filtfilt(b,a,rec_y_LVL);
filter_y_LVM=filtfilt(b,a,rec_y_LVM);
filter_y_LTA=filtfilt(b,a,rec_y_LTA);
filter_y_LGA=filtfilt(b,a,rec_y_LGA);
filter_y_LSO=filtfilt(b,a,rec_y_LSO);

figure(4);
subplot(9,1,1), plot(trig,'k'); ylabel('trigger'); hold off;
subplot(9,1,2), plot(filter_y_RRF); ylabel('Right RF'); ylim([0 0.1]); hold off;
subplot(9,1,3), plot(filter_y_RGM); ylabel('Right GMax'); ylim([0 0.2]); hold off;
subplot(9,1,4), plot(filter_y_RBF); ylabel('Right BF'); ylim([0 0.2]); hold off;
subplot(9,1,5), plot(filter_y_RVL); ylabel('Right VL'); ylim([0 0.2]); hold off;
subplot(9,1,6), plot(filter_y_RVM); ylabel('Right VM'); ylim([0 0.2]); hold off;
subplot(9,1,7), plot(filter_y_RTA); ylabel('Right TA'); ylim([0 0.5]); hold off;
subplot(9,1,8), plot(filter_y_RGA); ylabel('Right GA'); ylim([0 0.5]); hold off;
subplot(9,1,9), plot(filter_y_RSO); ylabel('Right Sol'); ylim([0 0.5]); 
xlabel('Low Pass Filtered EMG samples'); hold off;

figure(5);
subplot(9,1,1), plot(ptrig); ylabel('p trigger'); hold off;
subplot(9,1,2), plot(filter_y_LRF); ylabel('Left RF'); ylim([0 0.2]); hold off;
subplot(9,1,3), plot(filter_y_LGM); ylabel('Left GMax'); ylim([0 0.2]); hold off;
subplot(9,1,4), plot(filter_y_LBF); ylabel('Left BF'); ylim([0 0.3]); hold off;
subplot(9,1,5), plot(filter_y_LVL); ylabel('Left VL'); ylim([0 0.2]); hold off;
subplot(9,1,6), plot(filter_y_LVM); ylabel('Left VM'); ylim([0 0.2]); hold off;
subplot(9,1,7), plot(filter_y_LTA); ylabel('Left TA'); ylim([0 0.5]); hold off;
subplot(9,1,8), plot(filter_y_LGA); ylabel('Left GA'); ylim([0 0.5]); hold off;
subplot(9,1,9), plot(filter_y_LSO); ylabel('Left Sol'); ylim([0 0.5]); 
xlabel('Low Pass Filtered EMG samples'); hold off;

%% processing hip angle
% input:    LHA
%           RHA
%           loading eg x6L, x6R

% from psensor
% f_step, l_step

%% Convert HA to N*steps frames
% Output: C(interpolated data) and x(time)
% input:    RHA

%% Right hip angle
Y=RHA; %input
N=length(Y);
f=2000*20; %change to 2000 frames

x1 = linspace(1,N-1,f);
x=1:N;

x2 = 1:1:f;

C = interp1(x,Y,x1);

% plot(x, Y(:,2:4),'o',x1, C(:,2:4));
RHA2=C;

%% Left hip angle
Y=LHA; %input
N=length(Y);
f=2000*20; %change to 2000 frames

x1 = linspace(1,N-1,f);
x=1:N;

x2 = 1:1:f;

C = interp1(x,Y,x1);

plot(x, Y(:,2:4),'o',x1, C(:,2:4));
LHA2=C;

%% processing HA 

Cfreq=2000; % converted freq from above
% emg-(2*m)=qua*m,  % emg starts first, by 20 ms. 
m=0.02/(1/Cfreq); % m: no of frames slower, m=40

RHA2=RHA2((m+1):end,:); 
LHA2=LHA2((m+1):end,:);

% determine HA with strides
% f_step, l_step

HA_R=RHA2(i+RStep(f_step,1):i+RStep(l_step,1),2); %(:,2)-sagittal plane flex/ext
HA_L=LHA2(i+LStep(f_step,1):i+LStep(l_step,1),2);

figure(6); plot(HA_R,'k'); hold on; plot(HA_L); ylabel('Hip Angle'); hold off;
disp(sum(isnan(HA_R))); disp('RHA - BLack');
disp(sum(isnan(HA_L))); disp('LHA - blUE');