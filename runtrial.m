% For loading all processed data, running SPG model, and getting predicted results
% from the model

clear all
close all

trial=[16,17,19];   % trial no.
BW=92*9.81;         %body weight

for c=1 %:10; %;  % condition of subject for which to run analysis [1:5]
	t = 1;			% counter  in loop
		for tr = 1:3    %trial no.
            eval(['load pp',num2str(c),'/tr',num2str(trial(tr))]);	% load data 
			disp(['pp',num2str(c),'/tr',num2str(trial(tr))]);	% display data 
            eval(['F_R' num2str(trial(tr)) '= Force_R(:,1);']);
            eval(['F_L' num2str(trial(tr)) '= Force_L(:,1);']);  
            eval(['HA_R' num2str(trial(tr)) '= HA_R;']);
            eval(['HA_L' num2str(trial(tr)) '= HA_L;']);  

            eval(['filter_y_RTA' num2str(trial(tr)) '= filter_y_RTA;']);
            eval(['filter_y_RSO' num2str(trial(tr)) '= filter_y_RSO;']);
            eval(['filter_y_LTA' num2str(trial(tr)) '= filter_y_LTA;']);
            eval(['filter_y_LSO' num2str(trial(tr)) '= filter_y_LSO;']);
            t = t + 1;
        end;
end;  

a = input('which side? R - 1, L - 2 ');

        if a == 1;
            eval(['Force=[F_R' num2str(trial(tr-2)) '; F_R' num2str(trial(tr-1)) '; F_R' num2str(trial(tr)) ']/BW;';]); % Force. normalised.
            eval(['HA=([HA_R' num2str(trial(tr-2)) '; HA_R' num2str(trial(tr-1)) '; HA_R' num2str(trial(tr)) ']*pi)/180;';]); % hip angle. convert to radians.
            N=length(Force);
            x1 = linspace(0,4.9999e-004*N,N)';
            data2=[x1 Force HA];
            eval([ 'data2(:,4)=[filter_y_RSO' num2str(trial(tr-2)) '; filter_y_RSO' num2str(trial(tr-1)) '; filter_y_RSO' num2str(trial(tr)) '];' ]);
            eval([ 'data2(:,5)=[filter_y_RTA' num2str(trial(tr-2)) '; filter_y_RTA' num2str(trial(tr-1)) '; filter_y_RTA' num2str(trial(tr)) '];' ]);
        else
            eval(['Force=[F_L' num2str(trial(tr-2)) '; F_L' num2str(trial(tr-1)) '; F_L' num2str(trial(tr)) ']/BW;';]); % Force. normalised.
            eval(['HA=([HA_L' num2str(trial(tr-2)) '; HA_L' num2str(trial(tr-1)) '; HA_L' num2str(trial(tr)) ']*pi)/180;';]); % hip angle. convert to radians.          
            N=length(Force);
            x1 = linspace(0,4.9999e-004*N,N)';
            data2=[x1 Force HA];
            eval([ 'data2(:,4)=[filter_y_LSO' num2str(trial(tr-2)) '; filter_y_LSO' num2str(trial(tr-1)) '; filter_y_LSO' num2str(trial(tr)) '];' ]);
            eval([ 'data2(:,5)=[filter_y_LTA' num2str(trial(tr-2)) '; filter_y_LTA' num2str(trial(tr-1)) '; filter_y_LTA' num2str(trial(tr)) '];' ]);
        end;
        
disp(sum(isnan(data2)))

%% initial parameters

w1=0.16;
w2=0.27;
wv1=0.12; 
wv2=0.37; 
wa3=-1.64;
wa4=-4.49; 
wva3=98.76; 
wva4=284.60;
r1=2.2952;
r2=0.0054;
a11=12.64;
a12=-10.35; 
Ta1=0.090;
Ta2=0.061; 
b1=14.88; 
b2=46.03;
% 
para0(1)=w1;
para0(2)=w2;
para0(3)=wv1;
para0(4)=wv2;
para0(5)=wa3;
para0(6)=wa4;
para0(7)=wva3;
para0(8)=wva4;
para0(9)=r1;
para0(10)=r2;
para0(11)=a11; 
para0(12)=a12;
para0(13)=Ta1;
para0(14)=Ta2;
para0(15)=b1;
para0(16)=b2;

%% spg
daten={data2};
para_point=1:16;

[cf_,gof,out] = spgfit(daten, para0, para_point);

parameter=coeffvalues(cf_)
cf=[gof.rmse gof.sse]
simplot(daten, parameter);
