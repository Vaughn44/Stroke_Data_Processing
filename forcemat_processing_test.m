clear all
close all
clc
%%
subject_number= 2;   
subject= ['Subject ' num2str(subject_number)];
folder= ['Data/' subject];
dur= [10 12 12 10 12 12]*60; %s

% Labels
labels= {'frame','subframe','lasi_x','lasi_y','lasi_z','rasi_x','rasi_y','rasi_z','lpsi_x','lpsi_y','lpsi_z','rpsi_x','rpsi_y','rpsi_z','lthia_x','lthia_y','lthia_z','lthi_x','lthi_y','lthi_z','lkne_x','lkne_y','lkne_z','ltiba_x','ltiba_y','ltiba_z','ltib_x','ltib_y','ltib_z','lank_x','lank_y','lank_z','lhee_x','lhee_y','lhee_z','ltoe_x','ltoe_y','ltoe_z','rthia_x','rthia_y','rthia_z','rthi_x','rthi_y','rthi_z','rkne_x','rkne_y','rkne_z','rtiba_x','rtiba_y','rtiba_z','rtib_x','rtib_y','rtib_z','rank_x','rank_y','rank_z','rhee_x','rhee_y','rhee_z','rtoe_x','rtoe_y','rtoe_z','lhip_x','lhip_y','lhip_z','rhip_x','rhip_y','rhip_z'};

% Load data
traj_data = readtable([folder '/' 'traj.csv'],'HeaderLines',4);
fm_data= load([folder '/fm.txt']);
data= [traj_data];
data.Properties.VariableNames= labels;
data= data(1:dur(subject_number)*100,:);
frame_total= height(data);
%%
fm_time= fm_data(:,1)-fm_data(1,1);
t_end= dur(subject_number)*1000;
ind= find(min(abs(fm_time - t_end)) == abs(fm_time - t_end)); % Find closest value to end of VICON data
fm_raw_shifted= fm_data(:,5:end);

fm_data_truncated= fm_raw_shifted(1:ind,:);
fm_time_truncated= fm_time(1:ind);

fm_data_upsampled= interp1(fm_time_truncated,fm_data_truncated,linspace(0,fm_time_truncated(end),frame_total));

fm_data_backup= fm_data_upsampled;
data_backup= data;
fm_data= fm_data_upsampled;

lank_x= data.lank_x;
lank_y= data.lank_y;
lhee_x= data.lhee_x;
lhee_y= data.lhee_y;
ltoe_x= data.ltoe_x;
ltoe_y= data.ltoe_y;
rank_x= data.rank_x;
rank_y= data.rank_y;
rhee_x= data.rhee_x;
rhee_y= data.rhee_y;
rtoe_x= data.rtoe_x;
rtoe_y= data.rtoe_y;
%% Run function
[a b]= forcemat_processing(fm_data,data.lank_x,data.lank_y,data.lhee_x,data.lhee_y,data.ltoe_x,data.ltoe_y,data.rank_x,data.rank_y,data.rhee_x,data.rhee_y,data.rtoe_x,data.rtoe_y);

