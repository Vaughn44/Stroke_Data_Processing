%% VST Experiment Pre-processing
% Import data X
% Filter Emgs X
% Downsample Forcemats X
% Add EMGs to Data X
% Add gap to forcemats X
% Translate Forcemats
% Filter Forcemats
% Upsample Forcemats
% Combine data
% FVESPA
% Cut wrt. gait cycle
    
clear all
close all;
clc


dur= [10 12 12 10 12 12]*60; %s
paretic_side= {'L' 'R' 'L' 'R' 'L' 'R'};
muscle_number= [4 10 10 10 10 10];

%%%%%% Vaughn's thoughts %%%%%%%%
% Need these 4 files in a folder called "Subject_Name" (example: "Subject 5")
% vst_file= 'vst_output.txt';   (vst data)
% traj_file= 'traj.csv';        (marker data)
% kin_file= 'kin.csv';          (kinematic data from model output)
% emg_file= 'emg.csv';          (emg data)
%%%
% 2 files will be added to the folder upon running this code: the marker
% data & kinematic data in 1 table, & the emg data in another table

% Description of the VST Data Columns:
% Column 1:     Time stamp                                  (sec)
% Column 2:     Treadmill Speed PWM Value                   (0-255)
% Column 3:     Treadmill Speed                             (cm/s)
% Column 4:     (Sent) Desired Position of Linear Track     (cm)
% Column 5:     Current Position of Linear Track            (cm)
% Column 6:     Treadmill's Inclination Angle               (deg)
% Column 7:     Boolean value for experiment                (0/1)
% Column 8:     Frame Number                                (integer)
% Column 9:     X coordinate of Left Heel Marker            (mm)
% Column 10:    Y coordinate of Left Heel Marker            (mm)
% Column 11:    Z coordinate of Left Heel Marker            (mm)
% Column 12:    Gait Cycle Counter                          (integer)
% Column 13:    Frame of Last Heel-Strike Event             (integer)
% Column 14:    Desired Stiffness of Treadmill              (N/m)
% Column 15:    (Raw) Desired Position of Linear Track      (cm)

% Column 4 has the values of the desired position of the linear track that
% were actually sent to the linear track. Column 15 has the raw value of
% the desired position of the linear track as they are calculated. However,
% if the values are less than 0 or higher than 38, e.g. Inf, they won't be
% sent to the linear track. Column 15 is for debugging purposes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0
%% Load Data
% Specify the desired name of the final output file
subject= ['Subject ' num2str(subject_number)];
folder= ['Data/' subject];
output_file_name = 'data.mat';

% Labels
labels= {'frame','subframe','lasi_x','lasi_y','lasi_z','rasi_x','rasi_y','rasi_z','lpsi_x','lpsi_y','lpsi_z','rpsi_x','rpsi_y','rpsi_z','lthia_x','lthia_y','lthia_z','lthi_x','lthi_y','lthi_z','lkne_x','lkne_y','lkne_z','ltiba_x','ltiba_y','ltiba_z','ltib_x','ltib_y','ltib_z','lank_x','lank_y','lank_z','lhee_x','lhee_y','lhee_z','ltoe_x','ltoe_y','ltoe_z','rthia_x','rthia_y','rthia_z','rthi_x','rthi_y','rthi_z','rkne_x','rkne_y','rkne_z','rtiba_x','rtiba_y','rtiba_z','rtib_x','rtib_y','rtib_z','rank_x','rank_y','rank_z','rhee_x','rhee_y','rhee_z','rtoe_x','rtoe_y','rtoe_z','lhip_x','lhip_y','lhip_z','rhip_x','rhip_y','rhip_z','lank_absrx','lank_absry','lank_absrz','lank_rx','lank_ry','lank_rz','lfoot_prox','lfoot_proy','lfoot_proz','lhip_rx','lhip_ry','lhip_rz','lkne_rx','lkne_ry','lkne_rz','lpel_rx','lpel_ry','lpel_rz','rank_absrx','rank_absry','rank_absrz','rank_rx','rank_ry','rank_rz','rfoot_prox','rfoot_proy','rfoot_proz','rhip_rx','rhip_ry','rhip_rz','rkne_rx','rkne_ry','rkne_rz','rpel_rx','rpel_ry','rpel_rz'};
emg_labels= {'frame' 'subframe' 'lta' 'rta' 'lga' 'rga' 'lva' 'rva' 'lrf' 'rrf' 'lbf' 'rbf'};

% Load the file generated from the Vicon Nexus
traj_data = readtable([folder '/' 'traj.csv'],'HeaderLines',4);
kin_data= readtable([folder '/' 'kin.csv'],'HeaderLines',4);
kin_data= kin_data(:,3:end);
emg_data= readtable([folder '/' 'emg.csv'],'HeaderLines',4);
fm_data= readtable([folder '/' 'fm.txt']);
if size(fm_data,2)>8452
    fm_data(:,8453:end)= [];
end
hr_data = readtable([folder '/' 'hr.csv'],'HeaderLines',3);
hr_data = hr_data{1:dur(subject_number),3};
hr_data = repelem(hr_data,100);

% Combine and add labels
data= [traj_data kin_data];
data.Properties.VariableNames= labels;
emg_data.Properties.VariableNames= emg_labels(1:muscle_number(subject_number)+2);

clear kin_data traj_data
%% Truncate Data
data= data(1:dur(subject_number)*100,:);
data= addvars(data,hr_data,'NewVariableNames','hr');
emg_data= emg_data(1:dur(subject_number)*2000,:);
fm_time= fm_data(:,1);
t_end= fm_time{1,1}+dur(subject_number)*1000;
[ind,~]= find(fm_time{:,1}>t_end);
fm_data= fm_data(1:ind(1),:);
frame_total= height(data);

clear t_end ind fm_time hr_data
%% FVESPA & Toe Off
close all
clear lhs rhs lto rto
threshold= [0 5 6.005 0 10 0];
lhs(:,1)= fvespa(data.lhee_z,data.lhee_y,data.frame,threshold(subject_number));
rhs(:,1)= fvespa(data.rhee_z,data.rhee_y,data.frame,threshold(subject_number));
figure; subplot(1,2,1); plot(diff(lhs)); subplot(1,2,2); plot(diff(rhs));
[~,lto(:,1)]= findpeaks(-data.ltoe_y,'MinPeakWidth',0.65,'MinPeakDistance',30);
[~,rto(:,1)]= findpeaks(-data.rtoe_y,'MinPeakWidth',20,'MinPeakDistance',30);
% alert;
% return
if subject_number==1
    lhs(1)= 1;
    rhs(1)= 1;
elseif subject_number==2
    rhs(1)= 1;
elseif subject_number==3
    lhs(1)= 1;
    rhs(1)= 1;
elseif subject_number==4
    lhs(1)= 1;
    rhs(1)= 1;
elseif subject_number==5
    rhs(1)= 1;
elseif subject_number==6
    lhs(1)= 1;
    rhs(1)= 1; 
end

lcontact= 1;
rcontact= 1;
hs_counter = 1;
to_counter = 1;
state= 1;
lcontact= zeros(frame_total,1);
for i= 2:frame_total
    if any(lhs==i)
        hs_counter = hs_counter + 1;
    end
    if any(lto==i)
        to_counter = to_counter + 1;
    end
    
    if hs_counter == to_counter
        state = 1;
    elseif to_counter == hs_counter + 1
        state = 0;
    else
        'error: heel strike & toe off'
        break
    end
    lcontact(i,1) = state;
end
hs_counter = 1;
to_counter = 1;
state= 1;
rcontact= zeros(frame_total,1);
for i= 2:frame_total
    if any(rhs==i)
        hs_counter = hs_counter + 1;
    end
    if any(rto==i)
        to_counter = to_counter + 1;
    end
    
    if hs_counter == to_counter
        state = 1;
    elseif to_counter == hs_counter + 1
        state = 0;
    else
        'error: heel strike & toe off'
        break
    end
    rcontact(i,1) = state;
end
data= addvars(data,lcontact,rcontact,'NewVariableNames',{'lcontact','rcontact'});

%% Process EMGs and add to data
fs= 2000;
fnyq= fs/2;
fcuthigh= 30;
fcutlow= 300;
j= muscle_number(subject_number);
window= 200; % was 400
% Raw -> Signal -> Rect -> Env -> Filt -> mvc

% Could add preallocating for speed

for i=1:j
    emg_raw(:,i)= emg_data{:,i+2};
end
L= length(emg_raw);
f= fs*(0:(L/2))/L;
[b,a]=butter(4,[fcuthigh,fcutlow]/fnyq,'bandpass');
for i= 1:j
    emg_signal(:,i)= filtfilt(b,a,emg_raw(:,i));
end
for i= 1:j
    emg_rect(:,i)= abs(emg_signal(:,i));
end
for i= 1:j
    emg_env(:,i)= sqrt(movmean((emg_rect(:,i).^2),window));
end
for i= 1:j
    mvc(i)= max(emg_env(:,i));
end
for i= 1:j/2
    mvc(2*i-1:2*i)= max([mvc(2*i-1) mvc(2*i)]);
end
[b,a]= butter(4,5/fnyq,'low');
for i= 1:j
    emg_filt(:,i)= filtfilt(b,a,emg_env(:,i));
end
for i= 1:j
    emg_mvc(:,i)= (emg_filt(:,i)./mvc(i)).*100;
end
emg_data_filtered= emg_data;
for i= 1:j
    emg_data_filtered{:,i+2}= emg_mvc(:,i);
end
clear f fnyq fs fcuthigh fcutlow i j L mvc emg_raw emg_signal emg_rect emg_env emg_filt emg_mvc a b

% Downsample EMG & Add to DATA
temp= emg_data_filtered{:,3:end};
temp= interp1(1:length(temp),temp,linspace(1,length(temp),length(temp)/20));
if muscle_number(subject_number)==4
    data= addvars(data, temp(:,1), temp(:,2), temp(:,3), temp(:,4),'NewVariableNames',emg_labels(3:muscle_number(subject_number)+2));
elseif muscle_number(subject_number)==10
    data= addvars(data, temp(:,1), temp(:,2), temp(:,3), temp(:,4), temp(:,5), temp(:,6), temp(:,7), temp(:,8), temp(:,9), temp(:,10),'NewVariableNames',emg_labels(3:muscle_number(subject_number)+2));
else
    'error: number of EMG channels'
    return
end
%% Process Forcemats (add gap) & Add to data
% Calibrate & Upsample Forcemats
C= 0.1729; % scales forcemats to N 
temp= fm_data{:,5:end}.*C;
fm_data_upsampled= interp1(1:length(temp),temp,linspace(1,length(temp),frame_total));

for i= 1:frame_total
    fm_grid{i,1}= reshape(fm_data_upsampled(i,:),96,88);
    fm_grid_left{i,1}= fm_grid{i}(:,1:44);
    fm_grid_right{i,1}= fm_grid{i}(:,45:88);
end

c= 10.21; % distance between cells in mm
x_1_offset= -7.4; % x offset of left force mat origin in Vicon frame
y_1_offset= -54.43; % y offset of left force mat origin in Vicon frame
x_2_offset= 467.5; % x offset of right force mat origin in Vicon frame
y_2_offset= -57.55; % y offset of right force mat origin in Vicon frame

% filter left forcemat data
for i= 1:frame_total
    if data.lcontact(i)==0
        fm_grid_left{i} = zeros(96,44);
    else
        x_limits= [min([data.lank_x(i) data.ltoe_x(i)])-50 max([data.lhee_x(i) data.ltoe_x(i)])+75];
        y_limits= [data.lhee_y(i)-10 data.ltoe_y(i)+100];
        x_limits= round((x_limits - x_1_offset) / c + 1);
        y_limits= round((y_limits - y_1_offset) / c + 1);
        fm_grid_left{i}(:,1:x_limits(1)) = 0; % all cells too far left
        fm_grid_left{i}(:,x_limits(2):end) = 0; % all cells too far right
        fm_grid_left{i}(1:y_limits(1),:) = 0; % all cells too far backward
        fm_grid_left{i}(y_limits(2):end,:) = 0; % all cells too far forward
    end
end

% filter right forcemat data
for i= 1:frame_total
    if data.rcontact(i)==0
        fm_grid_right{i} = zeros(96,44);
    else
        x_limits= [min([data.rhee_x(i) data.rtoe_x(i)])-75 max([data.rank_x(i) data.rtoe_x(i)])+50];
        y_limits= [data.rhee_y(i)-10 data.rtoe_y(i)+75];
        if x_limits(1)<x_2_offset
            x_limits(1)= x_2_offset;
        end
        x_limits= round((x_limits - x_2_offset) / c + 1);
        y_limits= round((y_limits - y_2_offset) / c + 1);
        fm_grid_right{i}(:,1:x_limits(1)) = 0; % all cells too far left
        fm_grid_right{i}(:,x_limits(2):end) = 0; % all cells too far right
        fm_grid_right{i}(1:y_limits(1),:) = 0; % all cells too far backward
        fm_grid_right{i}(y_limits(2):end,:) = 0; % all cells too far forward
    end
end
    
for i= 1:frame_total
    fm_grid{i}= [fm_grid_left{i} fm_grid_right{i}];
end

%% CoP Calcs
% calculate new total cop
for t= 1:frame_total
    cop_x = 0;
    cop_y = 0;
    total_force = 0;
    cop_x_left = 0;
    cop_y_left = 0;
    cop_x_right = 0;
    cop_y_right = 0;
    total_force_left = 0;
    total_force_right = 0;
    for y= 1:96
        % total
        for x= 1:88
            force = fm_grid{t}(y,x);
            if x<44.5
                cop_x = cop_x + (c * (x-1) + x_1_offset) * force;
                cop_y = cop_y + (c * (y-1) + y_1_offset) * force;
            elseif x>44.5
                cop_x = cop_x + (c * (x-45) + x_2_offset) * force;
                cop_y = cop_y + (c * (y-1) + y_2_offset) * force;
            end
            total_force = total_force + force;
        end
        
        % left & right
        for x= 1:44
            force_left= fm_grid_left{t}(y,x);
            cop_x_left = cop_x_left + (c * (x-1) + x_1_offset) * force_left;
            cop_y_left = cop_y_left + (c * (y-1) + y_1_offset) * force_left;
            total_force_left = total_force_left + force_left;
            
            force_right= fm_grid_right{t}(y,x);
            cop_x_right = cop_x_right + (c * (x-1) + x_2_offset) * force_right;
            cop_y_right = cop_y_right + (c * (y-1) + y_2_offset) * force_right;
            total_force_right = total_force_right + force_right;
        end
    end
    COP_x(t,1) = cop_x / total_force;
    COP_y(t,1) = cop_y / total_force;
    fm_force(t,1) = total_force;
    COP_x_left(t,1) = cop_x_left / total_force_left;
    COP_y_left(t,1) = cop_y_left / total_force_left;
    fm_force_left(t,1) = total_force_left;
    COP_x_right(t,1) = cop_x_right / total_force_right;
    COP_y_right(t,1) = cop_y_right / total_force_right;
    fm_force_right(t,1) = total_force_right;
end
% fm_cop= [COP_x COP_y];
fm_add= [COP_x COP_y fm_force COP_x_left COP_y_left fm_force_left COP_x_right COP_y_right fm_force_right];
data= addvars(data,fm_add(:,1),fm_add(:,2),fm_add(:,3),fm_add(:,4),fm_add(:,5),fm_add(:,6),fm_add(:,7),fm_add(:,8),fm_add(:,9),'NewVariableNames',{'cop_x','cop_y','force','cop_x_left','cop_y_left','force_left','cop_x_right','cop_y_right','force_right'});
clear cop_x cop_y total_force COP_x COP_y fm_force fm_raw fm_grid_left fm_grid_right fm_data_upsampled fm_data
%% Export
save([folder '/' output_file_name],'data')
save([folder '/fm_' output_file_name],'fm_grid')
subject_number
end
%% Combine data into one .mat file
clear all
close all
clc

for i= 1:6
    load(['Data/Subject ' num2str(i) '/data.mat'])
    Data{i,1}= data;
    clear data
end

data= Data;
clear Data i;
save('Data/data.mat')


%% Combine fm_data into one .mat file
clear all
close all
clc

for i= 1:6
    load(['Data/Subject ' num2str(i) '/fm_data.mat'])
    temp{i,1}= fm_grid;
    clear data fm_grid
end

fm_data= temp;
clear Data i;
save('Data/fm_data.mat')

alert;