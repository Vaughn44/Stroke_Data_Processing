clear all
close all
clc

% Load data
fm_data= load('fm.txt');
traj_data = readtable('traj.csv','HeaderLines',4);
kin_data= readtable('kin.csv','HeaderLines',4);
kin_data= kin_data(:,3:end);
dur= 12*60; % experiment duration in seconds

% Labels
labels= {'frame','subframe','lasi_x','lasi_y','lasi_z','rasi_x','rasi_y','rasi_z','lpsi_x','lpsi_y','lpsi_z','rpsi_x','rpsi_y','rpsi_z','lthia_x','lthia_y','lthia_z','lthi_x','lthi_y','lthi_z','lkne_x','lkne_y','lkne_z','ltiba_x','ltiba_y','ltiba_z','ltib_x','ltib_y','ltib_z','lank_x','lank_y','lank_z','lhee_x','lhee_y','lhee_z','ltoe_x','ltoe_y','ltoe_z','rthia_x','rthia_y','rthia_z','rthi_x','rthi_y','rthi_z','rkne_x','rkne_y','rkne_z','rtiba_x','rtiba_y','rtiba_z','rtib_x','rtib_y','rtib_z','rank_x','rank_y','rank_z','rhee_x','rhee_y','rhee_z','rtoe_x','rtoe_y','rtoe_z','lhip_x','lhip_y','lhip_z','rhip_x','rhip_y','rhip_z','lank_absrx','lank_absry','lank_absrz','lank_rx','lank_ry','lank_rz','lfoot_prox','lfoot_proy','lfoot_proz','lhip_rx','lhip_ry','lhip_rz','lkne_rx','lkne_ry','lkne_rz','lpel_rx','lpel_ry','lpel_rz','rank_absrx','rank_absry','rank_absrz','rank_rx','rank_ry','rank_rz','rfoot_prox','rfoot_proy','rfoot_proz','rhip_rx','rhip_ry','rhip_rz','rkne_rx','rkne_ry','rkne_rz','rpel_rx','rpel_ry','rpel_rz'};

% Combine and add labels
data= [traj_data kin_data];
data.Properties.VariableNames= labels;

% Truncate Vicon Data
data= data(1:dur*100,:);
frame_total= height(data);

clear kin_data traj_data
%% Heel Strike Detection
close all
clear lhs rhs lto rto lhs_vespa rhs_vespa lhs_fm rhs_fm

% Heel strike using VESPA
threshold= 5;
lhs_vespa(:,1)= fvespa(data.lhee_z,data.lhee_y,data.frame,threshold);
rhs_vespa(:,1)= fvespa(data.rhee_z,data.rhee_y,data.frame,threshold);

% Upsample Forcemat Data
fm_time= fm_data(:,1)-fm_data(1,1);
t_end= dur*1000;
ind= find(min(abs(fm_time - t_end)) == abs(fm_time - t_end)); % Find closest value to end
fm_data_truncated= fm_data(1:ind,5:end);
fm_time_truncated= fm_time(1:ind);
fm_data_upsampled= interp1(fm_time_truncated,fm_data_truncated,linspace(0,fm_time_truncated(end),frame_total));

% Reshape Forcemats into grid
for i= 1:frame_total
    fm_grid{i,1}= reshape(fm_data_upsampled(i,:),96,88);
    fm_grid_left{i,1}= fm_grid{i}(:,1:44);
    fm_grid_right{i,1}= fm_grid{i}(:,45:88);
end

% Find "force" on left & right belts
for i= 1:frame_total
    force_left(i,1)= sum(fm_grid_left{i},'all');
    force_right(i,1)= sum(fm_grid_right{i},'all');
end

% Heel Strike using Forcemats
clear ind lhs_fm
signal= smooth(force_left,10);
[valleys(:,1), ind(:,1)]= findpeaks(-signal,'MinPeakHeight',-350,'MinPeakWidth',0.65,'MinPeakDistance',100);
valleys= -valleys;
[peaks(:,1),~]= findpeaks(signal,'MinPeakHeight',3000,'MinPeakWidth',0.65,'MinPeakDistance',100);
threshold= (mean(peaks)/1.1-mean(valleys))*.05 + mean(valleys); % finding 5% threshold
for i= 1:length(ind)
    temp= ind(i);
    while signal(temp)<threshold && temp < frame_total
        temp= temp+1;
    end
    lhs_fm(i,1)= temp;
end

clear ind rhs_fm
signal= smooth(force_right,10);
[valleys(:,1), ind(:,1)]= findpeaks(-signal,'MinPeakHeight',-350,'MinPeakWidth',0.65,'MinPeakDistance',100);
valleys= -valleys;
[peaks(:,1),~]= findpeaks(signal,'MinPeakHeight',3000,'MinPeakWidth',0.65,'MinPeakDistance',100);
threshold= (mean(peaks)/1.1-mean(valleys))*.05 + mean(valleys); % finding 5% threshold
for i= 1:length(ind)
    temp= ind(i);
    while signal(temp)<threshold && temp < frame_total
        temp= temp+1;
    end
    rhs_fm(i,1)= temp;
end

%%Validaion plots & comparing vespa & force mats
% figure; set(gcf,'color','w'); hold on;
% subplot(2,2,1); hold on;
% plot(data.lhee_y);
% plot(lhs_vespa,data.lhee_y(lhs_vespa),'rx');
% plot(lhs_fm,data.lhee_y(lhs_fm),'gx');
% legend('Sagittal Postion','F-VESPA','Force Mats')
% subplot(2,2,2); hold on;
% plot(data.rhee_y);
% plot(rhs_vespa,data.rhee_y(rhs_vespa),'rx');
% plot(rhs_fm,data.rhee_y(rhs_fm),'gx');
% legend('Sagittal Postion','F-VESPA','Force Mats')
% subplot(2,2,3); hold on;
% plot(diff(lhs_fm))
% subplot(2,2,4); hold on;
% plot(diff(rhs_fm))

% Use the force mats heel strikes moving forward
lhs= lhs_fm;
rhs= rhs_fm;

%% Toe Off Detection
[~,lto(:,1)]= findpeaks(-data.ltoe_y,'MinPeakWidth',20,'MinPeakDistance',30);
[~,rto(:,1)]= findpeaks(-data.rtoe_y,'MinPeakWidth',20,'MinPeakDistance',30);
%% Determining Contact (Stance & Swing)

% Calculate left contact
if lhs(1)>lto(1) % check if the left leg starts in stance or swing
    hs_counter = 1;
else
    hs_counter = 0;
end
to_counter = 1;
state= 1; % set a starting state (will be changed)
lcontact= zeros(frame_total,1); % preallocate
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
lcontact(1)= lcontact(2); % correct first value

% Calculate right contact
if rhs(1)>rto(1) % check if the right leg starts in stance or swing
    hs_counter = 1;
else
    hs_counter = 0;
end
to_counter = 1;
state= 0; % set a starting state (will be changed)
rcontact= zeros(frame_total,1); % preallocate
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
rcontact(1)= rcontact(2); % correct first value

data= addvars(data,lcontact,rcontact,'NewVariableNames',{'lcontact','rcontact'});

%% Process Forcemats & Add to data
% Translation constants
c= 10.21; % distance between cells in mm
x_1_offset= -7.4; % x offset of left force mat origin in Vicon frame
y_1_offset= -54.43; % y offset of left force mat origin in Vicon frame
x_2_offset= 467.5; % x offset of right force mat origin in Vicon frame
y_2_offset= -57.55; % y offset of right force mat origin in Vicon frame

% filter left forcemat data
for i= 1:frame_total
    if data.lcontact(i)==0 % if in swing set all values to 0
        fm_grid_left{i} = zeros(96,44);
    else % if in stance set all values far from the foot to 0
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
    if data.rcontact(i)==0 % if in swing set all values to 0
        fm_grid_right{i} = zeros(96,44);
    else % if in stance set all values far from the foot to 0
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
%% Scale Forcemat Data

weight= 88.45*9.81; % subject weight in N

% Recalculate "force" on left and right
for i= 1:length(fm_grid)
    force_left(i,1)= sum(fm_grid_left{i},'all');
    force_right(i,1)= sum(fm_grid_right{i},'all');
end

% Recalculate peak values
[left_force_peaks,~]= findpeaks(smooth(force_left),'MinPeakHeight',3000,'MinPeakWidth',0.65,'MinPeakDistance',100);
[right_force_peaks,~]= findpeaks((force_right),'MinPeakHeight',3000,'MinPeakWidth',0.65,'MinPeakDistance',100);

mean_force_peak= mean([left_force_peaks; right_force_peaks]); % mean between left and right
C= weight/mean_force_peak; % scaling factor

% Scale forcemat data
for i= 1:frame_total
    fm_grid_left{i}= fm_grid_left{i}*C;
    fm_grid_right{i}= fm_grid_right{i}*C;
end

% Recombine left and right into single grid 
for i= 1:frame_total
    fm_grid{i}= [fm_grid_left{i} fm_grid_right{i}];
end

%% Calculate New CoP
for t= 1:frame_total
    % set all variables to 0 to start
    moment_x = 0;
    moment_y = 0;
    total_force = 0;
    moment_x_left = 0;
    moment_y_left = 0;
    moment_x_right = 0;
    moment_y_right = 0;
    total_force_left = 0;
    total_force_right = 0;
    % Total CoP
    for y= 1:96 % for each row
        for x= 1:88 % for each column
            force = fm_grid{t}(y,x); % force value of a cell
            if x<44.5 % left belt transformation & calculation
                moment_x = moment_x + (c * (x-1) + x_1_offset) * force; % cacluate moments
                moment_y = moment_y + (c * (y-1) + y_1_offset) * force; % cacluate moments
            elseif x>44.5 % right belt transformation & calculation
                moment_x = moment_x + (c * (x-45) + x_2_offset) * force; % cacluate moments
                moment_y = moment_y + (c * (y-1) + y_2_offset) * force; % cacluate moments
            end
            total_force = total_force + force; % add forces
        end
        
        % Individual CoP
        for x= 1:44 % for each column
            force_left= fm_grid_left{t}(y,x); % force value of a cell
            moment_x_left = moment_x_left + (c * (x-1) + x_1_offset) * force_left; % cacluate moments
            moment_y_left = moment_y_left + (c * (y-1) + y_1_offset) * force_left; % cacluate moments
            total_force_left = total_force_left + force_left; % add forces
            
            force_right= fm_grid_right{t}(y,x); % force value of a cell
            moment_x_right = moment_x_right + (c * (x-1) + x_2_offset) * force_right; % cacluate moments
            moment_y_right = moment_y_right + (c * (y-1) + y_2_offset) * force_right; % cacluate moments
            total_force_right = total_force_right + force_right; % add forces
        end
    end
    COP_x(t,1) = moment_x / total_force; % calculate total cop
    COP_y(t,1) = moment_y / total_force; % calculate total cop
    fm_force(t,1) = total_force; % total force
    COP_x_left(t,1) = moment_x_left / total_force_left; % calculate left side cop
    COP_y_left(t,1) = moment_y_left / total_force_left; % calculate right side cop
    fm_force_left(t,1) = total_force_left; % total force on left side
    COP_x_right(t,1) = moment_x_right / total_force_right; % calculate right side cop
    COP_y_right(t,1) = moment_y_right / total_force_right; % calculate right side cop
    fm_force_right(t,1) = total_force_right; % total force on right side
end

% Add CoP to data
fm_add= [COP_x COP_y fm_force COP_x_left COP_y_left fm_force_left COP_x_right COP_y_right fm_force_right];
data= addvars(data,fm_add(:,1),fm_add(:,2),fm_add(:,3),fm_add(:,4),fm_add(:,5),fm_add(:,6),fm_add(:,7),fm_add(:,8),fm_add(:,9),'NewVariableNames',{'cop_x','cop_y','force','cop_x_left','cop_y_left','force_left','cop_x_right','cop_y_right','force_right'});
clear cop_x cop_y total_force COP_x COP_y fm_force fm_raw fm_grid_left fm_grid_right fm_data_upsampled fm_data