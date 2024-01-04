clear all
close all
clc

% Import Data
fm1= load("Data\Saturation Test\Saturation.txt");
fm2= load("Data\Saturation Test\NoSaturation.txt");
traj1 = readtable('Data/Saturation Test/Saturation.csv','HeaderLines',4);
traj2= readtable('Data/Saturation Test/NoSaturation.csv','HeaderLines',4);
%%
labels= {'frame','subframe','lasi_x','lasi_y','lasi_z','rasi_x','rasi_y','rasi_z','lpsi_x','lpsi_y','lpsi_z','rpsi_x','rpsi_y','rpsi_z','lthia_x','lthia_y','lthia_z','lthi_x','lthi_y','lthi_z','lkne_x','lkne_y','lkne_z','ltiba_x','ltiba_y','ltiba_z','ltib_x','ltib_y','ltib_z','lank_x','lank_y','lank_z','lhee_x','lhee_y','lhee_z','ltoe_x','ltoe_y','ltoe_z','rthia_x','rthia_y','rthia_z','rthi_x','rthi_y','rthi_z','rkne_x','rkne_y','rkne_z','rtiba_x','rtiba_y','rtiba_z','rtib_x','rtib_y','rtib_z','rank_x','rank_y','rank_z','rhee_x','rhee_y','rhee_z','rtoe_x','rtoe_y','rtoe_z'};
traj1.Properties.VariableNames= labels;
traj2.Properties.VariableNames= labels;

% Process Forcemats
[data{1}, fm_grid{1}, lhs{1}, rhs{1}, lto{1}, rto{1}]= process_force_mats(traj1,fm1);
[data{2}, fm_grid{2}, lhs{2}, rhs{2}, lto{2}, rto{2}]= process_force_mats(traj2,fm2);

% Assign Gait Cycles
for j= 1:2
    gc= ones(height(data{j}),1);
    for i= 2:height(data{j})
        if (data{j}.lcontact(i)-data{j}.lcontact(i-1))==1
            gc(i:end)= gc(end)+1;
        end
    end
    data{j}= addvars(data{j},gc,'After','rcontact','NewVariableNames','gc');
end
clear gc

% Section into gait cycles
clear data_gc
for j= 1:2
    counter = 1;
    for i= 1:data{j}.('gc')(end)
        temp= find(data{j}.('gc')==i);
        data_gc{j}{counter,1} = data{j}(temp,:);
        counter = counter + 1;
    end
end
for j= 1:2
    data_gc{j}(1:2,:)= [];
    data_gc{j}(end-1:end)= [];
end

% Normalize wrt Time
for j= 1:2
    labels= data{j}.Properties.VariableNames;
    for i= 1:length(data_gc{j})
        temp= interp1(linspace(0,1,height(data_gc{j}{i})),data_gc{j}{i}{:,:},linspace(0,1,100));
        temp= array2table(temp);
        temp.Properties.VariableNames= labels;
        data_gc_normalized{j}{i,1}= temp;
    end 
end 

% Average GRF Data
for j= 1:2
    temp= [];
    for i= 1:30
        temp= [temp data_gc_normalized{j}{end-i+1}.force_left];
    end
    left_grf(:,j)= mean(temp,2);
end

% Plot GRF
c= 8;
figure; set(gcf,'color','w'); hold on;
plot(left_grf(:,1)/c,'LineWidth',2)
plot(left_grf(:,2),'LineWidth',2)
legend('Scaled Saturated','Unsaturated')
ylabel('Force Units')
xlabel('Percent Gait Cycle')

% Average CoP Data
for j= 1:2
    temp_x= [];
    temp_y= [];
    for i= 1:30
        temp_x= [temp_x data_gc_normalized{j}{end-i+1}.cop_x];
        temp= mean([data_gc_normalized{j}{end-i+i}.lasi_y data_gc_normalized{j}{end-i+i}.rasi_y data_gc_normalized{j}{end-i+i}.lpsi_y data_gc_normalized{j}{end-i+i}.rpsi_y ],2);
        temp_y= [temp_y data_gc_normalized{j}{end-i+1}.cop_y-temp]; 
    end
    cop_x(:,j)= mean(temp_x,2);
    cop_y(:,j)= mean(temp_y,2);
end

% Plot CoP
figure; set(gcf,'color','w'); hold on;
plot(cop_x(:,1),cop_y(:,1),'LineWidth',2)
plot(cop_x(:,2),cop_y(:,2),'LineWidth',2)
legend('Saturated','Unsaturated')
ylabel('Ant/Post')
xlabel('Med/Lat')

% Plot all CoP
figure; set(gcf,'color','w'); hold on;
for i= 1:30
    plot(data_gc_normalized{1}{end-i+1}.cop_x,data_gc_normalized{1}{end-i+1}.cop_y,'b')
end
for i= 1:30
    plot(data_gc_normalized{2}{end-i+1}.cop_x,data_gc_normalized{2}{end-i+1}.cop_y,'r')
end

% Average Ankle Position
for j= 1:2
    temp_lx= [];
    temp_ly= [];
    temp_rx= [];
    temp_ry= [];
    for i= 1:30
        temp_lx= [temp_lx data_gc_normalized{j}{end-i+1}.lank_x];
        temp_ly= [temp_ly data_gc_normalized{j}{end-i+1}.lank_y];
        temp_rx= [temp_rx data_gc_normalized{j}{end-i+1}.rank_x];
        temp_ry= [temp_ry data_gc_normalized{j}{end-i+1}.rank_y];
    end
    lank_x(:,j)= mean(temp_lx,2);
    lank_y(:,j)= mean(temp_ly,2);
    rank_x(:,j)= mean(temp_rx,2);
    rank_y(:,j)= mean(temp_ry,2);
end

% Plot Average Ankle Pos
figure; set(gcf,'color','w'); hold on;
plot(lank_x(:,1),lank_y(:,1),'b','LineWidth',2)
plot(rank_x(:,1),rank_y(:,1),'b','LineWidth',2)
plot(lank_x(:,2),lank_y(:,2),'r','LineWidth',2)
plot(rank_x(:,2),rank_y(:,2),'r','LineWidth',2)
legend('Saturated','Saturated','Unsaturated','Unsaturated')
ylabel('Ant/Post')
xlabel('Med/Lat')

% Plot ankle positions
figure; set(gcf,'color','w'); hold on;
for i= 1:30
    plot(data_gc_normalized{1}{end-i+1}.lank_x,data_gc_normalized{1}{end-i+1}.lank_y,'b')
    plot(data_gc_normalized{1}{end-i+1}.rank_x,data_gc_normalized{1}{end-i+1}.rank_y,'b')
end
for i= 1:30
    plot(data_gc_normalized{2}{end-i+1}.lank_x,data_gc_normalized{2}{end-i+1}.lank_y,'r')
    plot(data_gc_normalized{2}{end-i+1}.rank_x,data_gc_normalized{2}{end-i+1}.rank_y,'r')
end

%% Fucntions
function [data_new, fm_grid, lhs, rhs, lto, rto]= process_force_mats(traj_data,fm_data)

% Find constants and set t(1) = 0
fm_data(:,1)= fm_data(:,1)-fm_data(1,1); % force mat system time starting at 0 in milliseconds
length_fm= fm_data(end,1);
length_traj= height(traj_data)*10; % length in milliseconds of experiment

% Truncate Data
if length_traj > length_fm
    traj_data= traj_data(1:round(length_fm/10),:);
    length_traj= height(traj_data)*10;
elseif length_traj < length_fm
    ind= find(min(abs(fm_data(:,1) - length_traj)) == abs(fm_data(:,1) - length_traj)); % Find closest value to end
    fm_data= fm_data(1:ind,:);
end

frame_total= height(traj_data);

% Upsample Forcemats
fm_data= interp1(fm_data(:,1),fm_data(:,5:end),linspace(0,fm_data(end,1),length_traj/10));

for i= 1:height(traj_data)
    fm_grid{i,1}= reshape(fm_data(i,:),96,88);
    fm_grid_left{i,1}= fm_grid{i}(:,1:44);
    fm_grid_right{i,1}= fm_grid{i}(:,45:88);
end

for i= 1:height(traj_data)
    force_left(i,1)= sum(fm_grid_left{i},'all');
    force_right(i,1)= sum(fm_grid_right{i},'all');
end

signal= smooth(force_left,10);
max_value= max(signal);
[~, ind(:,1)]= findpeaks(-signal,'MinPeakHeight',-350,'MinPeakWidth',0.65,'MinPeakDistance',100);
threshold= max_value*0.03;
for i= 1:length(ind)
    temp= ind(i);
    while signal(temp)<threshold && temp < height(traj_data)
        temp= temp+1;
    end
    lhs(i,1)= temp;
end

signal= smooth(force_right,10);
max_value= max(signal);
[~, ind(:,1)]= findpeaks(-signal,'MinPeakHeight',-350,'MinPeakWidth',0.65,'MinPeakDistance',100);
threshold= max_value*0.03;
for i= 1:length(ind)
    temp= ind(i);
    while signal(temp)<threshold && temp < height(traj_data)
        temp= temp+1;
    end
    rhs(i,1)= temp;
end

[~,lto(:,1)]= findpeaks(-traj_data.ltoe_y,'MinPeakWidth',20,'MinPeakDistance',30);
[~,rto(:,1)]= findpeaks(-traj_data.rtoe_y,'MinPeakWidth',20,'MinPeakDistance',30);

if lhs(1)>lto(1)
    hs_counter = 1;
else
    hs_counter = 0;
end
to_counter = 1;
state= 1;
lcontact= zeros(height(traj_data),1);
for i= 2:height(traj_data)
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
lcontact(1)= lcontact(2);

if rhs(1)>rto(1)
    hs_counter = 1;
else
    hs_counter = 0;
end
to_counter = 1;
state= 0;
rcontact= zeros(height(traj_data),1);
for i= 2:height(traj_data)
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
rcontact(1)= rcontact(2);

traj_data= addvars(traj_data,lcontact,rcontact,'NewVariableNames',{'lcontact','rcontact'});

c= 10.21; % distance between cells in mm
x_1_offset= -7.4; % x offset of left force mat origin in Vicon frame
y_1_offset= -54.43; % y offset of left force mat origin in Vicon frame
x_2_offset= 467.5; % x offset of right force mat origin in Vicon frame
y_2_offset= -57.55; % y offset of right force mat origin in Vicon frame

% filter left forcemat data
for i= 1:frame_total
    if traj_data.lcontact(i)==0
        fm_grid_left{i} = zeros(96,44);
    else
        x_limits= [min([traj_data.lank_x(i) traj_data.ltoe_x(i)])-50 max([traj_data.lhee_x(i) traj_data.ltoe_x(i)])+75];
        y_limits= [traj_data.lhee_y(i)-10 traj_data.ltoe_y(i)+100];
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
    if traj_data.rcontact(i)==0
        fm_grid_right{i} = zeros(96,44);
    else
        x_limits= [min([traj_data.rhee_x(i) traj_data.rtoe_x(i)])-75 max([traj_data.rank_x(i) traj_data.rtoe_x(i)])+50];
        y_limits= [traj_data.rhee_y(i)-10 traj_data.rtoe_y(i)+75];
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

% calculate new cop
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
fm_add= [COP_x COP_y fm_force COP_x_left COP_y_left fm_force_left COP_x_right COP_y_right fm_force_right];
data_new= addvars(traj_data,fm_add(:,1),fm_add(:,2),fm_add(:,3),fm_add(:,4),fm_add(:,5),fm_add(:,6),fm_add(:,7),fm_add(:,8),fm_add(:,9),'NewVariableNames',{'cop_x','cop_y','force','cop_x_left','cop_y_left','force_left','cop_x_right','cop_y_right','force_right'});

end