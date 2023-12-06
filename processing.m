%% Processing Code for Nov 2023 Stroke Experiments
% Analysis should be subject specific
% 
% Things to look at:
% - Step Length (Frozen) + Asymmetry
% - Stride Length
% - Anterior Step Length + Asymmetry
% - CoP Butterfly Path
% - Stance & Swing Asymmetry
% - Max Knee Flexion
% - Max Dorsiflexion
% - Joint Velocity During Swing
% - Hip Hiking
% - Hip Circumduction
% - Toe Height (min during swing)
% - Push Off Force
% - Knee Angle at HS
% - GA Push Off Activation
% - TA Swing Activation

% 

%% Import Data & Setup
clear all
close all
clc

load('data.mat')
dur= [10 12 12 10 12 12]*60; % s
paretic_side= {'L' 'R' 'L' 'R' 'L' 'R'}';
muscle_number= [4 10 10 10 10 10];
c= 10.07; % constant: mm/cell
gap= 25.84; % distance between mats (mm)
T= [12140 30100; % Transitions
    14510 36000;
    14900 36260;
    12250 30010;
    14560 36130;
    14580 36150];
TS= [0.4;0.5;0.5;0.3;0.45;0.35];
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
%% Calculate AoA & H2AS
g= [1 0]; % Ground reference
for j= 1:6
    for i=1:height(data{j})
        tempL= [data{j}.('lank_y')(i)-data{j}.('lhip_y')(i) data{j}.('lank_z')(i)-data{j}.('lhip_z')(i)];
        laoa(i,1)= acosd(dot(tempL,g)/(norm(tempL)*norm(g)));
        tempR= [data{j}.('rank_y')(i)-data{j}.('rhip_y')(i) data{j}.('rank_z')(i)-data{j}.('rhip_z')(i)];
        raoa(i,1)= acosd(dot(tempR,g)/(norm(tempR)*norm(g)));
        lh2as(i,1)= norm([data{j}.('lhip_x')(i)-data{j}.('lank_x')(i) data{j}.('lhip_y')(i)-data{j}.('lank_y')(i) data{j}.('lhip_z')(i)-data{j}.('lank_z')(i)]);
        rh2as(i,1)= norm([data{j}.('rhip_x')(i)-data{j}.('rank_x')(i) data{j}.('rhip_y')(i)-data{j}.('rank_y')(i) data{j}.('rhip_z')(i)-data{j}.('rank_z')(i)]);
    end
    temp= [laoa, lh2as, raoa, rh2as];
    data{j}= addvars(data{j}, laoa, lh2as, raoa, rh2as, 'After','rhip_z','NewVariableNames',{'laoa','lh2as','raoa','rh2as'});
    clear laoa raoa lh2as rh2as
end
%% Hip Angles Redo
v= [0 1]; % Vertical reference
for j= 1:6
    for i=1:height(data{j})
        temp= [data{j}.('lhip_y')(i)-data{j}.('lkne_y')(i) data{j}.('lhip_z')(i)-data{j}.('lkne_z')(i)];
        lhip_redo(i,1)= sign(temp(1))*acosd(dot(temp,v)/(norm(temp)*norm(v)));
    end
    for i=1:height(data{j})
        temp= [data{j}.('rhip_y')(i)-data{j}.('rkne_y')(i) data{j}.('rhip_z')(i)-data{j}.('rkne_z')(i)];
        rhip_redo(i,1)= sign(temp(1))*acosd(dot(temp,v)/(norm(temp)*norm(v)));
    end
    data{j}.('lhip_rx')= lhip_redo;
    data{j}.('rhip_rx')= rhip_redo;
    clear lhip_redo rhip_redo
end
%% Assign Gait Cycles
for j= 1:6
    gc= ones(height(data{j}),1);
    if paretic_side{j}=='L'
        for i= 2:height(data{j})
            if (data{j}.rcontact(i)-data{j}.rcontact(i-1))==1
                gc(i:end)= gc(end)+1;
            end
        end
    else
        for i= 2:height(data{j})
            if (data{j}.lcontact(i)-data{j}.lcontact(i-1))==1
                gc(i:end)= gc(end)+1;
            end
        end
    end
    data{j}= addvars(data{j},gc,'After','rcontact','NewVariableNames','gc');
end
clear gc
%% Sectioning into Gait Cycles
clear data_gc
for j= 1:6
    counter1 = 1;
    counter2 = 1;
    counter3 = 1;
    for i= 1:data{j}.('gc')(end)
        temp= find(data{j}.('gc')==i);
        if temp(1) < T(j,1)
            data_gc{j,1}{counter1,1} = data{j}(temp,:);
            counter1 = counter1 + 1;
        elseif temp(1) < T(j,2)
            data_gc{j,2}{counter2,1} = data{j}(temp,:);
            counter2 = counter2 + 1;
        else
            data_gc{j,3}{counter3,1} = data{j}(temp,:);
            counter3 = counter3 + 1;
        end
        
    end
end
%% Normalize wrt Time

% Baseline
for j= 1:6
    labels= data{j}.Properties.VariableNames;
    for i= 1:length(data_gc{j,1})
        temp= interp1(linspace(0,1,height(data_gc{j,1}{i})),data_gc{j,1}{i}{:,:},linspace(0,1,100));
        temp= array2table(temp);
        temp.Properties.VariableNames= labels;
        data_gc_normalized{j,1}{i,1}= temp;
    end 
end 
% Adaptation
for j= 1:6
    labels= data{j}.Properties.VariableNames;
    for i= 1:length(data_gc{j,2})
        temp= interp1(linspace(0,1,height(data_gc{j,2}{i})),data_gc{j,2}{i}{:,:},linspace(0,1,100));
        temp= array2table(temp);
        temp.Properties.VariableNames= labels;
        data_gc_normalized{j,2}{i,1}= temp;
    end 
end 
% Observation
for j= 1:6
    labels= data{j}.Properties.VariableNames;
    for i= 1:length(data_gc{j,3})
        temp= interp1(linspace(0,1,height(data_gc{j,3}{i})),data_gc{j,3}{i}{:,:},linspace(0,1,100));
        temp= array2table(temp);
        temp.Properties.VariableNames= labels;
        data_gc_normalized{j,3}{i,1}= temp;
    end 
end 
%% Outlier Detection
for j= 1:6
    for i=1:length(data_gc_normalized{j,1})
        A{1,1}(i,:)= data_gc_normalized{j,1}{i}.('laoa')';
        A{2,1}(i,:)= data_gc_normalized{j,1}{i}.('raoa')';
        A{3,1}(i,:)= data_gc_normalized{j,1}{i}.('lh2as')';
        A{4,1}(i,:)= data_gc_normalized{j,1}{i}.('rh2as')';
%         A{5,1}(i,:)= EMGDATA2000{j,1}{i}.('LTA')';
%         A{6,1}(i,:)= EMGDATA2000{j,1}{i}.('RTA')';
%         A{7,1}(i,:)= EMGDATA2000{j,1}{i}.('LGA')';
%         A{8,1}(i,:)= EMGDATA2000{j,1}{i}.('RGA')';
    end
    [~,~,temp]= outlier_detection(A);
    outliers{j,1}= find(temp);
    clear A
    for i=1:length(data_gc_normalized{j,2})
        A{1,1}(i,:)= data_gc_normalized{j,2}{i}.('laoa')';
        A{2,1}(i,:)= data_gc_normalized{j,2}{i}.('raoa')';
        A{3,1}(i,:)= data_gc_normalized{j,2}{i}.('lh2as')';
        A{4,1}(i,:)= data_gc_normalized{j,2}{i}.('rh2as')';
%         A{5,1}(i,:)= EMGDATA2000{j,2}{i}.('LTA')';
%         A{6,1}(i,:)= EMGDATA2000{j,2}{i}.('RTA')';
%         A{7,1}(i,:)= EMGDATA2000{j,2}{i}.('LGA')';
%         A{8,1}(i,:)= EMGDATA2000{j,2}{i}.('RGA')';
    end
    [~,~,temp]= outlier_detection(A);
    outliers{j,2}= find(temp);
    clear A
    for i=1:length(data_gc_normalized{j,3})
        A{1,1}(i,:)= data_gc_normalized{j,3}{i}.('laoa')';
        A{2,1}(i,:)= data_gc_normalized{j,3}{i}.('raoa')';
        A{3,1}(i,:)= data_gc_normalized{j,3}{i}.('lh2as')';
        A{4,1}(i,:)= data_gc_normalized{j,3}{i}.('rh2as')';
%         A{5,1}(i,:)= EMGDATA2000{j,3}{i}.('LTA')';
%         A{6,1}(i,:)= EMGDATA2000{j,3}{i}.('RTA')';
%         A{7,1}(i,:)= EMGDATA2000{j,3}{i}.('LGA')';
%         A{8,1}(i,:)= EMGDATA2000{j,3}{i}.('RGA')';
    end
    [~,~,temp]= outlier_detection(A);
    outliers{j,3}= find(temp);
    clear A
end

% Remove Outliers
for j= 1:6
    for k= 1:3
        data_gc{j,k}(outliers{j,k})= [];
        data_gc_normalized{j,k}(outliers{j,k})= [];
    end
end
%% Find Section Length
for j= 1:6
    for k= 1:3
        L(j,k)= length(data_gc{j,k}); 
    end
end
trans= L(:,1);
trans(:,2)= L(:,1)+L(:,2);
Len= sum(L,2);
        
%% Calculate Parameters of Interest
clear par
% Frozen Step Length
for k= 1:3
    for j= 1:6
        if paretic_side{j} == 'L'
            for i= 1:L(j,k)-1
                temp= find(diff(data_gc{j,k}{i}.lcontact)==1); % Left HS
                temp= temp(1);
                par.right_step_length{j,k}(i,1)= (data_gc{j,k}{i}.('rank_y')(end) - data_gc{j,k}{i}.('lank_y')(temp)) + TS(j)*(height(data_gc{j,k}{i})-temp)*10;
                
                temp= find(diff(data_gc{j,k}{i+1}.lcontact)==1);
                temp= temp(1);
                par.left_step_length{j,k}(i,1)= (data_gc{j,k}{i+1}.('lank_y')(temp) - data_gc{j,k}{i}.('rank_y')(end)) + TS(j)*temp*10;
                par.step_length_diff{j,k}(i,1)= par.left_step_length{j,k}(i,1)-par.right_step_length{j,k}(i,1);
            end
        elseif paretic_side{j} == 'R'
            for i= 1:L(j,k)-1
                temp= find(diff(data_gc{j,k}{i}.rcontact)==1); % Right HS
                temp= temp(1);
                par.left_step_length{j,k}(i,1)= (data_gc{j,k}{i}.('lank_y')(end) - data_gc{j,k}{i}.('rank_y')(temp)) + TS(j)*(height(data_gc{j,k}{i})-temp)*10;
                
                temp= find(diff(data_gc{j,k}{i+1}.rcontact)==1);
                temp= temp(1);
                par.right_step_length{j,k}(i,1)= (data_gc{j,k}{i+1}.('rank_y')(temp) - data_gc{j,k}{i}.('lank_y')(end)) + TS(j)*temp*10;
                par.step_length_diff{j,k}(i,1)= par.left_step_length{j,k}(i,1)-par.right_step_length{j,k}(i,1);
            end
        end
    end
end

% Stride Length
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            if paretic_side{j} == 'R'            
                par.stride_length{j,k}(i,1) = data_gc{j,k}{i}.lank_y(end) - data_gc{j,k}{i}.lank_y(1) + TS(j)*height(data_gc{j,k}{i})*10;
            elseif paretic_side{j} == 'L'
                par.stride_length{j,k}(i,1) = data_gc{j,k}{i}.rank_y(end) - data_gc{j,k}{i}.rank_y(1) + TS(j)*height(data_gc{j,k}{i})*10;            
            end
        end     
    end
end

% Anterior Step Length
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            if paretic_side{j} == 'R'   
                CoM= [data_gc{j,k}{i}{end,3:5}; data_gc{j,k}{i}{end,6:8}; data_gc{j,k}{i}{end,9:11}; data_gc{j,k}{i}{end,12:14}];
                CoM= mean(CoM,1);
                par.left_anterior{j,k}(i,1)= data_gc{j,k}{i}.('lhee_y')(end)-CoM(2);
                temp= find(diff(data_gc{j,k}{i}.rcontact)==1);
                temp= temp(1);
                CoM= [data_gc{j,k}{i}{temp,3:5}; data_gc{j,k}{i}{temp,6:8}; data_gc{j,k}{i}{temp,9:11}; data_gc{j,k}{i}{temp,12:14}];
                CoM= mean(CoM,1);
                par.right_anterior{j,k}(i,1)= data_gc{j,k}{i}.('rhee_y')(temp)-CoM(2);
                par.anterior_diff{j,k}(i,1)= par.left_anterior{j,k}(i,1)-par.right_anterior{j,k}(i,1);          
            elseif paretic_side{j} == 'L'
                CoM= [data_gc{j,k}{i}{end,3:5}; data_gc{j,k}{i}{end,6:8}; data_gc{j,k}{i}{end,9:11}; data_gc{j,k}{i}{end,12:14}];
                CoM= mean(CoM,1);
                par.right_anterior{j,k}(i,1)= data_gc{j,k}{i}.('rhee_y')(end)-CoM(2);
                temp= find(diff(data_gc{j,k}{i}.lcontact)==1);
                temp= temp(1);
                CoM= [data_gc{j,k}{i}{temp,3:5}; data_gc{j,k}{i}{temp,6:8}; data_gc{j,k}{i}{temp,9:11}; data_gc{j,k}{i}{temp,12:14}];
                CoM= mean(CoM,1);
                par.left_anterior{j,k}(i,1)= data_gc{j,k}{i}.('lhee_y')(temp)-CoM(2);
                par.anterior_diff{j,k}(i,1)= par.left_anterior{j,k}(i,1)-par.right_anterior{j,k}(i,1);          
            end
        end
    end
end

% Stance Time
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            l_gc= height(data_gc{j,k}{i});
            par.left_stance{j,k}(i,1)= length(find(data_gc{j,k}{i}.lcontact==1))/l_gc;
            par.right_stance{j,k}(i,1)= length(find(data_gc{j,k}{i}.rcontact==1))/l_gc;
            par.stance_diff{j,k}(i,1)= par.left_stance{j,k}(i)-par.right_stance{j,k}(i);
        end
    end
end

% Push Off
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            temp= find(data_gc{j,k}{i}.rcontact == 1); % right stance
            par.right_ga{j,k}(i,1)= max(data_gc{j,k}{i}.rga(temp));
            par.right_force{j,k}(i,1)= mean(data_gc{j,k}{i}.force_right(temp));
            temp= find(data_gc{j,k}{i}.lcontact == 1); % left stance
            par.left_ga{j,k}(i,1)= max(data_gc{j,k}{i}.lga(temp));
            par.left_force{j,k}(i,1)= mean(data_gc{j,k}{i}.force_left(temp));
        end
    end
end

% Heel Strike Parameters
for k= 1:3
    for j= 1:6
        if paretic_side{j} == 'L'
            for i= 1:L(j,k)
                temp= find(diff(data_gc{j,k}{i}.lcontact)==1); % Left HS
                temp= temp(1);
                par.left_knee_hs{j,k}(i,1)= data_gc{j,k}{i}.lkne_rx(temp);
                par.right_knee_hs{j,k}(i,1)= data_gc{j,k}{i}.rkne_rx(end);
            end
        elseif paretic_side{j} == 'R'
            for i= 1:L(j,k)-1
                temp= find(diff(data_gc{j,k}{i}.rcontact)==1); % Right HS
                temp= temp(1);
                par.right_knee_hs{j,k}(i,1)= data_gc{j,k}{i}.rkne_rx(temp);
                par.left_knee_hs{j,k}(i,1)= data_gc{j,k}{i}.lkne_rx(end);
            end
        end
    end
end

% Misc Parameters
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            par.max_left_knee_flexion{j,k}(i,1)= max(data_gc{j,k}{i}.lkne_rx);
            par.max_right_knee_flexion{j,k}(i,1)= max(data_gc{j,k}{i}.rkne_rx);
            
            % During Left Swing
            temp= find(data_gc{j,k}{i}.lcontact == 0); % left swing
            par.max_left_dorsiflexion{j,k}(i,1)= max(data_gc{j,k}{i}.lank_rx(temp));
            par.max_left_hip_hike{j,k}(i,1)= max(data_gc{j,k}{i}.lhip_z(temp));
            [~,ind]= min(abs(data_gc{j,k}{i}.lhip_y(temp)-data_gc{j,k}{i}.lank_y(temp)));
            par.left_circumduction{j,k}(i,1)= data_gc{j,k}{i}.lhip_x(temp(ind))-data_gc{j,k}{i}.lank_x(temp(ind));
            par.left_ta_activation{j,k}(i,1)= mean(data_gc{j,k}{i}.lta(temp));
            par.left_hip_velo{j,k}(i,1)= max(diff(data_gc{j,k}{i}.lhip_rx(temp))*100);
            par.left_knee_velo{j,k}(i,1)= min(diff(data_gc{j,k}{i}.lkne_rx(temp))*100);
            par.left_toe_clear{j,k}(i,1)= max(data_gc{j,k}{i}.ltoe_z(temp));
            
            % During Right Swing            
            temp= find(data_gc{j,k}{i}.rcontact == 0); % right swing
            par.max_right_dorsiflexion{j,k}(i,1)= max(data_gc{j,k}{i}.rank_rx(temp));
            par.max_right_hip_hike{j,k}(i,1)= max(data_gc{j,k}{i}.rhip_z(temp));
            [~,ind]= min(abs(data_gc{j,k}{i}.rhip_y(temp)-data_gc{j,k}{i}.rank_y(temp)));
            par.right_circumduction{j,k}(i,1)= data_gc{j,k}{i}.rank_x(temp(ind))-data_gc{j,k}{i}.rhip_x(temp(ind));
            par.right_ta_activation{j,k}(i,1)= mean(data_gc{j,k}{i}.rta(temp));
            par.right_hip_velo{j,k}(i,1)= max(diff(data_gc{j,k}{i}.rhip_rx(temp))*100);
            par.right_knee_velo{j,k}(i,1)= min(diff(data_gc{j,k}{i}.rkne_rx(temp))*100);
            par.right_toe_clear{j,k}(i,1)= max(data_gc{j,k}{i}.rtoe_z(temp));
        end
    end
end

% CoP Path
gc_window= 30;
for j= 1:6
    clear cop
    for i= 1:gc_window
        % Late Baseline
        com_x= mean([data_gc_normalized{j,1}{end-gc_window+i}.lasi_x data_gc_normalized{j,1}{end-gc_window+i}.rasi_x data_gc_normalized{j,1}{end-gc_window+i}.lpsi_x data_gc_normalized{j,1}{end-gc_window+i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,1}{end-gc_window+i}.lasi_y data_gc_normalized{j,1}{end-gc_window+i}.rasi_y data_gc_normalized{j,1}{end-gc_window+i}.lpsi_y data_gc_normalized{j,1}{end-gc_window+i}.rpsi_y ],2);
        cop.x1(:,i)= (c * data_gc_normalized{j,1}{end-gc_window+i}.cop_x - 11.1) - com_x;
        cop.y1(:,i)= (c * data_gc_normalized{j,1}{end-gc_window+i}.cop_y - 53.7) - com_y;
        % Early Adaptation
        com_x= mean([data_gc_normalized{j,2}{i}.lasi_x data_gc_normalized{j,2}{i}.rasi_x data_gc_normalized{j,2}{i}.lpsi_x data_gc_normalized{j,2}{i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,2}{i}.lasi_y data_gc_normalized{j,2}{i}.rasi_y data_gc_normalized{j,2}{i}.lpsi_y data_gc_normalized{j,2}{i}.rpsi_y ],2);
        cop.x2(:,i)= (c * data_gc_normalized{j,2}{i}.cop_x - 11.1) - com_x;
        cop.y2(:,i)= (c * data_gc_normalized{j,2}{i}.cop_y - 53.7) - com_y;
        % Late Adaptation
        com_x= mean([data_gc_normalized{j,2}{end-gc_window+i}.lasi_x data_gc_normalized{j,2}{end-gc_window+i}.rasi_x data_gc_normalized{j,2}{end-gc_window+i}.lpsi_x data_gc_normalized{j,2}{end-gc_window+i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,2}{end-gc_window+i}.lasi_y data_gc_normalized{j,2}{end-gc_window+i}.rasi_y data_gc_normalized{j,2}{end-gc_window+i}.lpsi_y data_gc_normalized{j,2}{end-gc_window+i}.rpsi_y ],2);
        cop.x3(:,i)= (c * data_gc_normalized{j,2}{end-gc_window+i}.cop_x - 11.1) - com_x;
        cop.y3(:,i)= (c * data_gc_normalized{j,2}{end-gc_window+i}.cop_y - 53.7) - com_y;
        % Early Observation
        com_x= mean([data_gc_normalized{j,3}{i}.lasi_x data_gc_normalized{j,3}{i}.rasi_x data_gc_normalized{j,3}{i}.lpsi_x data_gc_normalized{j,3}{i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,3}{i}.lasi_y data_gc_normalized{j,3}{i}.rasi_y data_gc_normalized{j,3}{i}.lpsi_y data_gc_normalized{j,3}{i}.rpsi_y ],2);
        cop.x4(:,i)= (c * data_gc_normalized{j,3}{i}.cop_x - 11.1) - com_x;
        cop.y4(:,i)= (c * data_gc_normalized{j,3}{i}.cop_y - 53.7) - com_y;
        % Late Observation
        com_x= mean([data_gc_normalized{j,3}{end-gc_window+i}.lasi_x data_gc_normalized{j,3}{end-gc_window+i}.rasi_x data_gc_normalized{j,3}{end-gc_window+i}.lpsi_x data_gc_normalized{j,3}{end-gc_window+i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,3}{end-gc_window+i}.lasi_y data_gc_normalized{j,3}{end-gc_window+i}.rasi_y data_gc_normalized{j,3}{end-gc_window+i}.lpsi_y data_gc_normalized{j,3}{end-gc_window+i}.rpsi_y ],2);
        cop.x5(:,i)= (c * data_gc_normalized{j,3}{end-gc_window+i}.cop_x - 11.1) - com_x;
        cop.y5(:,i)= (c * data_gc_normalized{j,3}{end-gc_window+i}.cop_y - 53.7) - com_y;
    end
    par.cop_path{j,1}= [mean(cop.x1,2) mean(cop.y1,2)];
    par.cop_path{j,1}= [par.cop_path{j,1};par.cop_path{j,1}(1,:)];
    par.cop_path{j,2}= [mean(cop.x2,2) mean(cop.y2,2)];
    par.cop_path{j,2}= [par.cop_path{j,2};par.cop_path{j,2}(1,:)];
    par.cop_path{j,3}= [mean(cop.x3,2) mean(cop.y3,2)];
    par.cop_path{j,3}= [par.cop_path{j,3};par.cop_path{j,3}(1,:)];
    par.cop_path{j,4}= [mean(cop.x4,2) mean(cop.y4,2)];
    par.cop_path{j,4}= [par.cop_path{j,4};par.cop_path{j,4}(1,:)];
    par.cop_path{j,5}= [mean(cop.x5,2) mean(cop.y5,2)];
    par.cop_path{j,5}= [par.cop_path{j,5};par.cop_path{j,5}(1,:)];
end


%% Last Thing Before Plots
return

%% Plot
close all
clc

% Setup Plot
j= 1; % Subject Number
figure; set(gcf,'color','w','Position',[1 74 1440 723]); hold on;
subplot(4,4,1); hold on;
if paretic_side{j} == 'L'
    main_title= sgtitle(['Subject ' num2str(j)  ' | \color[rgb]{0 .4471 .7647}Left Side Paretic \color{black}| ' '\color[rgb]{.8 .2078 .1451}Right Side Perturbed \color{black}| '  num2str(dur(j)/60) ' minutes | ' num2str(TS(j)) ' m/s']);
else
    main_title= sgtitle(['Subject ' num2str(j)  ' | \color[rgb]{.8 .2078 .1451}Right Side Paretic \color{black}| ' '\color[rgb]{0 .4471 .7647}Left Side Perturbed \color{black}| '  num2str(dur(j)/60) ' minutes | ' num2str(TS(j)) ' m/s']);
end
set(main_title,'FontWeight','bold','FontSize',22);

title_input= 'Frozen Step Length';
left_data_input= par.left_step_length(j,:);
right_data_input= par.right_step_length(j,:);
diff_data_input= par.step_length_diff(j,:);
ylabel_input= 'Ankle to Ankle (mm)';
plot_number_input= 1;
double_axis_plot(left_data_input,right_data_input,diff_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Stride Length';
data_input= par.stride_length(j,:);
ylabel_input= 'Ankle to Ankle (mm)';
plot_number_input= 2;
single_plot(data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Anterior Step Length';
left_data_input= par.left_anterior(j,:);
right_data_input= par.right_anterior(j,:);
diff_data_input= par.anterior_diff(j,:);
ylabel_input= 'CoM to Ankle (mm)';
plot_number_input= 3;
double_axis_plot(left_data_input,right_data_input,diff_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Stance Time';
left_data_input= par.left_stance(j,:);
right_data_input= par.right_stance(j,:);
diff_data_input= par.stance_diff(j,:);
ylabel_input= '% Gait Cycle';
plot_number_input= 4;
double_axis_plot(left_data_input,right_data_input,diff_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Max Knee Flexion during Swing';
left_data_input= par.max_left_knee_flexion(j,:);
right_data_input= par.max_right_knee_flexion(j,:);
ylabel_input= 'Flexion Angle (degrees)';
plot_number_input= 5;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Max Dorsiflexion during Swing';
left_data_input= par.max_left_dorsiflexion(j,:);
right_data_input= par.max_right_dorsiflexion(j,:);
ylabel_input= 'Flexion Angle (degrees)';
plot_number_input= 6;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Knee Angle at Heel Strike';
left_data_input= par.left_knee_hs(j,:);
right_data_input= par.right_knee_hs(j,:);
ylabel_input= 'Flexion Angle (degrees)';
plot_number_input= 7;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Center of Pressure';
data_input= par.cop_path(j,:);
plot_number_input= 8;
cop_plot(data_input,title_input,plot_number_input)

title_input= 'Hip Hiking';
left_data_input= par.max_left_hip_hike(j,:);
right_data_input= par.max_right_hip_hike(j,:);
ylabel_input= 'Hip Height (mm)';
plot_number_input= 9;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Hip Circumduction';
left_data_input= par.left_circumduction(j,:);
right_data_input= par.right_circumduction(j,:);
ylabel_input= 'Distance from Hip (mm)';
plot_number_input= 10;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Toe Clearance';
left_data_input= par.left_toe_clear(j,:);
right_data_input= par.right_toe_clear(j,:);
ylabel_input= 'Toe Height (mm)';
plot_number_input= 11;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Stance Force';
left_data_input= par.left_force(j,:);
right_data_input= par.right_force(j,:);
ylabel_input= 'Avg Force (units)';
plot_number_input= 12;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'GA Activation at Push Off';
left_data_input= par.left_ga(j,:);
right_data_input= par.right_ga(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 13;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'TA Activation during Swing';
left_data_input= par.left_ta_activation(j,:);
right_data_input= par.right_ta_activation(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 14;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Hip Velocity during Swing';
left_data_input= par.left_hip_velo(j,:);
right_data_input= par.right_hip_velo(j,:);
ylabel_input= 'Angular Velocity (deg/s)';
plot_number_input= 15;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Knee Velocity during Swing';
left_data_input= par.left_knee_velo(j,:);
right_data_input= par.right_knee_velo(j,:);
ylabel_input= 'Angular Velocity (deg/s)';
plot_number_input= 16;
normal_plot(left_data_input,right_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

%% CoP Path
close all
span = 6000; % (frames) 60sec
figure; hold on; set(gcf,'color','w');
for j= 1:6
    subplot(2,3,j); hold on;
    com_x= mean([data{j}.lasi_x data{j}.rasi_x data{j}.lpsi_x data{j}.rpsi_x],2);
    com_y= mean([data{j}.lasi_y data{j}.rasi_y data{j}.lpsi_y data{j}.rpsi_y],2);
    cop_x= c*data{j}.cop_x-11.1;
    cop_y= c*data{j}.cop_y-53.7;
    cop_com= [cop_x-com_x cop_y-com_y];
    temp= T(j,1)-span:T(j,1);
    plot(cop_com(temp,1),cop_com(temp,2)) % baseline
%     temp= T(j,1):T(j,1)+span;
%     plot(cop_com(temp,1),cop_com(temp,2)) % early adapt
%     temp= T(j,2)-span:T(j,2);
%     plot(cop_com(temp,1),cop_com(temp,2)) % late adapt
    temp= T(j,2):T(j,2)+span;
    plot(cop_com(temp,1),cop_com(temp,2)) % early obs
    temp= length(com_x)-span:length(com_x);
    plot(cop_com(temp,1),cop_com(temp,2)) % late obs
    legend('Baseline','Early Obs','Late Obs')
    title(['Subject ' num2str(j) ' CoP Path (' paretic_side{j} ' Side Affected)'])
    ylabel('Position wrt. CoM (mm)')
    xlabel('Position wrt. CoM (mm)')
end

k= 20;
figure;
% plot(cop_com(:,1),cop_com(:,2))
plot(movmean(cop_com(:,1),20),movmean(cop_com(:,2),20))
%% Functions
function cop_plot(only_data,plot_title,subplot_number)
% normal_plot(left_data,right_data,transitions,title,ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 50;

par1x= only_data{1}(:,1);
par2x= only_data{2}(:,1);
par3x= only_data{3}(:,1);
par4x= only_data{4}(:,1);
par5x= only_data{5}(:,1);
par1y= only_data{1}(:,2);
par2y= only_data{2}(:,2);
par3y= only_data{3}(:,2);
par4y= only_data{4}(:,2);
par5y= only_data{5}(:,2);

par1xf= smooth(par1x,SPAN,'rloess');
par2xf= smooth(par2x,SPAN,'rloess');
par3xf= smooth(par3x,SPAN,'rloess');
par4xf= smooth(par4x,SPAN,'rloess');
par5xf= smooth(par5x,SPAN,'rloess');
par1yf= smooth(par1y,SPAN,'rloess');
par2yf= smooth(par2y,SPAN,'rloess');
par3yf= smooth(par3y,SPAN,'rloess');
par4yf= smooth(par4y,SPAN,'rloess');
par5yf= smooth(par5y,SPAN,'rloess');

colororder('default');

subplot(4,4,subplot_number); hold on;
% plot(par1xf,par1yf,'LineWidth',2)
% plot(par2xf,par2yf,'LineWidth',2)
% plot(par3xf,par3yf,'LineWidth',2)
% plot(par4xf,par4yf,'LineWidth',2)
% plot(par5xf,par5yf,'LineWidth',2)
plot(par1x,par1y,'LineWidth',2)
plot(par2x,par2y,'LineWidth',2)
plot(par3x,par3y,'LineWidth',2)
plot(par4x,par4y,'LineWidth',2)
plot(par5x,par5y,'LineWidth',2)
axis equal
legend('Baseline','Early Adapt','Late Adapt','Early Obs','Late Obs')
title(plot_title)
end

function single_plot(only_data,transitions,x_limit,plot_title,plot_ylabel,subplot_number)
% normal_plot(left_data,right_data,transitions,title,ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 50;

par1= only_data{1};
par2= only_data{2};
par3= only_data{3};
par1f= smooth(par1,SPAN,'rloess');
par2f= smooth(par2,SPAN,'rloess');
par3f= smooth(par3,SPAN,'rloess');
parf= [par1f; par2f; par3f];
par= [par1; par2; par3];

subplot(4,4,subplot_number); hold on;
plot(1:length(parf),parf,'LineWidth',2,'color',purple,'HandleVisibility','off')
a= axis;
axis([0 x_limit a(3) a(4)])
axis manual
plot([transitions(1) transitions(1)],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([transitions(2) transitions(2)],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot(1:length(par),par,'LineWidth',1,'color',purplet,'HandleVisibility','off')
plot(1:length(parf),parf,'LineWidth',2,'color',purple)
% legend('Left','Right')
title(plot_title)
% xlabel('Gait Cycle Number')
ylabel(plot_ylabel)
end

function double_axis_plot(left_data,right_data,diff_data,transitions,x_limit,plot_title,plot_ylabel,subplot_number)
% normal_plot(left_data,right_data,transitions,title,ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 50;

par1L= left_data{1};
par2L= left_data{2};
par3L= left_data{3};
par1Lf= smooth(par1L,SPAN,'rloess');
par2Lf= smooth(par2L,SPAN,'rloess');
par3Lf= smooth(par3L,SPAN,'rloess');
parLf= [par1Lf; par2Lf; par3Lf];
parL= [par1L; par2L; par3L];

par1R= right_data{1};
par2R= right_data{2};
par3R= right_data{3};
par1Rf= smooth(par1R,SPAN,'rloess');
par2Rf= smooth(par2R,SPAN,'rloess');
par3Rf= smooth(par3R,SPAN,'rloess');
parRf= [par1Rf; par2Rf; par3Rf];
parR= [par1R; par2R; par3R];

par1D= diff_data{1};
par2D= diff_data{2};
par3D= diff_data{3};
par1Df= smooth(par1D,SPAN,'rloess');
par2Df= smooth(par2D,SPAN,'rloess');
par3Df= smooth(par3D,SPAN,'rloess');
parDf= [par1Df; par2Df; par3Df];
parD= [par1D; par2D; par3D];

colororder([0 0 0;
            purple(1:3)]) % axis colors

subplot(4,4,subplot_number); hold on;
yyaxis left
hold on
plot(1:length(parLf),parLf,'LineWidth',2,'LineStyle','-','color',blue,'Marker','none','HandleVisibility','off')
plot(1:length(parRf),parRf,'LineWidth',2,'LineStyle','-','color',red,'Marker','none','HandleVisibility','off')
a= axis;
axis([0 x_limit a(3) a(4)])
axis manual
yyaxis right
plot(1:length(parDf),parDf,'LineWidth',2,'LineStyle','-.','color',purple,'Marker','none','HandleVisibility','off')
axis manual
yyaxis left
plot([transitions(1) transitions(1)],[-10000 10000],'LineWidth',2,'LineStyle','-','color',[.3 .3 .3],'HandleVisibility','off')
plot([transitions(2) transitions(2)],[-10000 10000],'LineWidth',2,'LineStyle','-','color',[.3 .3 .3],'HandleVisibility','off')
plot(1:length(parL),parL,'LineWidth',1,'LineStyle','-','color',bluet,'Marker','none','HandleVisibility','off')
plot(1:length(parR),parR,'LineWidth',1,'LineStyle','-','color',redt,'Marker','none','HandleVisibility','off')
yyaxis right
plot(1:length(parR),parD,'LineWidth',1,'LineStyle','-','color',purplet,'Marker','none','HandleVisibility','off')
yyaxis left
plot(1:length(parLf),parLf,'LineWidth',2,'LineStyle','-','color',blue,'Marker','none')
plot(1:length(parRf),parRf,'LineWidth',2,'LineStyle','-','color',red,'Marker','none')
yyaxis right
plot(1:length(parDf),parDf,'LineWidth',2,'LineStyle','-.','color',purple,'Marker','none')
yyaxis left
if subplot_number==4
    legend('Left','Right','L-R','location','NorthEastOutside')
end
title(plot_title)
% xlabel('Gait Cycle Number')
ylabel(plot_ylabel)
end

function normal_plot(left_data,right_data,transitions,x_limit,plot_title,plot_ylabel,subplot_number)
% normal_plot(left_data,right_data,transitions,title,ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 50;

par1L= left_data{1};
par2L= left_data{2};
par3L= left_data{3};
par1Lf= smooth(par1L,SPAN,'rloess');
par2Lf= smooth(par2L,SPAN,'rloess');
par3Lf= smooth(par3L,SPAN,'rloess');
parLf= [par1Lf; par2Lf; par3Lf];
parL= [par1L; par2L; par3L];

par1R= right_data{1};
par2R= right_data{2};
par3R= right_data{3};
par1Rf= smooth(par1R,SPAN,'rloess');
par2Rf= smooth(par2R,SPAN,'rloess');
par3Rf= smooth(par3R,SPAN,'rloess');
parRf= [par1Rf; par2Rf; par3Rf];
parR= [par1R; par2R; par3R];


subplot(4,4,subplot_number); hold on;
plot(1:length(parLf),parLf,'LineWidth',2,'color',blue,'HandleVisibility','off')
plot(1:length(parRf),parRf,'LineWidth',2,'color',red,'HandleVisibility','off')
a= axis;
axis([0 x_limit a(3) a(4)])
axis manual
plot([transitions(1) transitions(1)],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([transitions(2) transitions(2)],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot(1:length(parL),parL,'LineWidth',1,'color',bluet,'HandleVisibility','off')
plot(1:length(parR),parR,'LineWidth',1,'color',redt,'HandleVisibility','off')
plot(1:length(parLf),parLf,'LineWidth',2,'color',blue)
plot(1:length(parRf),parRf,'LineWidth',2,'color',red)
% legend('Left','Right')
title(plot_title)
if ismember(subplot_number,[13 14 15 16])
    xlabel('Gait Cycle Number')
end
ylabel(plot_ylabel)
end