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
%% Import Data & Setup
clear all
close all
clc

load('Data/data.mat')
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
SPAN= 30;
weight= [56 88.45 81.64 72.57 117 74.84]*9.81;
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
%% Calculate Limb Orientation Angle
v= [0 -1]; % Vertical reference
for j= 1:6
    for i=1:height(data{j})
        tempL= [data{j}.('ltoe_y')(i)-data{j}.('lhip_y')(i) data{j}.('ltoe_z')(i)-data{j}.('lhip_z')(i)];
        lorientation(i,1)= sign(tempL(1))*acosd(dot(tempL,v)/(norm(tempL)*norm(v)));
        tempR= [data{j}.('rtoe_y')(i)-data{j}.('rhip_y')(i) data{j}.('rtoe_z')(i)-data{j}.('rhip_z')(i)];
        rorientation(i,1)= sign(tempR(1))*acosd(dot(tempR,v)/(norm(tempR)*norm(v)));
    end
    data{j}= addvars(data{j},lorientation,rorientation,'After','rhip_z','NewVariableNames',{'lorientation','rorientation'});
    clear lorientation rorientation
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
%% Knee Angles Redo
for j= 1:6
    for i= 1:height(data{j})
        hip= [data{j}.('lhip_y')(i) data{j}.('lhip_z')(i)];
        knee= [data{j}.('lkne_y')(i) data{j}.('lkne_z')(i)];
        ankle= [data{j}.('lank_y')(i) data{j}.('lank_z')(i)];
        h2k= hip-knee;
        k2a= knee-ankle;
        temp1= atan2d(k2a(1),k2a(2));
        temp2= atan2d(h2k(1),h2k(2));
        lknee_redo(i,1)= temp1-temp2;

        hip= [data{j}.('rhip_y')(i) data{j}.('rhip_z')(i)];
        knee= [data{j}.('rkne_y')(i) data{j}.('rkne_z')(i)];
        ankle= [data{j}.('rank_y')(i) data{j}.('rank_z')(i)];
        h2k= hip-knee;
        k2a= knee-ankle;
        temp1= atan2d(k2a(1),k2a(2));
        temp2= atan2d(h2k(1),h2k(2));
        rknee_redo(i,1)= temp1-temp2;
    end
end
%% Ankle Angles Redo
for j= 1:6
    for i= 1:height(data{j})
        knee= [data{j}.('lkne_y')(i) data{j}.('lkne_z')(i)];
        ankle= [data{j}.('lank_y')(i) data{j}.('lank_z')(i)];
        heel= [data{j}.('lhee_y')(i) data{j}.('lhee_z')(i)];
        toe= [data{j}.('ltoe_y')(i) data{j}.('ltoe_z')(i)];
        k2a= knee-ankle;
        h2t= toe-heel;
        temp1= atan2d(k2a(1),k2a(2));
        temp2= atan2d(h2t(1),h2t(2));
        lankle_redo(i,1)= 90+temp1-temp2;

        knee= [data{j}.('rkne_y')(i) data{j}.('rkne_z')(i)];
        ankle= [data{j}.('rank_y')(i) data{j}.('rank_z')(i)];
        heel= [data{j}.('rhee_y')(i) data{j}.('rhee_z')(i)];
        toe= [data{j}.('rtoe_y')(i) data{j}.('rtoe_z')(i)];
        k2a= knee-ankle;
        h2t= toe-heel;
        temp1= atan2d(k2a(1),k2a(2));
        temp2= atan2d(h2t(1),h2t(2));
        rankle_redo(i,1)= 90+temp1-temp2;
    end
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
for k = 1:3
    for j= 1:6
        for i= 1:length(data_gc{j,k})
            if height(data_gc{j,k}{i})<10
                data_gc{j,k}(i)= [];
            end
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
                par.step_length_healthy{j,k}(i,1)= (data_gc{j,k}{i}.('rank_y')(end) - data_gc{j,k}{i}.('lank_y')(temp)) + TS(j)*(height(data_gc{j,k}{i})-temp)*10;
                
                temp= find(diff(data_gc{j,k}{i+1}.lcontact)==1);
                temp= temp(1);
                par.step_length_paretic{j,k}(i,1)= (data_gc{j,k}{i+1}.('lank_y')(temp) - data_gc{j,k}{i}.('rank_y')(end)) + TS(j)*temp*10;
                par.step_length_diff{j,k}(i,1)= par.step_length_healthy{j,k}(i,1)-par.step_length_paretic{j,k}(i,1);
            end
        elseif paretic_side{j} == 'R'
            for i= 1:L(j,k)-1
                temp= find(diff(data_gc{j,k}{i}.rcontact)==1); % Right HS
                temp= temp(1);
                par.step_length_healthy{j,k}(i,1)= (data_gc{j,k}{i}.('lank_y')(end) - data_gc{j,k}{i}.('rank_y')(temp)) + TS(j)*(height(data_gc{j,k}{i})-temp)*10;
                
                temp= find(diff(data_gc{j,k}{i+1}.rcontact)==1);
                temp= temp(1);
                par.step_length_paretic{j,k}(i,1)= (data_gc{j,k}{i+1}.('rank_y')(temp) - data_gc{j,k}{i}.('lank_y')(end)) + TS(j)*temp*10;
                par.step_length_diff{j,k}(i,1)= par.step_length_healthy{j,k}(i,1)-par.step_length_paretic{j,k}(i,1);
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
                par.anterior_healthy{j,k}(i,1)= data_gc{j,k}{i}.('lank_y')(end)-CoM(2);
                temp= find(diff(data_gc{j,k}{i}.rcontact)==1);
                temp= temp(1);
                CoM= [data_gc{j,k}{i}{temp,3:5}; data_gc{j,k}{i}{temp,6:8}; data_gc{j,k}{i}{temp,9:11}; data_gc{j,k}{i}{temp,12:14}];
                CoM= mean(CoM,1);
                par.anterior_paretic{j,k}(i,1)= data_gc{j,k}{i}.('rank_y')(temp)-CoM(2);
                par.anterior_diff{j,k}(i,1)= par.anterior_healthy{j,k}(i,1)-par.anterior_paretic{j,k}(i,1);          
            elseif paretic_side{j} == 'L'
                CoM= [data_gc{j,k}{i}{end,3:5}; data_gc{j,k}{i}{end,6:8}; data_gc{j,k}{i}{end,9:11}; data_gc{j,k}{i}{end,12:14}];
                CoM= mean(CoM,1);
                par.anterior_healthy{j,k}(i,1)= data_gc{j,k}{i}.('rank_y')(end)-CoM(2);
                temp= find(diff(data_gc{j,k}{i}.lcontact)==1);
                temp= temp(1);
                CoM= [data_gc{j,k}{i}{temp,3:5}; data_gc{j,k}{i}{temp,6:8}; data_gc{j,k}{i}{temp,9:11}; data_gc{j,k}{i}{temp,12:14}];
                CoM= mean(CoM,1);
                par.anterior_paretic{j,k}(i,1)= data_gc{j,k}{i}.('lank_y')(temp)-CoM(2);
                par.anterior_diff{j,k}(i,1)= par.anterior_healthy{j,k}(i,1)-par.anterior_paretic{j,k}(i,1); 
            end
        end
    end
end

% Stance Time
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            l_gc= height(data_gc{j,k}{i});
            if paretic_side{j} == 'R'
                par.stance_healthy{j,k}(i,1)= length(find(data_gc{j,k}{i}.lcontact==1))/l_gc*100;
                par.stance_paretic{j,k}(i,1)= length(find(data_gc{j,k}{i}.rcontact==1))/l_gc*100;
                par.stance_diff{j,k}(i,1)= par.stance_healthy{j,k}(i)-par.stance_paretic{j,k}(i);
            elseif paretic_side{j} == 'L'
                par.stance_paretic{j,k}(i,1)= length(find(data_gc{j,k}{i}.lcontact==1))/l_gc*100;
                par.stance_healthy{j,k}(i,1)= length(find(data_gc{j,k}{i}.rcontact==1))/l_gc*100;
                par.stance_diff{j,k}(i,1)= par.stance_healthy{j,k}(i)-par.stance_paretic{j,k}(i);
            end
        end
    end
end

% Push Off
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            if paretic_side{j} == 'R'
                temp= find(data_gc{j,k}{i}.rcontact == 1); % right stance
                [~,ind]= max(diff(temp));
                temp= [temp(ind+1:end); temp(1:ind)];
                temp= temp(round(length(temp)/2):end);% 2nd half of right stance
                par.ga_pushoff_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.rga(temp));
                par.force_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.force_right(temp))/weight(j)*100;
                temp= find(data_gc{j,k}{i}.lcontact == 1); % left stance
                temp= temp(round(length(temp)/2):end); % 2nd half of left stance
                par.ga_pushoff_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.lga(temp));
                par.force_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.force_left(temp))/weight(j)*100;
                par.force_diff{j,k}(i,1)= par.force_healthy{j,k}(i)-par.force_paretic{j,k}(i);
            elseif paretic_side{j} == 'L'
                temp= find(data_gc{j,k}{i}.rcontact == 1); % right stance
                temp= temp(round(length(temp)/2):end); % 2nd half of right stance
                par.ga_pushoff_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.rga(temp));
                par.force_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.force_right(temp))/weight(j)*100;
                temp= find(data_gc{j,k}{i}.lcontact == 1); % left stance
                [~,ind]= max(diff(temp));
                temp= [temp(ind+1:end); temp(1:ind)];
                temp= temp(round(length(temp)/2):end); % 2nd half of left stance
                par.ga_pushoff_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.lga(temp));
                par.force_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.force_left(temp))/weight(j)*100;
                par.force_diff{j,k}(i,1)= par.force_healthy{j,k}(i)-par.force_paretic{j,k}(i);
            end
        end
    end
end

% Muscle Acitivity
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            if paretic_side{j} == 'R'
                par.ta_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.lta);
                par.ga_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.lga);
                par.ta_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.rta);
                par.ga_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.rga);
                if ismember(j,[2 3 4 5 6])
                    par.va_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.lva);
                    par.rf_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.lrf);
                    par.bf_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.lbf);
                    par.va_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.rva);
                    par.rf_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.rrf);
                    par.bf_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.rbf);
                else
                    par.va_healthy{j,k}(i,1)= nan;
                    par.rf_healthy{j,k}(i,1)= nan;
                    par.bf_healthy{j,k}(i,1)= nan;
                    par.va_paretic{j,k}(i,1)= nan;
                    par.rf_paretic{j,k}(i,1)= nan;
                    par.bf_paretic{j,k}(i,1)= nan;
                end
            elseif paretic_side{j} == 'L'
                par.ta_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.lta);
                par.ga_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.lga);
                par.ta_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.rta);
                par.ga_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.rga);
                if ismember(j,[2 3 4 5 6])
                    par.va_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.lva);
                    par.rf_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.lrf);
                    par.bf_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.lbf);
                    par.va_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.rva);
                    par.rf_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.rrf);
                    par.bf_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.rbf);
                else
                    par.va_healthy{j,k}(i,1)= nan;
                    par.rf_healthy{j,k}(i,1)= nan;
                    par.bf_healthy{j,k}(i,1)= nan;
                    par.va_paretic{j,k}(i,1)= nan;
                    par.rf_paretic{j,k}(i,1)= nan;
                    par.bf_paretic{j,k}(i,1)= nan;
                end
            end
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
                par.knee_hs_paretic{j,k}(i,1)= data_gc{j,k}{i}.lkne_rx(temp);
                par.knee_hs_healthy{j,k}(i,1)= data_gc{j,k}{i}.rkne_rx(end);
            end
        elseif paretic_side{j} == 'R'
            for i= 1:L(j,k)-1
                temp= find(diff(data_gc{j,k}{i}.rcontact)==1); % Right HS
                temp= temp(1);
                par.knee_hs_paretic{j,k}(i,1)= data_gc{j,k}{i}.rkne_rx(temp);
                par.knee_hs_healthy{j,k}(i,1)= data_gc{j,k}{i}.lkne_rx(end);
            end
        end
    end
end

% Misc Parameters
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            if paretic_side{j} == 'R'

                par.max_knee_flexion_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.lkne_rx);
                par.max_knee_flexion_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.rkne_rx);

                % During Left Swing
                temp= find(data_gc{j,k}{i}.lcontact == 0); % left swing
                par.max_dorsiflexion_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.lank_rx(temp));
                par.max_hip_hike_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.lhip_z(temp));
                [~,ind]= min(abs(data_gc{j,k}{i}.lhip_y(temp)-data_gc{j,k}{i}.lank_y(temp)));
                par.circumduction_healthy{j,k}(i,1)= data_gc{j,k}{i}.lhip_x(temp(ind))-data_gc{j,k}{i}.lank_x(temp(ind));
                par.ta_swing_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.lta(temp));
                par.hip_velo_healthy{j,k}(i,1)= max(diff(data_gc{j,k}{i}.lhip_rx(temp))*100);
                par.knee_velo_healthy{j,k}(i,1)= min(diff(data_gc{j,k}{i}.lkne_rx(temp))*100);
                par.toe_clear_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.ltoe_z(temp));

                % During Right Swing
                temp= find(data_gc{j,k}{i}.rcontact == 0); % right swing
                par.max_dorsiflexion_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.rank_rx(temp));
                par.max_hip_hike_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.rhip_z(temp));
                [~,ind]= min(abs(data_gc{j,k}{i}.rhip_y(temp)-data_gc{j,k}{i}.rank_y(temp)));
                par.circumduction_paretic{j,k}(i,1)= data_gc{j,k}{i}.rank_x(temp(ind))-data_gc{j,k}{i}.rhip_x(temp(ind));
                par.ta_swing_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.rta(temp));
                par.hip_velo_paretic{j,k}(i,1)= max(diff(data_gc{j,k}{i}.rhip_rx(temp))*100);
                par.knee_velo_paretic{j,k}(i,1)= min(diff(data_gc{j,k}{i}.rkne_rx(temp))*100);
                par.toe_clear_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.rtoe_z(temp));

            elseif paretic_side{j} == 'L'

                par.max_knee_flexion_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.lkne_rx);
                par.max_knee_flexion_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.rkne_rx);

                % During Left Swing
                temp= find(data_gc{j,k}{i}.lcontact == 0); % left swing
                par.max_dorsiflexion_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.lank_rx(temp));
                par.max_hip_hike_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.lhip_z(temp));
                [~,ind]= min(abs(data_gc{j,k}{i}.lhip_y(temp)-data_gc{j,k}{i}.lank_y(temp)));
                par.circumduction_paretic{j,k}(i,1)= data_gc{j,k}{i}.lhip_x(temp(ind))-data_gc{j,k}{i}.lank_x(temp(ind));
                par.ta_swing_paretic{j,k}(i,1)= mean(data_gc{j,k}{i}.lta(temp));
                par.hip_velo_paretic{j,k}(i,1)= max(diff(data_gc{j,k}{i}.lhip_rx(temp))*100);
                par.knee_velo_paretic{j,k}(i,1)= min(diff(data_gc{j,k}{i}.lkne_rx(temp))*100);
                par.toe_clear_paretic{j,k}(i,1)= max(data_gc{j,k}{i}.ltoe_z(temp));

                % During Right Swing
                temp= find(data_gc{j,k}{i}.rcontact == 0); % right swing
                par.max_dorsiflexion_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.rank_rx(temp));
                par.max_hip_hike_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.rhip_z(temp));
                [~,ind]= min(abs(data_gc{j,k}{i}.rhip_y(temp)-data_gc{j,k}{i}.rank_y(temp)));
                par.circumduction_healthy{j,k}(i,1)= data_gc{j,k}{i}.rank_x(temp(ind))-data_gc{j,k}{i}.rhip_x(temp(ind));
                par.ta_swing_healthy{j,k}(i,1)= mean(data_gc{j,k}{i}.rta(temp));
                par.hip_velo_healthy{j,k}(i,1)= max(diff(data_gc{j,k}{i}.rhip_rx(temp))*100);
                par.knee_velo_healthy{j,k}(i,1)= min(diff(data_gc{j,k}{i}.rkne_rx(temp))*100);
                par.toe_clear_healthy{j,k}(i,1)= max(data_gc{j,k}{i}.rtoe_z(temp));

            end
        end
    end
end

% CoP Path (Left = Healthy, Right = Paretic)
gc_window= 30;
cop_x_offset= 452.16;
for j= 1:6
    clear cop
    for i= 1:gc_window
        % Late Baseline
        com_x= mean([data_gc_normalized{j,1}{end-gc_window+i}.lasi_x data_gc_normalized{j,1}{end-gc_window+i}.rasi_x data_gc_normalized{j,1}{end-gc_window+i}.lpsi_x data_gc_normalized{j,1}{end-gc_window+i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,1}{end-gc_window+i}.lasi_y data_gc_normalized{j,1}{end-gc_window+i}.rasi_y data_gc_normalized{j,1}{end-gc_window+i}.lpsi_y data_gc_normalized{j,1}{end-gc_window+i}.rpsi_y ],2);
        cop.x1(:,i)= data_gc_normalized{j,1}{end-gc_window+i}.cop_x - cop_x_offset;
        cop.y1(:,i)= data_gc_normalized{j,1}{end-gc_window+i}.cop_y - com_y;
        % Early Adaptation
        com_x= mean([data_gc_normalized{j,2}{i}.lasi_x data_gc_normalized{j,2}{i}.rasi_x data_gc_normalized{j,2}{i}.lpsi_x data_gc_normalized{j,2}{i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,2}{i}.lasi_y data_gc_normalized{j,2}{i}.rasi_y data_gc_normalized{j,2}{i}.lpsi_y data_gc_normalized{j,2}{i}.rpsi_y ],2);
        cop.x2(:,i)= data_gc_normalized{j,2}{i}.cop_x - cop_x_offset; 
        cop.y2(:,i)= data_gc_normalized{j,2}{i}.cop_y - com_y;
        % Late Adaptation
        com_x= mean([data_gc_normalized{j,2}{end-gc_window+i}.lasi_x data_gc_normalized{j,2}{end-gc_window+i}.rasi_x data_gc_normalized{j,2}{end-gc_window+i}.lpsi_x data_gc_normalized{j,2}{end-gc_window+i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,2}{end-gc_window+i}.lasi_y data_gc_normalized{j,2}{end-gc_window+i}.rasi_y data_gc_normalized{j,2}{end-gc_window+i}.lpsi_y data_gc_normalized{j,2}{end-gc_window+i}.rpsi_y ],2);
        cop.x3(:,i)= data_gc_normalized{j,2}{end-gc_window+i}.cop_x - cop_x_offset;
        cop.y3(:,i)= data_gc_normalized{j,2}{end-gc_window+i}.cop_y - com_y;
        % Early Observation
        com_x= mean([data_gc_normalized{j,3}{i}.lasi_x data_gc_normalized{j,3}{i}.rasi_x data_gc_normalized{j,3}{i}.lpsi_x data_gc_normalized{j,3}{i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,3}{i}.lasi_y data_gc_normalized{j,3}{i}.rasi_y data_gc_normalized{j,3}{i}.lpsi_y data_gc_normalized{j,3}{i}.rpsi_y ],2);
        cop.x4(:,i)= data_gc_normalized{j,3}{i}.cop_x - cop_x_offset;
        cop.y4(:,i)= data_gc_normalized{j,3}{i}.cop_y - com_y;
        % Late Observation
        com_x= mean([data_gc_normalized{j,3}{end-gc_window+i}.lasi_x data_gc_normalized{j,3}{end-gc_window+i}.rasi_x data_gc_normalized{j,3}{end-gc_window+i}.lpsi_x data_gc_normalized{j,3}{end-gc_window+i}.rpsi_x ],2);
        com_y= mean([data_gc_normalized{j,3}{end-gc_window+i}.lasi_y data_gc_normalized{j,3}{end-gc_window+i}.rasi_y data_gc_normalized{j,3}{end-gc_window+i}.lpsi_y data_gc_normalized{j,3}{end-gc_window+i}.rpsi_y ],2);
        cop.x5(:,i)= data_gc_normalized{j,3}{end-gc_window+i}.cop_x - cop_x_offset;
        cop.y5(:,i)= data_gc_normalized{j,3}{end-gc_window+i}.cop_y - com_y;
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
for j= 1:6
    if paretic_side{j} == 'L'
        for i= 1:5
            par.cop_path{j,i}(:,1)= -par.cop_path{j,i}(:,1) + min(par.cop_path{j,i}(:,1)) + max(par.cop_path{j,i}(:,1));
        end
    end
end

% Individual CoP Path (Left = Healthy, Right = Paretic)
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            com_y= mean([data_gc_normalized{j,k}{i}.lasi_y data_gc_normalized{j,k}{i}.rasi_y data_gc_normalized{j,k}{i}.lpsi_y data_gc_normalized{j,k}{i}.rpsi_y],2);
            cop_x= data_gc_normalized{j,k}{i}.cop_x;
            cop_y= data_gc_normalized{j,k}{i}.cop_y;
            par.cop_path_ind{j,k}{i,1}= [cop_x cop_y-com_y];
        end
    end
end
for j= 1:6
    if paretic_side{j} == 'L'
        for k= 1:3
            for i= 1:L(j,k)
                par.cop_path_ind{j,k}{i,1}(:,1)= -par.cop_path_ind{j,k}{i}(:,1) + min(par.cop_path_ind{j,k}{i}(:,1)) + max(par.cop_path_ind{j,k}{i}(:,1));
            end
        end
    end
end


% CoP Parameters
for j= 1:6
    for k= 1:5
        tempx= par.cop_path{j,k}(:,1);
        tempy= par.cop_path{j,k}(:,2);
%         close all
%         figure; hold on; axis equal
%         plot(tempx(1:50),tempy(1:50))
%         plot(tempx(51:101),tempy(51:101))
        intersection= InterX([tempx(1:50)'; tempy(1:50)'],[tempx(51:100)'; tempy(51:100)']);
        left= tempx < intersection(1);
        right= tempx > intersection(1);
        healthy_anterior= max(tempy(left));
        paretic_anterior= max(tempy(right));
        healthy_posterior= min(tempy(left));
        paretic_posterior= min(tempy(right));
        healthy_lateral= min(tempx);
        paretic_lateral= max(tempx);
        
        par.healthy_c1(j,k)= healthy_anterior-healthy_posterior;
        par.paretic_c1(j,k)= paretic_anterior-paretic_posterior;
        par.healthy_c2(j,k)= healthy_anterior-paretic_posterior;
        par.paretic_c2(j,k)= paretic_anterior-healthy_posterior;
        par.lateral_symmetry(j,k)= mean([paretic_lateral healthy_lateral]) - intersection(1);
    end
end

% Individual CoP Parameters
for k= 1:3
    for j= 1:6
        for i= 1:L(j,k)
            tempx= par.cop_path_ind{j,k}{i}(:,1);
            tempy= par.cop_path_ind{j,k}{i}(:,2);
            intersection= InterX([tempx(1:50)'; tempy(1:50)'],[tempx(51:100)'; tempy(51:100)']);
            intersection= intersection(:,1);
            left= tempx < intersection(1);
            right= tempx > intersection(1);
            healthy_anterior= max(tempy(left));
            paretic_anterior= max(tempy(right));
            healthy_posterior= min(tempy(left));
            paretic_posterior= min(tempy(right));
            healthy_lateral= min(tempx);
            paretic_lateral= max(tempx);
            par.healthy_c1_ind{j,k}(i,1)= healthy_anterior-healthy_posterior;
            par.paretic_c1_ind{j,k}(i,1)= paretic_anterior-paretic_posterior;
            par.healthy_c2_ind{j,k}(i,1)= healthy_anterior-paretic_posterior;
            par.paretic_c2_ind{j,k}(i,1)= paretic_anterior-healthy_posterior;
            par.lateral_symmetry_ind{j,k}(i,1)= mean([paretic_lateral healthy_lateral]) - intersection(1);
        end
    end
end
%% Limb Phasing
for j= 1:6
    [l_flexion_values,l_flexion_ind]= findpeaks(data{j}.lorientation,'MinPeakHeight',14);
    [r_flexion_values,r_flexion_ind]= findpeaks(data{j}.rorientation,'MinPeakHeight',14);
    [l_extension_values,l_extension_ind]= findpeaks(-data{j}.lorientation,'MinPeakHeight',-10,'MinPeakDistance',100);
    l_extension_values= -l_extension_values;
    [r_extension_values,r_extension_ind]= findpeaks(-data{j}.rorientation,'MinPeakHeight',-10,'MinPeakDistance',100);
    r_extension_values= -r_extension_values;
    if paretic_side{j} == 'L'
        for i= 1:length(r_extension_ind)-1
            h_e2= r_extension_ind(i+1);
            h_e1= r_extension_ind(i);
            if r_flexion_ind(1) < r_extension_ind(1)
                h_f= r_flexion_ind(i+1);
            else
                h_f= r_flexion_ind(i);
            end
            p_e= l_extension_ind(i+1);
            P_e= (p_e-h_e1)/(h_e2-h_e1);
            if p_e > h_f
                phase_shift{j}(i,1)= (P_e - 0.5)/2 * (1-0.5) + 0.5;
            elseif p_e <= h_f
                phase_shift{j}(i,1)= (P_e - 0)/2 * (0.5 - 0);
            end
        end
    elseif paretic_side{j} == 'R'
        for i= 1:length(l_extension_ind)-1
            h_e2= l_extension_ind(i+1);
            h_e1= l_extension_ind(i);
            h_f= l_flexion_ind(i);
            if r_extension_ind(1) > l_extension_ind(1)
                p_e= r_extension_ind(i);
            else
                p_e= r_extension_ind(i+1);
            end
            P_e= (p_e-h_e1)/(h_e2-h_e1);
            if p_e > h_f
                phase_shift{j}(i,1)= (P_e - 0.5)/2 * (1-0.5) + 0.5;
            elseif p_e <= h_f
                phase_shift{j}(i,1)= (P_e - 0)/2 * (0.5 - 0);
            end
        end 
    end
end

figure; hold on;
for j= 1:6
    subplot(2,3,j)
    plot(phase_shift{j})
end
%% Return
return
%% CoP Metrics
for k= 1:5
    for j= 1:6
        tempx= par.cop_path{j,k}(:,1);
        tempy= par.cop_path{j,k}(:,2);
        close all
        figure; hold on
        axis([-200 200 -200 200])
        plot(tempx(1:50),tempy(1:50))
        plot(tempx(51:101),tempy(51:101))
        inter= InterX([tempx(1:50)'; tempy(1:50)'],[tempx(51:101)'; tempy(51:101)']);
        x_mean= mean([min(tempx) max(tempx)]);
        [~,ind]= min(abs(x_mean-(inter(1,:))));
        temp= inter(:,ind)';
        par.cop_cross{j,k}= temp;
        par.cop_cross_to_center(j,k)= hypot(temp(1),temp(2));
    end
end
%% Last Thing Before Plots
return
%% Plot
close all
clc

% Setup Plot
j= 1; % Subject Number
figure; set(gcf,'color','w','Position',[1 74 1440 723]); hold on;
subplot(3,4,1); hold on;
if paretic_side{j} == 'L'
    main_title= sgtitle(['Subject ' num2str(j)  ' | \color[rgb]{0 .4471 .7647}Left Side Paretic \color{black}| ' '\color[rgb]{.8 .2078 .1451}Right Side Perturbed \color{black}| '  num2str(dur(j)/60) ' minutes | ' num2str(TS(j)) ' m/s']);
else
    main_title= sgtitle(['Subject ' num2str(j)  ' | \color[rgb]{.8 .2078 .1451}Right Side Paretic \color{black}| ' '\color[rgb]{0 .4471 .7647}Left Side Perturbed \color{black}| '  num2str(dur(j)/60) ' minutes | ' num2str(TS(j)) ' m/s']);
end
set(main_title,'FontWeight','bold','FontSize',22);

title_input= 'Frozen Step Length';
healthy_data_input= par.left_step_length(j,:);
paretic_data_input= par.right_step_length(j,:);
diff_data_input= par.step_length_diff(j,:);
ylabel_input= 'Ankle to Ankle (mm)';
plot_number_input= 1;
double_axis_plot(healthy_data_input,paretic_data_input,diff_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Stride Length';
data_input= par.stride_length(j,:);
ylabel_input= 'Ankle to Ankle (mm)';
plot_number_input= 2;
single_plot(data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Anterior Step Length';
healthy_data_input= par.left_anterior(j,:);
paretic_data_input= par.right_anterior(j,:);
diff_data_input= par.anterior_diff(j,:);
ylabel_input= 'CoM to Ankle (mm)';
plot_number_input= 3;
double_axis_plot(healthy_data_input,paretic_data_input,diff_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Stance Time';
healthy_data_input= par.left_stance(j,:);
paretic_data_input= par.right_stance(j,:);
diff_data_input= par.stance_diff(j,:);
ylabel_input= '% Gait Cycle';
plot_number_input= 4;
double_axis_plot(healthy_data_input,paretic_data_input,diff_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Max Knee Flexion during Swing';
healthy_data_input= par.max_left_knee_flexion(j,:);
paretic_data_input= par.max_right_knee_flexion(j,:);
ylabel_input= 'Flexion Angle (degrees)';
plot_number_input= 5;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Max Dorsiflexion during Swing';
healthy_data_input= par.max_left_dorsiflexion(j,:);
paretic_data_input= par.max_right_dorsiflexion(j,:);
ylabel_input= 'Flexion Angle (degrees)';
plot_number_input= 6;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Knee Angle at Heel Strike';
healthy_data_input= par.left_knee_hs(j,:);
paretic_data_input= par.right_knee_hs(j,:);
ylabel_input= 'Flexion Angle (degrees)';
plot_number_input= 7;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Center of Pressure';
data_input= par.cop_path(j,:);
plot_number_input= 8;
cop_plot(data_input,title_input,plot_number_input)

title_input= 'Hip Hiking';
healthy_data_input= par.max_left_hip_hike(j,:);
paretic_data_input= par.max_right_hip_hike(j,:);
ylabel_input= 'Hip Height (mm)';
plot_number_input= 9;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Hip Circumduction';
healthy_data_input= par.left_circumduction(j,:);
paretic_data_input= par.right_circumduction(j,:);
ylabel_input= 'Distance from Hip (mm)';
plot_number_input= 10;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Toe Clearance';
healthy_data_input= par.left_toe_clear(j,:);
paretic_data_input= par.right_toe_clear(j,:);
ylabel_input= 'Toe Height (mm)';
plot_number_input= 11;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Stance Force';
healthy_data_input= par.left_force(j,:);
paretic_data_input= par.right_force(j,:);
ylabel_input= 'Avg Force (units)';
plot_number_input= 12;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'GA Activation at Push Off';
healthy_data_input= par.left_ga(j,:);
paretic_data_input= par.right_ga(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 13;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'TA Activation during Swing';
healthy_data_input= par.left_ta_activation(j,:);
paretic_data_input= par.right_ta_activation(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 14;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Hip Velocity during Swing';
healthy_data_input= par.left_hip_velo(j,:);
paretic_data_input= par.right_hip_velo(j,:);
ylabel_input= 'Angular Velocity (deg/s)';
plot_number_input= 15;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)

title_input= 'Knee Velocity during Swing';
healthy_data_input= par.left_knee_velo(j,:);
paretic_data_input= par.right_knee_velo(j,:);
ylabel_input= 'Angular Velocity (deg/s)';
plot_number_input= 16;
normal_plot(healthy_data_input,paretic_data_input,trans(j,:),Len(j),title_input,ylabel_input,plot_number_input)
%% BoxPlot
close all
clc

% Setup Plot
j= 1; % Subject Number
figure; set(gcf,'color','w','Position',[-1452 86 1440 723]); hold on;
subplot(3,4,1); hold on;
main_title= sgtitle(['Subject ' num2str(j)  ' | \color[rgb]{0 .4471 .7647} Perturbed Side \color{black}| ' '\color[rgb]{.8 .2078 .1451} Paretic Side \color{black}| '  num2str(dur(j)/60) ' minutes | ' num2str(TS(j)) ' m/s']);
set(main_title,'FontWeight','bold','FontSize',22);

title_input= 'Frozen Step Length';
healthy_data_input= par.step_length_healthy(j,:);
paretic_data_input= par.step_length_paretic(j,:);
ylabel_input= 'Ankle to Ankle (mm)';
plot_number_input= 1;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Anterior Step Length';
healthy_data_input= par.anterior_healthy(j,:);
paretic_data_input= par.anterior_paretic(j,:);
ylabel_input= 'CoM to Ankle (mm)';
plot_number_input= 2;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Stance Time';
only_data_input= par.stance_diff(j,:);
ylabel_input= '% Gait Cycle';
plot_number_input= 3;
single_box_whisker_plot(only_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Stance Force';
healthy_data_input= par.force_healthy(j,:);
paretic_data_input= par.force_paretic(j,:);
ylabel_input= 'Avg Force (units)';
plot_number_input= 4;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average TA Activation';
healthy_data_input= par.ta_healthy(j,:);
paretic_data_input= par.ta_paretic(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 5;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average GA Activation';
healthy_data_input= par.ga_healthy(j,:);
paretic_data_input= par.ga_paretic(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 6;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

if j==2 | j==3 || j==4 || j==5 || j==6
    title_input= 'Average VA Activation';
    healthy_data_input= par.va_healthy(j,:);
    paretic_data_input= par.va_paretic(j,:);
    ylabel_input= 'Activation Level (%)';
    plot_number_input= 7;
    box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)
    
    title_input= 'Average RF Activation';
    healthy_data_input= par.rf_healthy(j,:);
    paretic_data_input= par.rf_paretic(j,:);
    ylabel_input= 'Activation Level (%)';
    plot_number_input= 8;
    box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)
    
    title_input= 'Average BF Activation';
    healthy_data_input= par.bf_healthy(j,:);
    paretic_data_input= par.bf_paretic(j,:);
    ylabel_input= 'Activation Level (%)';
    plot_number_input= 9;
    box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)
end

title_input= 'TA Activation during Swing';
healthy_data_input= par.ta_swing_healthy(j,:);
paretic_data_input= par.ta_swing_paretic(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 10;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'GA Activation at Push Off';
healthy_data_input= par.ga_pushoff_healthy(j,:);
paretic_data_input= par.ga_pushoff_paretic(j,:);
ylabel_input= 'Activation Level (%)';
plot_number_input= 11;
box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Center of Pressure';
data_input= par.cop_path(j,:);
plot_number_input= 12;
cop_plot(data_input,title_input,plot_number_input)

saveas(gcf,['Subject' num2str(j) '.png'])
%% Average BoxPlot
close all
clc

% Setup Plot
figure; set(gcf,'color','w','Position',[-1452 86 1440 723]); hold on;
subplot(3,4,1); hold on;
main_title= sgtitle('Average Data | \color[rgb]{0 .4471 .7647} Perturbed Side \color{black}| \color[rgb]{.8 .2078 .1451} Paretic Side \color{black}');
set(main_title,'FontWeight','bold','FontSize',22);

title_input= 'Frozen Step Length';
healthy_data_input= par.step_length_healthy;
paretic_data_input= par.step_length_paretic;
ylabel_input= 'Ankle to Ankle (mm)';
plot_number_input= 1;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Anterior Step Length';
healthy_data_input= par.anterior_healthy;
paretic_data_input= par.anterior_paretic;
ylabel_input= 'CoM to Ankle (mm)';
plot_number_input= 2;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Stance Time';
only_data_input= par.stance_diff;
ylabel_input= '% Gait Cycle';
plot_number_input= 3;
avg_single_box_whisker_plot(only_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Stance Force';
healthy_data_input= par.force_healthy;
paretic_data_input= par.force_paretic;
ylabel_input= 'Avg Force (units)';
plot_number_input= 4;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average TA Activation';
healthy_data_input= par.ta_healthy;
paretic_data_input= par.ta_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 5;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average GA Activation';
healthy_data_input= par.ga_healthy;
paretic_data_input= par.ga_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 6;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average VA Activation';
healthy_data_input= par.va_healthy;
paretic_data_input= par.va_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 7;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average RF Activation';
healthy_data_input= par.rf_healthy;
paretic_data_input= par.rf_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 8;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'Average BF Activation';
healthy_data_input= par.bf_healthy;
paretic_data_input= par.bf_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 9;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'TA Activation during Swing';
healthy_data_input= par.ta_swing_healthy;
paretic_data_input= par.ta_swing_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 10;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

title_input= 'GA Activation at Push Off';
healthy_data_input= par.ga_pushoff_healthy;
paretic_data_input= par.ga_pushoff_paretic;
ylabel_input= 'Activation Level (%)';
plot_number_input= 11;
avg_box_whisker_plot(healthy_data_input,paretic_data_input,title_input,ylabel_input,plot_number_input)

% title_input= 'Center of Pressure';
% data_input= par.cop_path(j,:);
% plot_number_input= 12;
% cop_plot(data_input,title_input,plot_number_input)
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
%% **** FINAL FIGURES ****
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 30;
%% RESULTS AVERAGE
%% Step Length & Anterior Average - Box & Whisker
% Step Length
plot_title= 'Step Length';
healthy_data_input= par.step_length_healthy;
paretic_data_input= par.step_length_paretic;
plot_ylabel= 'Ankle to Ankle (mm)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

figure; set(gcf,'color','w','position',[325 425 959 345]); hold on;
subplot(1,2,1); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline         Adaptation              Observation      ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color','interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% Anterior Step Length
plot_title= 'Anterior Step Length';
healthy_data_input= par.anterior_healthy;
paretic_data_input= par.anterior_paretic;
plot_ylabel= 'CoM to Ankle (mm)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(1,2,2); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline        Adaptation            Observation     ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

saveas(gcf,'Fig_StepLength.png')
%% Stance Symmetry - Box & Whisker
plot_title= 'Stance Phase Asymmetry';
data_input= par.stance_diff;
plot_ylabel= 'Gait Cycle (%)';

poi= data_input;
for j= 1:6
    temp_a(:,j)= poi{j,1}(end-SPAN+1:end);
    temp_b1(:,j)= poi{j,2}(1:SPAN);
    temp_b2(:,j)= poi{j,2}(end-SPAN+1:end);
    temp_c1(:,j)= poi{j,3}(1:SPAN);
    temp_c2(:,j)= poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(temp_a,2);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= mean(temp_b1,2); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= mean(temp_b2,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(temp_c1,2);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= mean(temp_c2,2);

poi= [box_data(:,1) box_data(:,3) box_data(:,5) box_data(:,7) box_data(:,9)];

colors = purple(1:3);

figure; set(gcf,'color','w'); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline                  Adaptation                          Observation          ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(poi);
if signif(1) == 1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(1) == -1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
end
if signif(2) == -1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(2) == 1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
end

saveas(gcf,'Fig_StanceSymmetry.png')
%% Push-off Asymmetry
plot_title= 'Vertical Push-Off GRF Asymmetry';
data_input= par.force_diff;
plot_ylabel= '% Body Weight';

poi= data_input;
for j= 1:6
    temp_a(:,j)= poi{j,1}(end-SPAN+1:end);
    temp_b1(:,j)= poi{j,2}(1:SPAN);
    temp_b2(:,j)= poi{j,2}(end-SPAN+1:end);
    temp_c1(:,j)= poi{j,3}(1:SPAN);
    temp_c2(:,j)= poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(temp_a,2);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= mean(temp_b1,2); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= mean(temp_b2,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(temp_c1,2);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= mean(temp_c2,2);

poi= [box_data(:,1) box_data(:,3) box_data(:,5) box_data(:,7) box_data(:,9)];

colors = purple(1:3);

figure; set(gcf,'color','w'); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline                  Adaptation                          Observation          ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(poi);
if signif(1) == 1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(1) == -1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
end
if signif(2) == -1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(2) == 1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
end

saveas(gcf,'Fig_PushOffSymmetry.png')
%% Step Length Symmetry
plot_title= 'Step Length Asymmetry';
data_input= par.step_length_diff;
plot_ylabel= 'Ankle to Ankle (mm)';

poi= data_input;
for j= 1:6
    temp_a(:,j)= poi{j,1}(end-SPAN+1:end);
    temp_b1(:,j)= poi{j,2}(1:SPAN);
    temp_b2(:,j)= poi{j,2}(end-SPAN+1:end);
    temp_c1(:,j)= poi{j,3}(1:SPAN);
    temp_c2(:,j)= poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(temp_a,2);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= mean(temp_b1,2); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= mean(temp_b2,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(temp_c1,2);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= mean(temp_c2,2);

poi= [box_data(:,1) box_data(:,3) box_data(:,5) box_data(:,7) box_data(:,9)];

colors = purple(1:3);

figure; set(gcf,'color','w','position',[281 499 1116 420]); hold on;
subplot(1,2,1); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline                  Adaptation                          Observation          ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(poi);
if signif(1) == 1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(1) == -1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
end
if signif(2) == -1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(2) == 1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
end

plot_title= 'Anterior Step Length Asymmetry';
data_input= par.anterior_diff;
plot_ylabel= 'CoM to Ankle (mm)';

poi= data_input;
for j= 1:6
    temp_a(:,j)= poi{j,1}(end-SPAN+1:end);
    temp_b1(:,j)= poi{j,2}(1:SPAN);
    temp_b2(:,j)= poi{j,2}(end-SPAN+1:end);
    temp_c1(:,j)= poi{j,3}(1:SPAN);
    temp_c2(:,j)= poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(temp_a,2);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= mean(temp_b1,2); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= mean(temp_b2,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(temp_c1,2);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= mean(temp_c2,2);

poi= [box_data(:,1) box_data(:,3) box_data(:,5) box_data(:,7) box_data(:,9)];

colors = purple(1:3);

subplot(1,2,2); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline                  Adaptation                          Observation          ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(poi);
if signif(1) == 1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(1) == -1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
end
if signif(2) == -1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(2) == 1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
end

saveas(gcf,'Fig_StepLengthSymmetry.png')
%% Anterior Step Length Symmetry
plot_title= 'Anterior Step Length Symmetry';
data_input= par.anterior_diff;
plot_ylabel= 'CoM to Ankle (mm)';

poi= data_input;
for j= 1:6
    temp_a(:,j)= poi{j,1}(end-SPAN+1:end);
    temp_b1(:,j)= poi{j,2}(1:SPAN);
    temp_b2(:,j)= poi{j,2}(end-SPAN+1:end);
    temp_c1(:,j)= poi{j,3}(1:SPAN);
    temp_c2(:,j)= poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(temp_a,2);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= mean(temp_b1,2); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= mean(temp_b2,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(temp_c1,2);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= mean(temp_c2,2);

poi= [box_data(:,1) box_data(:,3) box_data(:,5) box_data(:,7) box_data(:,9)];

colors = purple(1:3);

figure; set(gcf,'color','w'); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline                  Adaptation                          Observation          ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(poi);
if signif(1) == 1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(1) == -1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
end
if signif(2) == -1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(2) == 1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
end

saveas(gcf,'Fig_AnteriorStepLengthSymmetry.png')
%% Muscle Activity - Box & Whisker
% TA
plot_title= 'Tibialis Anterior';
healthy_data_input= par.ta_healthy;
paretic_data_input= par.ta_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

figure; set(gcf,'color','w','position',[-62 410 1924 420]); hold on;
subplot(1,5,1); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% GA
plot_title= 'Gastrocnemius';
healthy_data_input= par.ga_healthy;
paretic_data_input= par.ga_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(1,5,2); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% VA
plot_title= 'Vastus Medialis';
healthy_data_input= par.va_healthy;
paretic_data_input= par.va_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

h_temp_a(:,1)= [];
h_temp_b1(:,1)= [];
h_temp_b2(:,1)= [];
h_temp_c1(:,1)= [];
h_temp_c2(:,1)= [];
p_temp_a(:,1)= [];
p_temp_b1(:,1)= [];
p_temp_b2(:,1)= [];
p_temp_c1(:,1)= [];
p_temp_c2(:,1)= [];


box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(1,5,3); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% RF
plot_title= 'Rectus Femoris';
healthy_data_input= par.rf_healthy;
paretic_data_input= par.rf_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

h_temp_a(:,1)= [];
h_temp_b1(:,1)= [];
h_temp_b2(:,1)= [];
h_temp_c1(:,1)= [];
h_temp_c2(:,1)= [];
p_temp_a(:,1)= [];
p_temp_b1(:,1)= [];
p_temp_b2(:,1)= [];
p_temp_c1(:,1)= [];
p_temp_c2(:,1)= [];


box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(1,5,4); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% BF
plot_title= 'Biceps Femoris';
healthy_data_input= par.bf_healthy;
paretic_data_input= par.bf_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

h_temp_a(:,1)= [];
h_temp_b1(:,1)= [];
h_temp_b2(:,1)= [];
h_temp_c1(:,1)= [];
h_temp_c2(:,1)= [];
p_temp_a(:,1)= [];
p_temp_b1(:,1)= [];
p_temp_b2(:,1)= [];
p_temp_c1(:,1)= [];
p_temp_c2(:,1)= [];


box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(1,5,5); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

saveas(gcf,'Fig_MuscleActivity.png')
%% Push Off Force Average - Box & Whisker
% Step Length
plot_title= 'Vertical Push-Off GRF';
healthy_data_input= par.force_healthy;
paretic_data_input= par.force_paretic;
plot_ylabel= '% Body Weight';

h_poi= healthy_data_input;
for j= 1:6
    h_temp_a(:,j)= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1(:,j)= h_poi{j,2}(1:SPAN);
    h_temp_b2(:,j)= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1(:,j)= h_poi{j,3}(1:SPAN);
    h_temp_c2(:,j)= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 1:6
    p_temp_a(:,j)= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1(:,j)= p_poi{j,2}(1:SPAN);
    p_temp_b2(:,j)= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1(:,j)= p_poi{j,3}(1:SPAN);
    p_temp_c2(:,j)= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= mean(h_temp_a,2);
box_data(:,2)= mean(p_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(h_temp_b1,2);
box_data(:,5)= mean(p_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(h_temp_b2,2);
box_data(:,8)= mean(p_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(h_temp_c1,2);
box_data(:,11)= mean(p_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(h_temp_c2,2);
box_data(:,14)= mean(p_temp_c2,2);

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

figure; set(gcf,'color','w'); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel(' Baseline                  Adaptation                         Observation          ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end
saveas(gcf,'Fig_PushOffForce.png')
%% CoP Path
plot_title= 'CoP Path';
data_input= par.cop_path;
plot_xlabel= 'Medial/Lateral Position (mm)';
plot_ylabel= 'Anterior/Posterior Position (mm)';

poi= data_input;
for j= 1:6
    temp_ax(:,j)= poi{j,1}(:,1);
    temp_ay(:,j)= poi{j,1}(:,2);
    temp_c1x(:,j)= poi{j,4}(:,1);
    temp_c1y(:,j)= poi{j,4}(:,2);
    temp_c2x(:,j)= poi{j,5}(:,1);
    temp_c2y(:,j)= poi{j,5}(:,2);
end

plot_data= [];
plot_data(:,1)= mean(temp_ax,2);
plot_data(:,2)= mean(temp_ay,2);
plot_data(:,3)= mean(temp_c1x,2);
plot_data(:,4)= mean(temp_c1y,2);
plot_data(:,5)= mean(temp_c2x,2);
plot_data(:,6)= mean(temp_c2y,2);

for i= 1:6
    plot_data(:,i)= smooth(plot_data(:,i));
end

colororder('default');

figure; set(gcf,'color','w'); hold on;
plot(plot_data(:,1),plot_data(:,2),'LineWidth',2)
plot(plot_data(:,3),plot_data(:,4),'LineWidth',2)
plot(plot_data(:,5),plot_data(:,6),'LineWidth',2)
axis equal
axis manual
plot([-1000 1000],[0 0],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([0 0],[-1000 1000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% legend('Baseline','Early Adapt','Late Adapt','Early Obs','Late Obs','Position',[0.90827 0.23599 0.076389 0.11203])
legend('Baseline','Early Obs','Late Obs','Position',[0.5 0.15 0.186 0.126])
title(plot_title)
xlabel(plot_xlabel)
ylabel(plot_ylabel)

saveas(gcf,'Fig_CoP.png')
%% Stance Symmetry Average - Box & Whisker
% Average Muscle Activty for all 5 muscles - Box & Whisker
% Push Off Force Average - Box & Whisker
% 2 Representative (Good & Bad) CoP Path Plot
%% CoP Plot (Subject 3)
figure; set(gcf,'color','w','position',[680 558 1097 420]); hold on;
title_input= 'Center of Pressure (Subject 3)';
data_input= par.cop_path(3,:);
plot_number_input= 2;
double_cop_plot(data_input,title_input,plot_number_input)

par1x= data_input{5}(:,1);
par1y= data_input{5}(:,2);

subplot(1,2,1); hold on;
plot(par1x,par1y,'LineWidth',2)
box on
axis([-150 150 -150 150])
% axis equal
% axis manual
% plot([-1000 1000],[0 0],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% plot([0 0],[-1000 1000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% legend('Baseline','Early Adapt','Late Adapt','Early Obs','Late Obs','Position',[0.90827 0.23599 0.076389 0.11203])
% legend('Baseline','Early Obs','Late Obs','Position',[0.79879 0.14558 0.094804 0.12619])
title('Example Center of Pressure Path')
xlabel('Medial/Lateral Position (mm)')
ylabel('Anterior/Posterior Position (mm)')
saveas(gcf,'Fig_CoPPath.png')
%% Muscle Activitly (Subject 2)
% TA
plot_title= 'Tibialis Anterior (TA)';
healthy_data_input= par.ta_healthy;
paretic_data_input= par.ta_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 2
    h_temp_a= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1= h_poi{j,2}(1:SPAN);
    h_temp_b2= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1= h_poi{j,3}(1:SPAN);
    h_temp_c2= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 2
    p_temp_a= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1= p_poi{j,2}(1:SPAN);
    p_temp_b2= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1= p_poi{j,3}(1:SPAN);
    p_temp_c2= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= h_temp_a;
box_data(:,2)= p_temp_a;
box_data(1:SPAN,3)= nan; 
box_data(:,4)= h_temp_b1;
box_data(:,5)= p_temp_b1;
box_data(1:SPAN,6)= nan; 
box_data(:,7)= h_temp_b2;
box_data(:,8)= p_temp_b2;
box_data(1:SPAN,9)= nan; 
box_data(:,10)= h_temp_c1;
box_data(:,11)= p_temp_c1;
box_data(1:SPAN,12)= nan; 
box_data(:,13)= h_temp_c2;
box_data(:,14)= p_temp_c2;

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

figure; set(gcf,'color','w','position',[279 118 1426 788]); hold on;
subplot(2,3,1); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% GA
plot_title= 'Gastrocnemius (GA)';
healthy_data_input= par.ga_healthy;
paretic_data_input= par.ga_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 2
    h_temp_a= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1= h_poi{j,2}(1:SPAN);
    h_temp_b2= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1= h_poi{j,3}(1:SPAN);
    h_temp_c2= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 2
    p_temp_a= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1= p_poi{j,2}(1:SPAN);
    p_temp_b2= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1= p_poi{j,3}(1:SPAN);
    p_temp_c2= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= h_temp_a;
box_data(:,2)= p_temp_a;
box_data(1:SPAN,3)= nan; 
box_data(:,4)= h_temp_b1;
box_data(:,5)= p_temp_b1;
box_data(1:SPAN,6)= nan; 
box_data(:,7)= h_temp_b2;
box_data(:,8)= p_temp_b2;
box_data(1:SPAN,9)= nan; 
box_data(:,10)= h_temp_c1;
box_data(:,11)= p_temp_c1;
box_data(1:SPAN,12)= nan; 
box_data(:,13)= h_temp_c2;
box_data(:,14)= p_temp_c2;

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(2,3,2); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% VA
plot_title= 'Vastus Medialis (VA)';
healthy_data_input= par.va_healthy;
paretic_data_input= par.va_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 2
    h_temp_a= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1= h_poi{j,2}(1:SPAN);
    h_temp_b2= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1= h_poi{j,3}(1:SPAN);
    h_temp_c2= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 2
    p_temp_a= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1= p_poi{j,2}(1:SPAN);
    p_temp_b2= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1= p_poi{j,3}(1:SPAN);
    p_temp_c2= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= h_temp_a;
box_data(:,2)= p_temp_a;
box_data(1:SPAN,3)= nan; 
box_data(:,4)= h_temp_b1;
box_data(:,5)= p_temp_b1;
box_data(1:SPAN,6)= nan; 
box_data(:,7)= h_temp_b2;
box_data(:,8)= p_temp_b2;
box_data(1:SPAN,9)= nan; 
box_data(:,10)= h_temp_c1;
box_data(:,11)= p_temp_c1;
box_data(1:SPAN,12)= nan; 
box_data(:,13)= h_temp_c2;
box_data(:,14)= p_temp_c2;

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(2,3,3); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% RF
plot_title= 'Rectus Femoris (RF)';
healthy_data_input= par.rf_healthy;
paretic_data_input= par.rf_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 2
    h_temp_a= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1= h_poi{j,2}(1:SPAN);
    h_temp_b2= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1= h_poi{j,3}(1:SPAN);
    h_temp_c2= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 2
    p_temp_a= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1= p_poi{j,2}(1:SPAN);
    p_temp_b2= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1= p_poi{j,3}(1:SPAN);
    p_temp_c2= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= h_temp_a;
box_data(:,2)= p_temp_a;
box_data(1:SPAN,3)= nan; 
box_data(:,4)= h_temp_b1;
box_data(:,5)= p_temp_b1;
box_data(1:SPAN,6)= nan; 
box_data(:,7)= h_temp_b2;
box_data(:,8)= p_temp_b2;
box_data(1:SPAN,9)= nan; 
box_data(:,10)= h_temp_c1;
box_data(:,11)= p_temp_c1;
box_data(1:SPAN,12)= nan; 
box_data(:,13)= h_temp_c2;
box_data(:,14)= p_temp_c2;

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(2,3,4); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

% BF
plot_title= 'Biceps Femoris (BF)';
healthy_data_input= par.bf_healthy;
paretic_data_input= par.bf_paretic;
plot_ylabel= 'Muscle Activity (%)';

h_poi= healthy_data_input;
for j= 2
    h_temp_a= h_poi{j,1}(end-SPAN+1:end);
    h_temp_b1= h_poi{j,2}(1:SPAN);
    h_temp_b2= h_poi{j,2}(end-SPAN+1:end);
    h_temp_c1= h_poi{j,3}(1:SPAN);
    h_temp_c2= h_poi{j,3}(end-SPAN+1:end);
end

p_poi= paretic_data_input;
for j= 2
    p_temp_a= p_poi{j,1}(end-SPAN+1:end);
    p_temp_b1= p_poi{j,2}(1:SPAN);
    p_temp_b2= p_poi{j,2}(end-SPAN+1:end);
    p_temp_c1= p_poi{j,3}(1:SPAN);
    p_temp_c2= p_poi{j,3}(end-SPAN+1:end);
end

box_data= [];
box_data(:,1)= h_temp_a;
box_data(:,2)= p_temp_a;
box_data(1:SPAN,3)= nan; 
box_data(:,4)= h_temp_b1;
box_data(:,5)= p_temp_b1;
box_data(1:SPAN,6)= nan; 
box_data(:,7)= h_temp_b2;
box_data(:,8)= p_temp_b2;
box_data(1:SPAN,9)= nan; 
box_data(:,10)= h_temp_c1;
box_data(:,11)= p_temp_c1;
box_data(1:SPAN,12)= nan; 
box_data(:,13)= h_temp_c2;
box_data(:,14)= p_temp_c2;

healthy_poi= [box_data(:,1) box_data(:,4) box_data(:,7) box_data(:,10) box_data(:,13)];
paretic_poi= [box_data(:,2) box_data(:,5) box_data(:,8) box_data(:,11) box_data(:,14)];

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(2,3,5); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
% ylabel(plot_ylabel)
xlabel('Baseline     Adaptation          Observation    ')

inc= 0.05;
ax= axis;
axis([ax(1) ax(2) ax(3) (ax(4)-ax(3))*(1+inc*6)+ax(3)])
y1= (ax(4)-ax(3))*(1+inc*1)+ax(3); 
y2= (ax(4)-ax(3))*(1+inc*2)+ax(3);
y3= (ax(4)-ax(3))*(1+inc*3)+ax(3);
y4= (ax(4)-ax(3))*(1+inc*4)+ax(3);
y5= (ax(4)-ax(3))*(1+inc*5)+ax(3);

signif= significance_test_average(healthy_poi);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test_average(paretic_poi);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

saveas(gcf,'Fig_MuscleActivity2.png')
%% Extract Data for Supplementary Materials
poi= par.lateral_symmetry_ind;
poi_extract= [];
for j= 1:6
    significance(j,:)= significance_test(poi(j,:));
end
significance

for j= 1:6
    poi_extract{j,1}= num2str(round(mean(poi{j,1}(end-SPAN+1:end)),3,'significant'));
    poi_extract{j,2}= num2str(round(mean(poi{j,3}(1:SPAN)),3,'significant'));
    if significance(j,1) ~= 0
        poi_extract{j,2}= ['\textbf{' poi_extract{j,2} '}'];
    end
    poi_extract{j,3}= num2str(round(mean(poi{j,3}(end-SPAN+1:end)),3,'significant'));
    if significance(j,2) ~= 0
        poi_extract{j,3}= ['\textbf{' poi_extract{j,3} '}'];
    end
end
for_copy{1,1}= ...
    [' & ' poi_extract{1,1}...
    ' & ' poi_extract{1,2}...
    ' & ' poi_extract{1,3}...
    ' & ' poi_extract{2,1}...
    ' & ' poi_extract{2,2}...
    ' & ' poi_extract{2,3}...
    ' & ' poi_extract{3,1}...
    ' & ' poi_extract{3,2}...
    ' & ' poi_extract{3,3} ' \\'];
for_copy{2,1}= ...
    [' & ' poi_extract{4,1}...
    ' & ' poi_extract{4,2}...
    ' & ' poi_extract{4,3}...
    ' & ' poi_extract{5,1}...
    ' & ' poi_extract{5,2}...
    ' & ' poi_extract{5,3}...
    ' & ' poi_extract{6,1}...
    ' & ' poi_extract{6,2}...
    ' & ' poi_extract{6,3} ' \\']
%% Functions
function avg_single_box_whisker_plot(parameter_of_interest,plot_title,plot_ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 30;

poi= parameter_of_interest;
for j= 1:6
    temp_a(:,j)= poi{j,1}(end-SPAN+1:end);
    temp_b1(:,j)= poi{j,2}(1:SPAN);
    temp_b2(:,j)= poi{j,2}(end-SPAN+1:end);
    temp_c1(:,j)= poi{j,3}(1:SPAN);
    temp_c2(:,j)= poi{j,3}(end-SPAN+1:end);
end

box_data(:,1)= mean(temp_a,2);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= mean(temp_b1,2); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= mean(temp_b2,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(temp_c1,2);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= mean(temp_c2,2);

colors = purple(1:3);

subplot(3,4,subplot_number); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

% signif= significance_test(only_data);
% ax= axis;
% y1= ax(4)*1.05;
% y2= ax(4)*1.25;
% axis([ax(1) ax(2) ax(3) ax(4)*1.5])
% if signif(1) == 1
%     plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
%     text(mean([1 7]),mean([y1 y2]),'??','fontsize',8,'color',purple)
% elseif signif(1) == -1
%     plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
%     text(mean([1 7]),mean([y1 y2]),'??','fontsize',8,'color',purple)
% end
% if signif(2) == -1
%     plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
%     text(mean([1 9]),mean([y2 ax(4)*1.5]),'??','fontsize',8,'color',purple)
% elseif signif(2) == 1
%     plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
%     text(mean([1 9]),mean([y2 ax(4)*1.5]),'??','fontsize',8,'color',purple)
% end

xticks([])
title(plot_title)
ylabel(plot_ylabel)

end

function avg_box_whisker_plot(left_parameter_of_interest,right_parameter_of_interest,plot_title,plot_ylabel,subplot_number)
% box_whisker_plot(left_data,right_data,plot_title,plot_ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 30;

l_poi= left_parameter_of_interest;
for j= 1:6
    l_temp_a(:,j)= l_poi{j,1}(end-SPAN+1:end);
    l_temp_b1(:,j)= l_poi{j,2}(1:SPAN);
    l_temp_b2(:,j)= l_poi{j,2}(end-SPAN+1:end);
    l_temp_c1(:,j)= l_poi{j,3}(1:SPAN);
    l_temp_c2(:,j)= l_poi{j,3}(end-SPAN+1:end);
end

r_poi= right_parameter_of_interest;
for j= 1:6
    r_temp_a(:,j)= r_poi{j,1}(end-SPAN+1:end);
    r_temp_b1(:,j)= r_poi{j,2}(1:SPAN);
    r_temp_b2(:,j)= r_poi{j,2}(end-SPAN+1:end);
    r_temp_c1(:,j)= r_poi{j,3}(1:SPAN);
    r_temp_c2(:,j)= r_poi{j,3}(end-SPAN+1:end);
end
 
box_data(:,1)= mean(l_temp_a,2);
box_data(:,2)= mean(r_temp_a,2);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= mean(l_temp_b1,2);
box_data(:,5)= mean(r_temp_b1,2);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= mean(l_temp_b2,2);
box_data(:,8)= mean(r_temp_b2,2);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= mean(l_temp_c1,2);
box_data(:,11)= mean(r_temp_c1,2);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= mean(l_temp_c2,2);
box_data(:,14)= mean(r_temp_c2,2);

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(3,4,subplot_number); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)


% ax= axis;
% axis([ax(1) ax(2) ax(3) ax(4)*1.5])
% y1= ax(4)*1.05;
% y2= ax(4)*1.15;
% y3= ax(4)*1.25;
% y4= ax(4)*1.35;
% y5= ax(4)*1.45;
% 
% signif= significance_test(left_data);
% if signif(1) == 1
%     plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
%     text(mean([1 10]),mean([y1 y2]),'??','fontsize',8,'color',blue)
% elseif signif(1) == -1
%     plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
%     text(mean([1 10]),mean([y1 y2]),'??','fontsize',8,'color',blue)
% end
% if signif(2) == -1
%     plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
%     text(mean([1 13]),mean([y2 y3]),'??','fontsize',8,'color',blue)
% elseif signif(2) == 1
%     plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
%     text(mean([1 13]),mean([y2 y3]),'??','fontsize',8,'color',blue)
% end
% 
% signif= significance_test(right_data);
% if signif(1) == 1
%     plot([2 11],[y3 y3],'color',red,'LineWidth',1)
%     text(mean([2 11]),mean([y3 y4]),'??','fontsize',8,'color',red)
% elseif signif(1) == -1
%     plot([2 11],[y3 y3],'color',red,'LineWidth',1)
%     text(mean([2 11]),mean([y3 y4]),'??','fontsize',8,'color',red)
% end
% if signif(2) == -1
%     plot([2 14],[y4 y4],'color',red,'LineWidth',1)
%     text(mean([2 14]),mean([y4 y5]),'??','fontsize',8,'color',red)
% elseif signif(2) == 1
%     plot([2 14],[y4 y4],'color',red,'LineWidth',1)
%     text(mean([2 14]),mean([y4 y5]),'??','fontsize',8,'color',red)
% end

end

function single_box_whisker_plot(only_data,plot_title,plot_ylabel,subplot_number)
% box_whisker_plot(left_data,right_data,plot_title,plot_ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 30;
 
box_data(:,1)= only_data{1}(end-SPAN+1:end);
box_data(1:SPAN,2)= nan; 
box_data(:,3)= only_data{2}(1:SPAN); 
box_data(1:SPAN,4)= nan; 
box_data(:,5)= only_data{2}(end-SPAN+1:end);
box_data(1:SPAN,6)= nan; 
box_data(:,7)= only_data{3}(1:SPAN);
box_data(1:SPAN,8)= nan; 
box_data(:,9)= only_data{3}(end-SPAN+1:end);

colors = purple(1:3);

subplot(3,4,subplot_number); hold on;
plot([-100 100],[0 0],'LineWidth',1,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off')
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([2 2],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([6 6],[-10000 10000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')

signif= significance_test(only_data);
ax= axis;
y1= ax(4)*1.05;
y2= ax(4)*1.25;
axis([ax(1) ax(2) ax(3) ax(4)*1.5])
if signif(1) == 1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(1) == -1
    plot([1 7],[y1 y1],'color',purple,'LineWidth',1)
    text(mean([1 7]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
end
if signif(2) == -1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 ax(4)*1.5]),'\boldmath$\star \downarrow$','fontsize',8,'color',purple,'interpreter','latex')
elseif signif(2) == 1
    plot([1 9],[y2 y2],'color',purple,'LineWidth',1)
    text(mean([1 9]),mean([y2 ax(4)*1.5]),'\boldmath$\star \uparrow$','fontsize',8,'color',purple,'interpreter','latex')
end

xticks([])
title(plot_title)
ylabel(plot_ylabel)
end

function box_whisker_plot(left_data,right_data,plot_title,plot_ylabel,subplot_number)
% box_whisker_plot(left_data,right_data,plot_title,plot_ylabel,subplot_number)
blue= [0,114/255,195/255,1];
red= [204/255,53/255,37/255,1];
bluet= [0,114/255,195/255,.1];
redt= [204/255,53/255,37/255,.1];
bluem= [0,114/255,195/255];
redm= [204/255,53/255,37/255];
purple= [163/255 41/255 214/255 1];
purplet= [163/255 41/255 214/255 .1];
purplem= [163/255 41/255 214/255];
SPAN= 30;
 
box_data(:,1)= left_data{1}(end-SPAN+1:end);
box_data(:,2)= right_data{1}(end-SPAN+1:end);
box_data(1:SPAN,3)= nan; 
box_data(:,4)= left_data{2}(1:SPAN); 
box_data(:,5)= right_data{2}(1:SPAN); 
box_data(1:SPAN,6)= nan; 
box_data(:,7)= left_data{2}(end-SPAN+1:end);
box_data(:,8)= right_data{2}(end-SPAN+1:end);
box_data(1:SPAN,9)= nan; 
box_data(:,10)= left_data{3}(1:SPAN);
box_data(:,11)= right_data{3}(1:SPAN);
box_data(1:SPAN,12)= nan; 
box_data(:,13)= left_data{3}(end-SPAN+1:end);
box_data(:,14)= right_data{3}(end-SPAN+1:end);

for i= 1:14
    if ismember(i,[1 4 7 10 13])
        colors(i,:) = red(1:3);
    elseif ismember(i,[2 5 8 11 14])
        colors(i,:) = blue(1:3);
    end
end

subplot(3,4,subplot_number); hold on;
boxplot(box_data,'outliersize',4)
set(findobj(gca,'type','line'),'linew',1)
set(findobj(gca,'type','line'),'color','k');
% set(findobj(gcf,'tag','Outliers'),'MarkerSize',25);
h = findobj(gca,'Tag','Box');
h1= findobj(gca,'Tag','Outliers');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
set(h1,'MarkerEdgeColor','k');
plot([3 3],[-20000 20000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
plot([9 9],[-20000 20000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
xticks([])
title(plot_title)
ylabel(plot_ylabel)


ax= axis;
axis([ax(1) ax(2) ax(3) ax(4)*1.5])
y1= ax(4)*1.05;
y2= ax(4)*1.15;
y3= ax(4)*1.25;
y4= ax(4)*1.35;
y5= ax(4)*1.45;

signif= significance_test(left_data);
if signif(1) == 1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(1) == -1
    plot([1 10],[y1 y1],'color',blue,'LineWidth',1)
    text(mean([1 10]),mean([y1 y2]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
end
if signif(2) == -1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \downarrow$','fontsize',8,'color',blue,'interpreter','latex')
elseif signif(2) == 1
    plot([1 13],[y2 y2],'color',blue,'LineWidth',1)
    text(mean([1 13]),mean([y2 y3]),'\boldmath$\star \uparrow$','fontsize',8,'color',blue,'interpreter','latex')
end

signif= significance_test(right_data);
if signif(1) == 1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(1) == -1
    plot([2 11],[y3 y3],'color',red,'LineWidth',1)
    text(mean([2 11]),mean([y3 y4]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
end
if signif(2) == -1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \downarrow$','fontsize',8,'color',red,'interpreter','latex')
elseif signif(2) == 1
    plot([2 14],[y4 y4],'color',red,'LineWidth',1)
    text(mean([2 14]),mean([y4 y5]),'\boldmath$\star \uparrow$','fontsize',8,'color',red,'interpreter','latex')
end

end

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

colororder('default');

subplot(3,4,subplot_number); hold on;
plot(par1x,par1y,'LineWidth',2)
plot(par4x,par4y,'LineWidth',2)
plot(par5x,par5y,'LineWidth',2)

axis equal
% axis manual
% plot([-1000 1000],[0 0],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% plot([0 0],[-1000 1000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% legend('Baseline','Early Adapt','Late Adapt','Early Obs','Late Obs','Position',[0.90827 0.23599 0.076389 0.11203])
legend('Baseline','Early Obs','Late Obs','Position',[0.90827 0.23599 0.076389 0.11203])
title(plot_title)
end

function double_cop_plot(only_data,plot_title,subplot_number)
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

colororder('default');

subplot(1,2,subplot_number); hold on;
plot(par1x,par1y,'LineWidth',2)
plot(par4x,par4y,'LineWidth',2)
plot(par5x,par5y,'LineWidth',2)
box on
% axis equal
% axis manual
% plot([-1000 1000],[0 0],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% plot([0 0],[-1000 1000],'LineWidth',2,'color',[.3 .3 .3],'HandleVisibility','off')
% legend('Baseline','Early Adapt','Late Adapt','Early Obs','Late Obs','Position',[0.90827 0.23599 0.076389 0.11203])
legend('Baseline','Early Obs','Late Obs','Position',[0.79879 0.14558 0.094804 0.12619])
title(plot_title)
xlabel('Medial/Lateral Position (mm)')
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

subplot(3,4,subplot_number); hold on;
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

subplot(3,4,subplot_number); hold on;
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


subplot(3,4,subplot_number); hold on;
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

function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')
%   Author : NS
%   Version: 3.0, 21 Sept. 2010
%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.
    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

function h_output = significance_test(only_data)
SPAN = 30;
a = only_data{1}(end-SPAN+1:end);
c1 = only_data{3}(1:SPAN);
c2 = only_data{3}(end-SPAN+1:end);

[~, h(1,1)] = ranksum(a,c1,'alpha',0.05,'tail','left');
[~, h(1,2)] = ranksum(a,c2,'alpha',0.05,'tail','left');
[~, h(2,1)] = ranksum(a,c1,'alpha',0.05,'tail','right');
[~, h(2,2)] = ranksum(a,c2,'alpha',0.05,'tail','right');

h_output= h(1,:)-h(2,:);
end

function h_output = significance_test_average(only_data)
a = only_data(:,1);
c1 = only_data(:,4);
c2 = only_data(:,5);

[~, h(1,1)] = ranksum(a,c1,'alpha',0.05,'tail','left');
[~, h(1,2)] = ranksum(a,c2,'alpha',0.05,'tail','left');
[~, h(2,1)] = ranksum(a,c1,'alpha',0.05,'tail','right');
[~, h(2,2)] = ranksum(a,c2,'alpha',0.05,'tail','right');

h_output= h(1,:)-h(2,:);
end