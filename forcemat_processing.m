function [fm_data_processed,fm_parameters]= forcemat_processing(fm_data,lank_x,lank_y,lhee_x,lhee_y,ltoe_x,ltoe_y,rank_x,rank_y,rhee_x,rhee_y,rtoe_x,rtoe_y)
% [fm_data_processed,fm_parameters]= 
%  forcemat_processing(fm_data,lank_x,lank_y,lhee_x,lhee_y,ltoe_x,ltoe_y,rank_x,rank_y,rhee_x,rhee_y,rtoe_x,rtoe_y)
%
% VICON data should already be truncated and resampled at the same sampling
% rate (up-sampled or down-sampled)

% Check sampling of VICON & forcemat data
if length(lank_x) ~= length(fm_data)
    error('VICON data and forcemat data are not the same length. Resample them.')
else
    frame_total= length(lank_x);
end

% Reshape Forcemats into grid
for i= 1:frame_total
    fm_grid{i,1}= reshape(fm_data(i,:),96,88);
    fm_grid_left{i,1}= fm_grid{i}(:,1:44);
    fm_grid_right{i,1}= fm_grid{i}(:,45:88);
end

% Find "force" on left & right belts (not the final force value)
for i= 1:frame_total
    force_left(i,1)= sum(fm_grid_left{i},'all');
    force_right(i,1)= sum(fm_grid_right{i},'all');
end

% Left heel strike & toe off using forcemats 
clear ind peaks valleys
signal= smooth(force_left,10); % filter forcemat data
[valleys(:,1), ind(:,1)]= findpeaks(-signal,'MinPeakHeight',-350,'MinPeakWidth',1,'MinPeakDistance',100);
valleys= -valleys; % find valley values for forcemats (magnitude during swing)
[peaks(:,1),~]= findpeaks(signal,'MinPeakHeight',3000,'MinPeakWidth',1,'MinPeakDistance',50);
threshold= (mean(peaks)/1.1-mean(valleys))*.05 + mean(valleys); % finding rough 5% BW threshold
for i= 1:length(ind)
    temp_hs= ind(i); % temporary heel strike variable
    temp_to= ind(i); % temporary toe off variable
    while signal(temp_hs)<threshold && temp_hs < frame_total % edge case for end of experiment
        temp_hs= temp_hs+1; % if the value isn't above the threshold move forward 1 frame
    end
    while signal(temp_to)<threshold && temp_to > 0 % end case for beginning of experiment
        temp_to= temp_to-1; % if the value isn't above the threshold move backward 1 frame
    end
    lhs(i,1)= temp_hs;
    lto(i,1)= temp_to;
end

% Right heel strike & toe off using forcemats 
clear ind peaks valleys
signal= smooth(force_right,10); % filter forcemat data
[valleys(:,1), ind(:,1)]= findpeaks(-signal,'MinPeakHeight',-350,'MinPeakWidth',1,'MinPeakDistance',100);
valleys= -valleys; % find valley values for forcemats (magnitude during swing)
[peaks(:,1),~]= findpeaks(signal,'MinPeakHeight',3000,'MinPeakWidth',1,'MinPeakDistance',50);
threshold= (mean(peaks)/1.1-mean(valleys))*.05 + mean(valleys); % finding rough 5% BW threshold
for i= 1:length(ind)
    temp_hs= ind(i); % temporary heel strike variable
    temp_to= ind(i); % temporary toe off variable
    while signal(temp_hs)<threshold && temp_hs < frame_total % edge case for end of experiment
        temp_hs= temp_hs+1; % if the value isn't above the threshold move forward 1 frame
    end
    while signal(temp_to)<threshold && temp_to > 0 % end case for beginning of experiment
        temp_to= temp_to-1; % if the value isn't above the threshold move backward 1 frame
    end
    rhs(i,1)= temp_hs;
    rto(i,1)= temp_to;
end

% Correct for edge case of beginning or end of experiment
if lto(1)== 1
    lto(1)= [];
end
if rto(1)== 1
    rto(1)= [];
end
if lhs(end) == frame_total
    lhs(end)= [];
end
if rhs(end) == frame_total
    rhs(end)= [];
end
%% Determining Contact (Stance & Swing)

% Calculate left contact
if lhs(1)>lto(1) % starts in stance
    hs_counter = 1;
    lcontact= ones(frame_total,1); % preallocate
else % starts in swing
    hs_counter = 0;
    lcontact= zeros(frame_total,1); % preallocate
end
to_counter = 1;
for i= 2:frame_total
    if any(lhs==i)
        hs_counter = hs_counter + 1;
    end
    if any(lto==i)
        to_counter = to_counter + 1;
    end
    
    if hs_counter == to_counter
        state = 1; % stance detected
    elseif to_counter == hs_counter + 1
        state = 0; % swing detected
    else
        error('Error in left contact calculations. Heel strike & toe off out of sync.')
        break
    end
    lcontact(i,1) = state;
end

% Calculate right contact
if rhs(1)>rto(1) % starts in stance
    hs_counter = 1;
    rcontact= ones(frame_total,1); % preallocate
else % starts in swing
    hs_counter = 0;
    rcontact= zeros(frame_total,1); % preallocate
end
to_counter = 1;
for i= 2:frame_total
    if any(rhs==i)
        hs_counter = hs_counter + 1;
    end
    if any(rto==i)
        to_counter = to_counter + 1;
    end
    
    if hs_counter == to_counter
        state = 1; % stance detected
    elseif to_counter == hs_counter + 1
        state = 0; % swing detected
    else
        error('Error in right contact calculations. Heel strike & toe off out of sync.')
        break
    end
    rcontact(i,1) = state;
end

%% Process Forcemats & Add to data
% Translation constants
c= 10.21; % distance between cells in mm
x_1_offset= -7.4; % x offset of left force mat origin in Vicon frame
y_1_offset= -54.43; % y offset of left force mat origin in Vicon frame
x_2_offset= 467.5; % x offset of right force mat origin in Vicon frame
y_2_offset= -57.55; % y offset of right force mat origin in Vicon frame

% filter left forcemat data
for i= 1:frame_total
    if lcontact(i)==0 % if in swing set all values to 0
        fm_grid_left{i} = zeros(96,44);
    else % if in stance set all values far from the foot to 0
        x_limits= [min([lank_x(i) ltoe_x(i)])-50 max([lhee_x(i) ltoe_x(i)])+75];
        y_limits= [lhee_y(i)-10 ltoe_y(i)+100];
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
    if rcontact(i)==0 % if in swing set all values to 0
        fm_grid_right{i} = zeros(96,44);
    else % if in stance set all values far from the foot to 0
        x_limits= [min([rank_x(i) rtoe_x(i)])-50 max([rhee_x(i) rtoe_x(i)])+75];
        y_limits= [rhee_y(i)-10 rtoe_y(i)+100];
        x_limits= round((x_limits - x_2_offset) / c + 1);
        y_limits= round((y_limits - y_2_offset) / c + 1);
        fm_grid_right{i}(:,1:x_limits(1)) = 0; % all cells too far right
        fm_grid_right{i}(:,x_limits(2):end) = 0; % all cells too far right
        fm_grid_right{i}(1:y_limits(1),:) = 0; % all cells too far backward
        fm_grid_right{i}(y_limits(2):end,:) = 0; % all cells too far forward
    end
end
%% Scale Forcemat Data

% Recalculate "force" on left and right (accurate just not scaled)
for i= 1:length(fm_grid)
    force_left(i,1)= sum(fm_grid_left{i},'all');
    force_right(i,1)= sum(fm_grid_right{i},'all');
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

%% Set Outputs
fm_data_processed= fm_grid;
fm_parameters.lcontact= lcontact;
fm_parameters.rcontact= rcontact;
fm_parameters.lhs= lhs;
fm_parameters.rhs= rhs;
fm_parameters.lto= lto;
fm_parameters.rto= rto;
fm_parameters.cop_x= COP_x;
fm_parameters.cop_y= COP_y;
fm_parameters.total_force= fm_force;
fm_parameters.cop_x_left= COP_x_left;
fm_parameters.cop_y_left= COP_y_left;
fm_parameters.cop_x_rihgt= COP_x_right;
fm_parameters.cop_y_right= COP_y_right;
fm_parameters.force_left= fm_force_left;
fm_parameters.force_right= fm_force_right;
end