function offline_heel_strikes = fvespa(lhee_vert,lhee_sag,frames,velo_thld)
global real_time_heel_strikes;
%% Heel-strike Detection (FVESPA)
%Custom Algorithm FVESPA - Real-time Implementation
% The F-VESPA detects heel-strikes using the vertical and sagital position 
% of the posterior midsole of the heel(HEEL marker). Heel-strikes are defined
% as the frames that correspond to local minima of the vertical position 
% with a temporal prominence of three past and one future samples, and a 
% positive sagital velocity. For more information refer to the following
% paper: 10.1109/IROS51168.2021.9636335

% Creator: Chrysostomos Karakasis
% Last edit: 11/9/2023
% "Update F-VESPA parameters to accomodate the new Vicon origin on VST2.

% Fixed parameters-variables for the 2nd-order Butterworth Filter
Fs = 100; % (Hz)
fc = 20; % (Hz)
omega_c = 2*pi*fc;
T = 1/Fs;
b1 = (omega_c*T)^2;
b2 = 2*(omega_c*T)^2;
b3 = (omega_c*T)^2;

a1 = 4 + 2*sqrt(2)*(omega_c*T) + (omega_c*T)^2;
a2 = -8 + 2*(omega_c*T)^2;
a3 = 4 - 2*sqrt(2)*(omega_c*T) + (omega_c*T)^2;

% Initialize the initial values for the filtering so that it can be
% applied even for the first sample
lhee_vert_raw_past = [0 0];     % Values of the raw signal at the last two samples ((k-1) (k-2))
lhee_vert_filt_past = [0 0];    % Values of the filtered signal at the last two samples ((k-1) (k-2))
lhee_sag_raw_past = [0 0];      % Values of the raw signal at the last two samples ((k-1) (k-2))
lhee_sag_filt_past = [0 0];     % Values of the filtered signal at the last two samples ((k-1) (k-2))

y_filt_vert = [];               % Vector to store the filtered values of the vertical position of the left heel marker
y_filt_sag = [];                % Vector to store the filtered values of the sagittal position of the left heel marker

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization for the FVESPA
min_lhee=-1000; %Vertical component (in this case Y) of the LHEE marker
flag=0;
flag2=1;
mins=[];
maxs=[];
vel_prev_1=0;
vel_prev_2=0;
vel_prev_3=0;
computed_heel_strikes=[];
local_minima_fvespa = [];
sagittal_velocity = [];
lhee_sag_filtered = [];

iter_count = 1;
while(1) % Run loop forever until the experiment is over
    % Load new-fresh samples of the desired markers
    % The Left Heel marker will be utilized
    lhee_vert_new = lhee_vert(iter_count);  %Y axis of Left Heel Marker (motion on the vertical plane)
    lhee_sag_new = lhee_sag(iter_count);    %X axis of Left Heel Marker (motion in the sagittal plane)
    
    %We decide to use a 2nd order low pass butterworth filter at
    %cutoff 20Hz. The corresponding filter parameters were calculated in
    %the beginning. Since a 2nd order filter is applied, the last two
    %samples of both the raw and filtered signals are required.
    
    % Apply the filter to the vertical position of the left heel marker
    lhee_vert_new_f=(b1*lhee_vert_new + b2*lhee_vert_raw_past(1) + b3*lhee_vert_raw_past(2) - a2*lhee_vert_filt_past(1) - a3*lhee_vert_filt_past(2) )/a1;
    
    % Apply the filter to the sagittal position of the left heel marker
    lhee_sag_new_f=(b1*lhee_sag_new + b2*lhee_sag_raw_past(1) + b3*lhee_sag_raw_past(2) - a2*lhee_sag_filt_past(1) - a3*lhee_sag_filt_past(2) )/a1;
    
    lhee_sag_filtered = [lhee_sag_filtered; lhee_sag_new_f];
    
    %% Apply the F-VESPA algorithm
    % In the implementation of the proposed code on the data from STAR, the Y
    % axis was defined in a way, where as that foot was moving forward, the Y
    % axis was decreasing. Hence, just before heel strike, as the foot is moving
    % backwards, the derivative of Y was increasing. 
    
    vel_z = lhee_vert_new_f - lhee_vert_filt_past(1);                                       % Velocity in the vertical plane
    vel_s = lhee_sag_new_f - lhee_sag_filt_past(1);                                         % Velocity in the sagittal plane
    sagittal_velocity = [sagittal_velocity; vel_s];
    vel_z_thld = 0; % [0 0 0 0 0 0 0                                                                 % This is set to 0 in the real-time F-VESPA
    vel_s_thld = velo_thld; % [0 5 6.005 0 10 0]                                                 % This is set to 0 in the real-time F-VESPA
    if iter_count == 24710
        iter_count;
    end
    if vel_z>=vel_z_thld && vel_prev_1<=0 && vel_prev_2<=0 && vel_prev_3<=0 && flag2==0
        local_minima_fvespa = [local_minima_fvespa iter_count];
        if vel_s<=vel_s_thld && lhee_vert_new_f<500                                        % Necessary for Vicon H.S. and extra check that heel-strikes are not detected in swing phase
            % Here a new heel-strike was detected!!!
            min_lhee=lhee_vert_filt_past(1);
            %This was modified to stress that although we want to detect global minima, this detection has a delay of one frame
            computed_heel_strikes=[computed_heel_strikes,frames(iter_count)];
            flag2=1;
        end
    elseif iter_count>2 && vel_z<0 && vel_prev_1<=0 && vel_prev_2>=0 && ge(vel_prev_3,0) && ((lhee_vert_filt_past(2)-min_lhee)>100)
        maxs=[maxs; iter_count-2];
        flag2=0;
    end
    vel_prev_3=vel_prev_2;
    vel_prev_2=vel_prev_1;
    vel_prev_1=vel_z;
    %%% End of FVESPA algorithm for this iteration
    
    % Update the past samples vectors (raw and filtered) for the next incoming sample
    lhee_vert_raw_past = [lhee_vert_new lhee_vert_raw_past(1)];
    lhee_vert_filt_past = [lhee_vert_new_f lhee_vert_filt_past(1)];
    % Append the new filtered value into the total one for post processing
    y_filt_vert = [y_filt_vert lhee_vert_new_f];
    
    % Update the past samples vectors (raw and filtered) for the next incoming sample
    lhee_sag_raw_past = [lhee_sag_new lhee_sag_raw_past(1)];
    lhee_sag_filt_past = [lhee_sag_new_f lhee_sag_filt_past(1)];
    % Append the new filtered value into the total one for post processing
    y_filt_sag = [y_filt_sag lhee_sag_new_f];
    
    iter_count = iter_count + 1;
    if iter_count > length(lhee_vert)
        break
    end
end

vel_vert = [0 diff(y_filt_vert)];
vel_sag = [0 diff(y_filt_sag)];

%% Figures showing the detected heel-strikes for the vertical and sagital position and velocity

% offline_heel_strikes = real_time_heel_strikes';
% offline_heel_strikes(find(y_filt_sag(real_time_heel_strikes)>690)) = computed_heel_strikes(find(y_filt_sag(real_time_heel_strikes)>690));

% figure(2)
% ax(1) = subplot(4,1,1);
% hold on;
% plot(vel_vert)
% plot(computed_heel_strikes,vel_vert(computed_heel_strikes),'rx','Linewidth',2,'MarkerSize',10)
% plot(real_time_heel_strikes,vel_vert(real_time_heel_strikes),'kd','Linewidth',2,'MarkerSize',10)
% plot(1:length(vel_vert),zeros(size(vel_vert)),'k--')
% plot(1:length(vel_vert),vel_z_thld*ones(size(vel_vert)),'r--')
% legend('Vertical Velocity','Offline F-VESPA H.S.','Real-time F-VESPA H.S.','Threshold')
% title('Vertical Velocity of Left Heel Marker','interpreter','latex')
% ax(2) = subplot(4,1,2);
% hold on;
% plot(vel_sag)
% plot(computed_heel_strikes,vel_sag(computed_heel_strikes),'rx','Linewidth',2,'MarkerSize',10)
% plot(real_time_heel_strikes,vel_sag(real_time_heel_strikes),'kd','Linewidth',2,'MarkerSize',10)
% plot(1:length(vel_sag),0*ones(size(vel_sag)),'k--')
% plot(1:length(vel_sag),vel_s_thld*ones(size(vel_sag)),'r--')
% title('Sagital Velocity of Left Heel Marker','interpreter','latex')
% ax(3) = subplot(4,1,3);
% hold on;
% plot(y_filt_vert)
% plot(computed_heel_strikes,y_filt_vert(computed_heel_strikes),'rx')
% plot(real_time_heel_strikes,y_filt_vert(real_time_heel_strikes),'kd')
% plot(1:length(y_filt_vert),500*ones(size(y_filt_vert)),'r--')
% legend('Vertical Position','Offline F-VESPA H.S.','Real-time F-VESPA H.S.','Threshold')
% title('Vertical Position of Left Heel Marker','interpreter','latex')
% ax(4) = subplot(4,1,4);
% hold on;
% plot(y_filt_sag)
% plot(computed_heel_strikes,y_filt_sag(computed_heel_strikes),'rx')
% plot(real_time_heel_strikes,y_filt_sag(real_time_heel_strikes),'kd')
% title('Sagital Position of Left Heel Marker','interpreter','latex')
% linkaxes(ax,'x')

offline_heel_strikes = computed_heel_strikes;

% close all;
end