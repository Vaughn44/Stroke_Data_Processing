clear all
close all
clc

% Import Data
for j= 1:6
    file_name= ['Data/Subject ' num2str(j) '/fm_data.mat'];
    load(file_name);
    fm{j}= fm_data;
    clear fm_data;
end    

% Calculations
for j= 1:6
    for i= 1:length(fm{j})
    max_value{j}(i,1)= max(fm{j}{i},[],'all');
    num_saturated{j}(i,1)= sum(fm{j}{i}>254,'all');
    num_active{j}(i,1)= sum(fm{j}{i}>0,'all');
    percent_saturated{j}(i,1)= sum(fm{j}{i}>254,'all')/sum(fm{j}{i}>0,'all')*100;
    end
end

% Plots
figure; set(gcf,'color','w'); hold on;
for j= 1:6
    plot(num_saturated{j},'LineWidth',2)
end
legend('Subject 1','Subject 2','Subject 3','Subject 4','Subject 5','Subject 6')
xlabel('Time Step')
ylabel('Number of Cells Saturated')
axis([3140 3290 0 30])

figure; set(gcf,'color','w'); hold on;
for j= 1:6
    plot(percent_saturated{j},'LineWidth',2)
end
legend('Subject 1','Subject 2','Subject 3','Subject 4','Subject 5','Subject 6')
xlabel('Time Step')
ylabel('Percent of Cells Saturated')
axis([3140 3290 0 30])