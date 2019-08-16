%%This function takes in "ID0X_raw.mat" file, and evaluate 1 day data
%% 
% 1.seperate one day data into 4 time zone 12pm-6pm, 6pm-12am,
% 12am-6am,6am-12pm
% 2.Based on the acc data, and threshold value, find high, med, low
% activities
% 3.plot the activities distribution for 4 time zones
% 4. save the time and magnitude of high,med,low activities for 4 time
% zones
clc;clear
fileID='DATA0000'; % walking data unavailable at the moment
%addpath('Walk_analysis_code')
%addpath('Sleep_analysis_code')
load([fileID,'_raw.mat'])
%% check model
Label = vertcat(zeros(6399,1), ones(6000,1), ones(6500,1)*2, ones(6000,1)*3, zeros(13100,1), ones(6000,1), ones(6500,1)*2, ones(6200,1)*3, zeros(11750,1),ones(6600,1), ones(6200,1)*2, ones(6600,1)*3, zeros(29780,1));
Data = horzcat(MyData,Label);
%% parameters
delay=1;
sampling_rate=100;
one_minute_points=60*sampling_rate;
one_day_points=24*60*60*sampling_rate;
one_hour_points=60*60*sampling_rate;
total_days=ceil(length(MyData)/one_day_points);
 %%  get all data from 3days
    CurrentdayData=MyData;
    Times=CurrentdayData(:,1);
    Acc_X=CurrentdayData(:,2);    Acc_Y=CurrentdayData(:,3);     Acc_Z=CurrentdayData(:,4);
    Gyr_X=CurrentdayData(:,5);    Gyr_Y=CurrentdayData(:,6);     Gyr_Z=CurrentdayData(:,7);
    R_Acc_XZ=sqrt( Acc_X.^2+ Acc_Z.^2);
    R_Gyr_XYZ=sqrt(Gyr_X.^2+Gyr_Y.^2+Gyr_Z.^2);
    R_Acc_XYZ=sqrt(Acc_X.^2+Acc_Y.^2+Acc_Z.^2);
    R_Gyr=R_Gyr_XYZ;

    gyr_median_vector=zeros(length(R_Gyr_XYZ),1);
    gyr_std_vector=zeros(length(R_Gyr_XYZ),1);
    acc_mean_vector=zeros(length(R_Acc_XZ),1);
    acc_var_vector=zeros(length(R_Acc_XZ),1);
    %% plot
    figure(1); clf;
    subplot(3,1,1)%split the plot into 3 figures. This 3 x1 matrix. The tirhd index tells you wwhich box (1st box)
    plot(Times,R_Gyr) %walking is measured exclusively by gyroscope data
    datetick('x','HH:MM PM','keepticks');
    xlim([Times(1) Times(end)])
    %ylim([0 max(R_Gyr)])

window_size=sampling_rate*2;
    for j=1:window_size:length(R_Gyr_XYZ)
        gyr_to_analyze=R_Gyr_XYZ(j:min(length(R_Gyr_XYZ),j+window_size-1));
        acc_to_analyze=R_Acc_XZ(j:min(length(R_Acc_XZ),j+window_size-1));
        
        acc_mean_value=mean(acc_to_analyze);
        acc_var_value=var(acc_to_analyze);
        acc_mean_vector(j:min(length(R_Acc_XZ),j+window_size-1),1)=acc_mean_value*ones(length(acc_to_analyze),1);
        acc_var_vector(j:min(length(R_Acc_XZ),j+window_size-1),1)=acc_var_value*ones(length(acc_to_analyze),1);
   
        gyr_median_value=median(gyr_to_analyze);
        gyr_std_value=std(gyr_to_analyze);
        gyr_median_vector(j:min(length(R_Gyr_XYZ),j+window_size-1),1)=gyr_median_value*ones(length(gyr_to_analyze),1);
        gyr_std_vector(j:min(length(R_Gyr_XYZ),j+window_size-1),1)=gyr_std_value*ones(length(gyr_to_analyze),1);
     
    end
   
 subplot(3,1,2) %
 
 high=zeros(length(R_Gyr_XYZ),1);
 med=zeros(length(R_Gyr_XYZ),1);
 low=zeros(length(R_Gyr_XYZ),1);
 
 %% threshold 
 median_threshold_high_med=55;
 median_threshold_med_low=20;
 %std_threshold_high_med=40;
 %std_threshold_med_low=15;
for k=1:window_size:length(R_Gyr_XYZ)
    if gyr_median_vector(k)>=median_threshold_high_med%||std_vector(k)>=std_threshold_high_med
        high(k:min(length(R_Gyr_XYZ),k+window_size-1),1)=ones(min(window_size,length(R_Gyr_XYZ)-k+1),1);
    elseif  gyr_median_vector(k)<median_threshold_med_low || (acc_mean_vector(k)>0.8 && acc_var_vector(k)<0.1)%|| std_vector(k)<std_threshold_med_low
       low(k:min(length(R_Gyr_XYZ),k+window_size-1),1)=ones(min(window_size,length(R_Gyr_XYZ)-k+1),1);
    else
         med(k:min(length(R_Gyr_XYZ),k+window_size-1),1)=ones(min(window_size,length(R_Gyr_XYZ)-k+1),1);
    end
    
end
plot(Times,low*1,'r.',Times,med*2,'g.',Times,high*3,'b.')
title(['Time zone',num2str(i),' activities']);
legend('low','med','high')
datetick('x','HH:MM PM');
xlim([Times(1) Times(end)])
ylim([0.1 4])


subplot(3,1,3)  %plot percentage
name = {'low';'med';'high'};
bar([sum(low)/length(R_Gyr_XYZ),sum(med)/length(R_Gyr_XYZ),sum(high)/length(R_Gyr_XYZ)])
set(gca,'xticklabel',name)

%% save data
walk_acc_vector=med.*R_Acc_XYZ;
walk.Racc=walk_acc_vector(find(walk_acc_vector~=0));
walk_gyr_vector=med.*R_Gyr_XYZ;
walk.Rgyr=walk_gyr_vector(find(walk_gyr_vector~=0));
walk.Data = Data(find(walk_acc_vector~=0),:);
save([fileID,'_walk.mat'],'walk')
%% check Precision & Recall
P = (sum(walk.Data(:,8)==1) + sum(walk.Data(:,8)==3))/length(walk.Data);
P
R = (sum(walk.Data(:,8)==1) + sum(walk.Data(:,8)==3))/(sum(Label == 1)+sum(Label == 3));
R

figure(2),
plot(Times,R_Gyr);
hold on
plot(Times,gyr_median_vector,'Linewidth',2);
hold on
plot(Times,med*30,'Linewidth',2);
hold on
plot(Times,Label*10,'Linewidth',2);