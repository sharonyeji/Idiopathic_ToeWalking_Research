clc;clear
fileID='DATA0000'; %Dynamic field, input file name
load([fileID,'_raw.mat'])
%% check model
Label = vertcat(zeros(5799,1), ones(6600,1), ones(6400,1)*2, ones(6400,1)*3, zeros(12400,1), ones(6400,1), ones(6300,1)*2, ones(6700,1)*3, zeros(11450,1),ones(6600,1), ones(6200,1)*2, ones(6600,1)*3, zeros(29780,1));
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
    Acc=CurrentdayData(:,2:4); Acc_X=CurrentdayData(:,2);    Acc_Y=CurrentdayData(:,3);     Acc_Z=CurrentdayData(:,4); 
    Gyr=CurrentdayData(:,5:7); Gyr_X=CurrentdayData(:,5);    Gyr_Y=CurrentdayData(:,6);     Gyr_Z=CurrentdayData(:,7);
    
    % calculate the most sensitive axis, M
    M = zeros(length(CurrentdayData)-sampling_rate-1,1);
    for ii = 1:length(CurrentdayData)-sampling_rate-1
        m = sum(abs(Gyr(ii:ii+sampling_rate-1,:)));
        if min(m) == m(1)
            M(ii) = Gyr_X(ii);
        elseif min(m) == m(2)
            M(ii) = Gyr_Y(ii);
        else
            M(ii) = Gyr_Z(ii);
        end
    end
    
    Times = Times(1:length(CurrentdayData)-sampling_rate-1,1);
    %% plot
    figure(3); clf;
    subplot(3,1,1)%split the plot into 3 figures. This 3 x1 matrix. The tirhd index tells you wwhich box (1st box)
    plot(Times,M) %walking is measured exclusively by gyroscope data
    datetick('x','HH:MM PM','keepticks');
    xlim([Times(1) Times(end)])
    ylim([0 max(M)])
    
    %% perform Fast Fourier Transform
    high=zeros(length(M),1);
    med=zeros(length(M),1);
    low=zeros(length(M),1);
    X = zeros(length(M),1);
    
    window_size = sampling_rate; %window to be capture 1. seconds
        for jj=1:window_size:length(M)
            X(jj:min(length(M),jj+window_size-1)) = fft(M(jj:min(length(M),jj+window_size-1)));
            XX(jj:jj+4) = abs(X(jj:jj+4));
            c = mean(X(2:3)); %mean of amplitude at frequency 0.6 to 2Hz
            o = mean(X(1:2)); %mean of amplitude at frequency 0 to 0.6Hz
            if c < o
                low(jj:min(length(M),jj+window_size-1),1) = ones(min(window_size,length(M)-jj+1),1);
            elseif c > o && c < 10
                med(jj:min(length(M),jj+window_size-1),1) = ones(min(window_size,length(M)-jj+1),1);
            elseif c > o && c >= 10
                high(jj:min(length(M),jj+window_size-1),1) = ones(min(window_size,length(M)-jj+1),1);
            end
        end    

XXX=XX(find(XX~=0),:);
n = length(X);
f = (0:n-1)*100/n;
figure(4)
plot(f,abs(X(1:n)))
%hold on
%plot(f, abs(A(1:n)))
plot(XX(12300:12500))

figure(3),subplot(3,1,2)        
plot(Times,low*1,'r.',Times,med*2,'g.',Times,high*3,'b.')
title(['Time zone',num2str(1i),' activities']);
legend('low','med','high')
datetick('x','HH:MM PM');
xlim([Times(1) Times(end)])
ylim([0.1 4])

subplot(3,1,3)  %plot percentage
name = {'low';'med';'high'};
bar([sum(low)/length(M),sum(med)/length(M),sum(high)/length(M)])
set(gca,'xticklabel',name)
%% save data
walk_vector=high.*M;
walk.Data = Data(find(walk_vector~=0),:);
save([fileID,'_walk.mat'],'walk')
%% check Precision & Recall
P = (sum(walk.Data(:,8)==1) + sum(walk.Data(:,8)==3))/length(walk.Data);
P
R = (sum(walk.Data(:,8)==1) + sum(walk.Data(:,8)==3))/(sum(Label == 1)+sum(Label == 3));
R