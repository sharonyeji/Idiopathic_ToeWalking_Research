%% K Means Clustering
% K Means Clustering is a simple unsupervised learning method.
% K Means Clustering does not have training/testing data.  It uses
% centroids to calculate how to group data.  Labels are not used in this
% clustering method.  The purpose of K Means Clustering is to find groups
% and patterns in the data that have not yet been explicitly labeled. In
% our case, we will use K Means to see if there are distinct differences in
% toe walking/walking metrics.

load('DATA0000_raw.mat')

% subset 2 minute of regular walking & toe walking data
WalkData = horzcat(MyData(66875:78874,:));
ToeWalkData = horzcat(MyData(96975:108974,:));

%% Feature Extraction
% An accelerometer and gyroscope with a frequency of 100Hz were used to 
% collect data.  The magnitude of the accelerometer and the magnitude of the 
% gyroscope were calculated and normalized.  One gait cycle is on average 1
% second, so the mean, standard deviation, covariance, root mean square, 
% kurtosis, and skew were caculated for every second of normalized data.

%Feature Extraction
%Combine WalkData and ToeWalkData
D = vertcat(WalkData(:,2:7), ToeWalkData(:,2:7));

%Calculate magnitude of Accelerometer and Gyroscope
%Find magnitude of Accelerometer/Gyroscope
RA = zeros(24000,1);
RG = zeros(24000,1);
for i=1:size(D,1)
    RA(i,1) = norm(D(i,1:3));
    RG(i,1) = norm(D(i,4:6));
end
normRA = (RA - min(RA))/(max(RA)-min(RA));
normRG = (RG - min(RG))/(max(RG)-min(RG));
normD = horzcat(normRA, normRG);
%% Magnitude of Accelerometer and Gyroscope
figure(1),
subplot(2,1,1), plot(WalkData(:,1),RA(1:12000,:),'g');
title('Magnitude of Accelerometer Walking Data')
subplot(2,1,2), plot(ToeWalkData(:,1),RA(12001:24000,:));
title('Magnitude of Accelerometer Toe Walking Data')

figure(2),
subplot(2,1,1), plot(WalkData(:,1),RG(1:12000,:),'g');
title('Magnitude of Gyroscope Walking Data')
subplot(2,1,2), plot(ToeWalkData(:,1),RG(12001:24000,:));
title('Magnitude of Gyroscope Toe Walking Data')
%%
%Create matrix that contains the following features:
%Mean
%Standard Deviation
%Covariance
%Root Mean Square
%Kurtosis
%Skewness
ME = zeros(240,1);
ST = zeros(240,1);
C = zeros(240,1);
R = zeros(240,1);
K = zeros(240,1);
S = zeros(240,1);

for a = 1:(size(normD,1)/100)
    ME(a,1) = mean2(normD((((a-1)*100)+1):a*100,:));
    ST(a,1) = std2(normD((((a-1)*100)+1):a*100,:));
    C(a,1) = cov(vertcat(normD(((a-1)*100)+1:a*100,1),normD(((a-1)*100)+1:a*100,2)));
    R(a,1) = rms(vertcat(normD(((a-1)*100)+1:a*100,1),normD(((a-1)*100)+1:a*100,2)));
    K(a,1) = kurtosis(vertcat(normD(((a-1)*100)+1:a*100,1),normD(((a-1)*100)+1:a*100,2)));
    S(a,1) = skewness(vertcat(normD(((a-1)*100)+1:a*100,1),normD(((a-1)*100)+1:a*100,2)));
end

normME = (ME - min(ME))/(max(ME)-min(ME));
normST = (ST - min(ST))/(max(ST)-min(ST));
normC = (C - min(C))/(max(C)-min(C));
normR = (R - min(R))/(max(R)-min(R));
normK = (K - min(K))/(max(K)-min(K));
normS = (S - min(S))/(max(S)-min(S));

F = horzcat(normME, normST, normC, normR, normK, normS);
%% Feature Selection
% After calculating the features, the VIF factors for each feature were
% measured to check for multicollinearity.  The features that had high VIFs
% were removed, and the dataset was checked for multicollinearity again.
% This step was repeated until 3 features were left. They are the standard
% deviation, the RMS, and the Kurtosis.

%Check for multicollinearity
x = vif_function(F);
%VIF for mean is extremely high, so remove feature and run again.

F = horzcat(normST, normC, normR, normK, normS);
x = vif_function(F);
%VIF for covariance is high, so remove feature and run again.

F = horzcat(normST, normR, normK, normS);
x = vif_function(F);

F = horzcat(normST, normR, normK);
x = vif_function(F);
%All VIFs are < 10, continue with K Means Clustering
%Final Features are:
%1. Standard Deviation
%2. RMS
%3. Kurtosis

%Include Activity Label
    %Label 1 = regular walking
    %Label 2 = toe walking

Label = vertcat(ones(120,1), ones(120,1)*2);
%% Perform K Means Clustering on Feature (F) data
% Finally, k-means clustering was performed on the dataset.  Because there 
% are two groups, toe walkers and regular walkers, we know to set the cluster
% to 2. 'Cityblock distance' was used instead of 'Euclidean distance' because 
% the cityblock method allowed for the data to be separated by clusters 
% vertically instead of horizontally.  The k-means clustering method was 
% replicated 10 times to make sure the smallest total sum of distances was 
% chosen for the centroids.
rng(1);
opts = statset('Display','final');
[cidx2, cmeans2] = kmeans(F, 2,'Distance','cityblock','Replicates',10,'Options',opts);
figure(3), 
[silh2, h] = silhouette(F, cidx2, 'cityblock');
mean_Sil = mean(silh2);
%Average Silhouette score for 'cityblock' 0.31312
%% Plot results
% The results show the clusters that were formed for toe-walkers and normal
% walkers. On the left is the k-means cluster that was computed, and on the
% right is the true cluster. Because we have the true values of the
% classification, an accuracy score was computed to be around 82%.
figure(4), subplot(1,2,1),
ptsymb = {'b','r'};
for i = 1:2
    clust = find(cidx2==i);
    scatter3(F(clust,1),F(clust,2),F(clust,3),ptsymb{i});
    hold on
end
scatter3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko');
scatter3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'kx');
hold off
view(-20,20)
title('K means clustering');
xlabel('St. Dev');
ylabel('RMS');
zlabel('Kurtosis');

subplot(1,2,2),
ptsymb = {'b','r'};
for i = 1:2
    clust2 = find(Label==i);
    scatter3(F(clust2,1),F(clust2,2),F(clust2,3),ptsymb{i});
    hold on
end
scatter3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'ko');
scatter3(cmeans2(:,1),cmeans2(:,2),cmeans2(:,3),'kx');
hold off
view(-20,20)
title('True classification');
xlabel('St. Dev');
ylabel('RMS');
zlabel('Kurtosis');

accuracy_score = (sum(cidx2(1:120)==1)+sum(cidx2(121:240)==2))/240
sensitivity_score = 120/(120 + (sum(cidx2(121:240)==1)))
specificity_score = 120/(120 + (sum(cidx2(1:120)==2)))

%accuracy score of 81.67%
%% Conclusion
% The k means method is not a robust unsupervised learning method, mainly
% because determining the centroids/clusters highly depends on the data
% that is given.  However, based on the results that we obtained shows that
% other classification methods will provide great results.

%publish('Exploration.m','pdf')