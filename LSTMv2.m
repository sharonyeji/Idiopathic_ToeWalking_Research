opengl('save', 'software')
load('DATA0000_raw.mat');

%Extract walking and toe walking data from 3 subjects
SharonW = MyData(6400:12399,:); %SharonW
SharonT = MyData(18900:24899,:); %SharonT
DavidW = MyData(38000:43999,:); %DavidW
DavidT = MyData(50500:56699,:); %DavidT
JoeyW = MyData(68450:75049,:); %JoeyW
JoeyT = MyData(81250:87849,:); %JoeyT

list = {SharonW, SharonT, DavidW, DavidT};
XTrain = {0,0,0,0};
XVal = {0,0,0,0};
XTest = horzcat(JoeyW(:,2:7).',JoeyT(:,2:7).');
YTrain = {0,0,0,0};
YVal = {0,0,0,0};
[m1 n1] = size(JoeyW);
[m2 n2] = size(JoeyT);
YTest = horzcat(categorical(repmat({'Walking'},m1,1)).',categorical(repmat({'Toe-Walking'},m2,1)).');

i = 1;
for j = list
    data = j{1};
    [m n] = size(data);
    [trainInd, valInd, testInd] = divideblock(m,0.8,0.2,0);
    XTrain{1,i} = data(trainInd,2:7).';
    XVal{1,i} = data(valInd,2:7).';
        if mod(i,2)==1 %if odd, walk. if even toe-walk
            YTrain{1,i} = categorical(repmat({'Walking'},m*0.8,1)).';
            YVal{1,i} = categorical(repmat({'Walking'},m*0.2,1)).';
        else
            YTrain{1,i} = categorical(repmat({'Toe-Walking'},m*0.8,1)).';
            YVal{1,i} = categorical(repmat({'Toe-Walking'},m*0.2,1)).';
        end
    i = i + 1;
end

inputSize = 6; %number of features
numHiddenUnits = 100; %Frequency, Hz
numClasses = 2; 

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

%To train on a GPU, if available, set 'ExecutionEnvironment' to 'auto' (this is the default value).
options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'ValidationData',{XVal,YVal}, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',20, ...
    'Verbose',0, ...
    'Plots','training-progress');

rng(6);
net = trainNetwork(XTrain,YTrain,layers,options);

YPred = classify(net,XTest);
acc = sum(YPred == YTest)./numel(YTest)

figure
plot(YPred,'.-')
hold on
plot(YTest)
hold off

xlabel("Time Step")
ylabel("Activity")
title("Predicted Activities")
legend(["Predicted" "Test Data"])