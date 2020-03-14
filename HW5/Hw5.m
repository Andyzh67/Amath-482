%% Part 1 Fully-connected neural network
load('fashion_mnist.mat')

X_train = im2double(X_train);
X_test = im2double(X_test);

X_train = reshape(X_train,[60000 28 28 1]);
X_train = permute(X_train,[2 3 4 1]);

X_test = reshape(X_test,[10000 28 28 1]);
X_test = permute(X_test,[2 3 4 1]);

X_valid = X_train(:,:,:,1:5000);
X_train = X_train(:,:,:,5001:end);

y_valid = categorical(y_train(1:5000))';
y_train = categorical(y_train(5001:end))';
y_test = categorical(y_test)';


layers = [imageInputLayer([28 28 1])
        fullyConnectedLayer(500)
        reluLayer
        fullyConnectedLayer(500)
        reluLayer
        fullyConnectedLayer(500)
        reluLayer
        fullyConnectedLayer(100)
        reluLayer
        fullyConnectedLayer(10)
        softmaxLayer
        classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',8,...
    'Shuffle','every-epoch',...
    'InitialLearnRate',1e-4, ...
    'L2Regularization',1e-3, ...
    'ValidationData',{X_valid,y_valid}, ...
    'Verbose',false, ...
    'Plots','training-progress')

net = trainNetwork(X_train,y_train,layers,options);
%%
y_pred = classify(net,X_train);
plotconfusion(y_train,y_pred)
%%
y_pred = classify(net,X_valid);
plotconfusion(y_valid,y_pred)
%%
y_pred = classify(net,X_test);
plotconfusion(y_test,y_pred)



%%
clear; close all; clc;
load('fashion_mnist.mat')


X_train = im2double(X_train);
X_test = im2double(X_test);

X_train = reshape(X_train,[60000 28 28 1]);
X_train = permute(X_train,[2 3 4 1]);

X_test = reshape(X_test,[10000 28 28 1]);
X_test = permute(X_test,[2 3 4 1]);

X_valid = X_train(:,:,:,1:5000);
X_train = X_train(:,:,:,5001:end);

y_valid = categorical(y_train(1:5000))';
y_train = categorical(y_train(5001:end))';
y_test = categorical(y_test)';

layers = [
    imageInputLayer([28 28 1],"Name","imageinput")
    convolution2dLayer([3 3],32,"Name","conv_1","Padding","same")
    batchNormalizationLayer
    tanhLayer
    maxPooling2dLayer([3 3],"Name","avgpool2d_1","Padding","same","Stride",[2 2])
    convolution2dLayer([3 3],64,"Name","conv_2")
    batchNormalizationLayer
    tanhLayer
    maxPooling2dLayer([3 3],"Name","avgpool2d_2","Padding","same","Stride",[2 2])
    convolution2dLayer([3 3],128,"Name","conv_3")  
    batchNormalizationLayer
    tanhLayer
    fullyConnectedLayer(100)
    reluLayer
    fullyConnectedLayer(100)
    reluLayer
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer("Name","classoutput")];

options = trainingOptions('adam', ...
    'MaxEpochs',8,...
    'Shuffle','every-epoch', ...
    'ValidationFrequency',30, ...
    'InitialLearnRate',1e-3, ...
    'L2Regularization',1e-3, ...
    'ValidationData',{X_valid,y_valid}, ...
    'Verbose',false, ...
    'Plots','training-progress')

net = trainNetwork(X_train,y_train,layers,options);
%%
y_pred = classify(net,X_train);
plotconfusion(y_train,y_pred)
%%
y_pred = classify(net,X_valid);
plotconfusion(y_valid,y_pred)
%%
y_pred = classify(net,X_test);
plotconfusion(y_test,y_pred)
