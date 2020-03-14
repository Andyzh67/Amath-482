clear; close all; clc;

% Case 1 - ideal

load('cam1_1.mat')
load('cam2_1.mat') 
load('cam3_1.mat')


[X1(1), Y1(1), ~, frameNum1(1)] = size(vidFrames1_1);
[X1(2), Y1(2), ~, frameNum1(2)] = size(vidFrames2_1);
[X1(3), Y1(3), ~, frameNum1(3)] = size(vidFrames3_1);

%{
numFrames = size(vidFrames3_1,4);
for j = 1:numFrames
    X = vidFrames3_1(:,:,:,j);
    imshow(X); drawnow
end
%}
%%
% trace the flashlight for each frame for each camera angle
for j = 1 : frameNum1(1)
    M1 = double(rgb2gray(vidFrames1_1(:, :, :, j)));
    M1(:, [1:300, 400:end]) = 0;
    maxVal = max(M1(:));
    [yval, xval] = find(M1 >= 0.95*maxVal);
    x1(j, 1) = mean(xval);
    y1(j, 1) = mean(yval); 
    %{
    pcolor(M1)
    shading interp
    drawnow
    colormap(gray)
    %}
end

for j = 1 : frameNum1(2)
    M2 = double(rgb2gray(vidFrames2_1(:, :, :, j)));
    M2(:, [1:200, 400:end]) = 0;
    maxVal = max(M2(:));
    [yval, xval] = find(M2 >= 0.95*maxVal);
    x2(j, 1) = mean(xval);
    y2(j, 1) = mean(yval);
    %{
    pcolor(M2)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% rotate back by swapping xval and yval
for j = 1 : frameNum1(3)
    M3 = double(rgb2gray(vidFrames3_1(:, :, :, j)));
    M3(:, [1:240, 480:end]) = 0;
    M3([1:200, 350:end], :) = 0;
    maxVal = max(M3(:));
    [yval, xval] = find(M3 >= 0.95*maxVal);
    x3(j, 1) = mean(yval);
    y3(j, 1) = mean(xval); 
    %{
    pcolor(M3)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% cut each data so that the can
% starts from the highest point
x1 = x1(30 : end);
y1 = y1(30 : end);
x2 = x2(39 : end);
y2 = y2(39 : end);
x3 = x3(30 : end);
y3 = y3(30 : end);

% organize data with same length
minlength = min([length(y1), length(y2), length(y3)]);
x1 = x1(1 : minlength);
y1 = y1(1 : minlength);
x2 = x2(1 : minlength);
y2 = y2(1 : minlength);
x3 = x3(1 : minlength);
y3 = y3(1 : minlength);

Mx = [x1, y1, x2, y2, x3, y3];
Mx = Mx';
frame = 1:minlength;

% examine for the path of paint can and appropriate
% threshold values
figure(1)
for i = 1 : 3
    subplot(3, 1, i)
    plot(frame, Mx(2 * i - 1, :), frame, Mx(2 * i, :))
    legend('horizontal position', 'vertical position', 'FontSize', 12, 'Location', 'bestoutside')
    xlabel('frame', 'FontSize', 12)
    ylabel('position in camera', 'FontSize', 12)
    title(['path in ideal case (camera',num2str(i),')'], 'FontSize', 12)
end

% subtract the mean
[m, n] = size(Mx);
mu = mean(Mx, 2);
Mx = Mx - repmat(mu, 1, n);


% perform PCA
[U, S, V] = svd(Mx / sqrt(n-1), 'econ');
Y = U' * Mx;

figure(2)
plot(frame, Y)
xlabel('frame', 'FontSize', 12)
ylabel('normalized position of paint can', 'FontSize', 12)
title('principal component in ideal case', 'FontSize', 12)
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5', 'mode 6', 'FontSize', 12, 'Location', 'bestoutside')
grid on

figure(3)
sig = diag(S);
plot(sig.^2/sum(sig.^2),'ro','Linewidth',2)
xticks([0:1:6])
xlabel('energy modes', 'FontSize', 12)
ylabel('energy percentage (%)', 'FontSize', 12)
title('energy of modes in ideal case', 'FontSize', 12)


%% Case 2 - noisy
clear; close all; clc;

load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

[X2(1), Y2(1), ~, frameNum2(1)] = size(vidFrames1_2);
[X2(2), Y2(2), ~, frameNum2(2)] = size(vidFrames2_2);
[X2(3), Y2(3), ~, frameNum2(3)] = size(vidFrames3_2);
%{
numFrames = size(vidFrames3_2,4);
for j = 1:numFrames
    X = vidFrames3_2(:,:,:,j);
    imshow(X); drawnow
end
%}
%%
% trace the flashlight for each frame for each camera angle
for j = 1 : frameNum2(1)
    M1 = double(rgb2gray(vidFrames1_2(:, :, :, j)));
    M1(:, [1:300, 400:end]) = 0;
    M1(1:200, :) = 0;
    maxVal = max(M1(:));
    [yval, xval] = find(M1 >= 0.95*maxVal);
    x1(j, 1) = mean(xval);
    y1(j, 1) = mean(yval);
    %{
    pcolor(M1)
    shading interp
    drawnow
    colormap(gray)
    %}
end

for j = 1 : frameNum2(2)
    M2 = double(rgb2gray(vidFrames2_2(:, :, :, j)));
    M2(:, [1:160, 480:end]) = 0;
    M2(400:end, :) = 0;
    maxVal = max(M2(:));
    [yval, xval] = find(M2 >= 0.96*maxVal);
    x2(j, 1) = mean(xval);
    y2(j, 1) = mean(yval); 
    %{
    pcolor(M2)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% rotate back by swopping xval and yval
for j = 1 : frameNum2(3)
    M3 = double(rgb2gray(vidFrames3_2(:, :, :, j)));
    M3(:, [1:260, 500:end]) = 0;
    M3([1:160, 360:end], :) = 0;
    maxVal = max(M3(:));
    [yval, xval] = find(M3 >= 0.93*maxVal);
    x3(j, 1) = mean(yval);
    y3(j, 1) = mean(xval);  
    %{
    pcolor(M3)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% cut each data so that the can
% starts from the highest point
x1 = x1(12 : end);
y1 = y1(12 : end);
x2 = x2(5 : end);
y2 = y2(5 : end);
x3 = x3(15 : end);
y3 = y3(15 : end);

% organize data with same length
minlength = min([length(y1), length(y2), length(y3)]);
x1 = x1(1 : minlength);
y1 = y1(1 : minlength);
x2 = x2(1 : minlength);
y2 = y2(1 : minlength);
x3 = x3(1 : minlength);
y3 = y3(1 : minlength);

Mx = [x1, y1, x2, y2, x3, y3];
Mx = Mx';
frame = 1:minlength;

% examine for the path of paint can and appropriate
% threshold values
figure(1)
for i = 1 : 3
    subplot(3, 1, i)
    plot(frame, Mx(2 * i - 1, :), frame, Mx(2 * i, :))
    legend('horizontal position', 'vertical position', 'FontSize', 12, 'Location', 'bestoutside')
    xlabel('frame', 'FontSize', 12)
    ylabel('position in camera', 'FontSize', 12)
    title(['path in noisy case (camera',num2str(i),')'], 'FontSize', 12)
end

% subtract the mean
[m, n] = size(Mx);
mu = mean(Mx, 2);
Mx = Mx - repmat(mu, 1, n);


% perform PCA
[U, S, V] = svd(Mx / sqrt(n-1), 'econ');
Y = U' * Mx;

figure(2)
plot(frame, Y)
xlabel('frame', 'FontSize', 12)
ylabel('normalized position of paint can', 'FontSize', 12)
title('principal component in noisy case', 'FontSize', 12)
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5', 'mode 6', 'FontSize', 12, 'Location', 'bestoutside')
grid on

figure(3)
sig = diag(S);
plot(sig.^2/sum(sig.^2),'ro','Linewidth',2)
xticks([0:1:6])
xlabel('energy modes', 'FontSize', 12)
ylabel('energy percentage (%)', 'FontSize', 12)
title('energy of modes in noisy case', 'FontSize', 12)


%%
clear; close all; clc;

% Case 3 - horizontal displacement
load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')


[X3(1), Y3(1), ~, frameNum3(1)] = size(vidFrames1_3);
[X3(2), Y3(2), ~, frameNum3(2)] = size(vidFrames2_3);
[X3(3), Y3(3), ~, frameNum3(3)] = size(vidFrames3_3);

%%
numFrames = size(vidFrames1_3,4);
for j = 1:numFrames
    X = vidFrames1_3(:,:,:,j);
    imshow(X); drawnow
end
%%
% trace the flashlight for each frame for each camera angle
for j = 1 : frameNum3(1)
    M1 = double(rgb2gray(vidFrames1_3(:, :, :, j)));
    M1(:, [1:270, 400:end]) = 0;
    M1(1:230, :) = 0;
    maxVal = max(M1(:));
    [yval, xval] = find(M1 >= 0.95*maxVal);
    x1(j, 1) = mean(xval);
    y1(j, 1) = mean(yval);
    %{
    pcolor(M1)
    shading interp
    drawnow
    colormap(gray)
    %}
end

for j = 1 : frameNum3(2)
    M2 = double(rgb2gray(vidFrames2_3(:, :, :, j)));
    M2(:, [1:200, 440:end]) = 0;
    M2(1:120, :) = 0;
    maxVal = max(M2(:));
    [yval, xval] = find(M2 >= 0.98*maxVal);
    x2(j, 1) = mean(xval);
    y2(j, 1) = mean(yval); 
    %{
    pcolor(M2)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% rotate back by swopping xval and yval
for j = 1 : frameNum3(3)
    M3 = double(rgb2gray(vidFrames3_3(:, :, :, j)));
    M3(:, [1:200, 480:end]) = 0;
    M3([1:120, 360:end], :) = 0;
    maxVal = max(M3(:));
    [yval, xval] = find(M3 >= 0.95*maxVal);
    x3(j, 1) = mean(yval);
    y3(j, 1) = mean(xval);  
    %{
    pcolor(M3)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% cut each data so that the can
% starts from the highest point
x1 = x1(9 : end);
y1 = y1(9 : end);
x2 = x2(9 : end);
y2 = y2(9 : end);
x3 = x3(11 : end);
y3 = y3(11 : end);

% organize data with same length
minlength = min([length(y1), length(y2), length(y3)]);
x1 = x1(1 : minlength);
y1 = y1(1 : minlength);
x2 = x2(1 : minlength);
y2 = y2(1 : minlength);
x3 = x3(1 : minlength);
y3 = y3(1 : minlength);

frame = 1:minlength;
Mx = [x1, y1, x2, y2, x3, y3];
Mx = Mx';

% examine for the path of paint can and appropriate
% threshold values
figure(1)
for i = 1 : 3
    subplot(3, 1, i)
    plot(frame, Mx(2 * i - 1, :), frame, Mx(2 * i, :))
    legend('horizontal position', 'vertical position', 'FontSize', 12, 'Location', 'bestoutside')
    xlabel('frame', 'FontSize', 12)
    ylabel('position in camera', 'FontSize', 12)
    title(['path in horizontal displacement case (camera',num2str(i),')'], 'FontSize', 12)
end

% subtract the mean
[m, n] = size(Mx);
mu = mean(Mx, 2);
Mx = Mx - repmat(mu, 1, n);


% perform PCA
[U, S, V] = svd(Mx / sqrt(n-1), 'econ');
Y = U' * Mx;

figure(2)
plot(frame, Y)
xlabel('frame', 'FontSize', 12)
ylabel('normalized position of paint can', 'FontSize', 12)
title('principal component in horizontal displacement case', 'FontSize', 12)
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5', 'mode 6', 'FontSize', 12, 'Location', 'bestoutside')
grid on

figure(3)
sig = diag(S);
plot(sig.^2/sum(sig.^2),'ro','Linewidth',2)
xticks([0:1:6])
xlabel('energy modes', 'FontSize', 12)
ylabel('energy percentage (%)', 'FontSize', 12)
title('energy of modes in horizontal displacement case', 'FontSize', 12)

%%
clear; close all; clc;

% Case 4 - horizontal displacement and rotation
load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')


[X4(1), Y4(1), ~, frameNum4(1)] = size(vidFrames1_4);
[X4(2), Y4(2), ~, frameNum4(2)] = size(vidFrames2_4);
[X4(3), Y4(3), ~, frameNum4(3)] = size(vidFrames3_4);

%%
numFrames = size(vidFrames2_4,4);
for j = 1:numFrames
    X = vidFrames2_4(:,:,:,j);
    imshow(X); drawnow
end
%%
% trace the flashlight for each frame for each camera angle
for j = 1 : frameNum4(1)
    M1 = double(rgb2gray(vidFrames1_4(:, :, :, j)));
    M1(:, [1:300, 480:end]) = 0;
    M1([1:200, 420:end], :) = 0;
    maxVal = max(M1(:));
    [yval, xval] = find(M1 >= 0.95*maxVal);
    x1(j, 1) = mean(xval);
    y1(j, 1) = mean(yval); 
    %{
    pcolor(M1)
    shading interp
    drawnow
    colormap(gray)
    %}
end

for j = 1 : frameNum4(2)
    M2 = double(rgb2gray(vidFrames2_4(:, :, :, j)));
    M2(:, [1:210, 420:end]) = 0;
    M2([1:60, 400:end], :) = 0;
    maxVal = max(M2(:));
    [yval, xval] = find(M2 >= 0.98*maxVal);
    x2(j, 1) = mean(xval);
    y2(j, 1) = mean(yval);  
    %{
    pcolor(M2)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% rotate back by swopping xval and yval
for j = 1 : frameNum4(3)
    M3 = double(rgb2gray(vidFrames3_4(:, :, :, j)));
    M3(:, [1:300, 480:end]) = 0;
    M3([1:120, 320:end], :) = 0;
    maxVal = max(M3(:));
    [yval, xval] = find(M3 >= 0.92*maxVal);
    x3(j, 1) = mean(yval);
    y3(j, 1) = mean(xval);   
    %{
    pcolor(M3)
    shading interp
    drawnow
    colormap(gray)
    %}
end

% cut each data so that the can
% starts from the highest point
x1 = x1(7 : end);
y1 = y1(7 : end);
x2 = x2(5 : end);
y2 = y2(5 : end);
x3 = x3(1 : end);
y3 = y3(1 : end);

% organize data with same length
minlength = min([length(y1), length(y2), length(y3)]);
x1 = x1(1 : minlength);
y1 = y1(1 : minlength);
x2 = x2(1 : minlength);
y2 = y2(1 : minlength);
x3 = x3(1 : minlength);
y3 = y3(1 : minlength);

frame = 1:minlength;
Mx = [x1, y1, x2, y2, x3, y3];
Mx = Mx';

% examine for the path of paint can and appropriate
% threshold values
figure(1)
for i = 1 : 3
    subplot(3, 1, i)
    plot(frame, Mx(2 * i - 1, :), frame, Mx(2 * i, :))
    legend('horizontal position', 'vertical position', 'FontSize', 12, 'Location', 'bestoutside')
    xlabel('frame', 'FontSize', 12)
    ylabel('position in camera', 'FontSize', 12)
    title(['path (horizontal displacement & rotation) (camera',num2str(i),')'], 'FontSize', 12)
end

% subtract the mean
[m, n] = size(Mx);
mu = mean(Mx, 2);
Mx = Mx - repmat(mu, 1, n);

% perform PCA
[U, S, V] = svd(Mx / sqrt(n-1), 'econ');
Y = U' * Mx;

figure(2)
plot(frame, Y)
xlabel('frame', 'FontSize', 12)
ylabel('normalized position of paint can', 'FontSize', 12)
title('principal component in horizontal displacement and rotation case', 'FontSize', 12)
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5', 'mode 6', 'FontSize', 12, 'Location', 'bestoutside')
grid on

figure(3)
sig = diag(S);
plot(sig.^2/sum(sig.^2),'ro','Linewidth',2)
xticks([0:1:6])
xlabel('energy modes', 'FontSize', 12)
ylabel('energy percentage (%)', 'FontSize', 12)
title('energy of modes in horizontal displacement and rotation case', 'FontSize', 12)