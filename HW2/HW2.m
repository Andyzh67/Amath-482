clear all; close all; clc;

% Part 1

load handel
v = y';
%{
figure(1)
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');
hold off
%}

% prepare data
L = 9;
v = v(1 : length(v) - 1);
n = length(v);
t = (1:length(v)) / Fs;
k= (2*pi/L)*[0:(n/2-1), -n/2:-1];
ks = fftshift(k);

% create Gabor Filters

% explore on the effect of different gabor window width (delta t of 0.1 is 
% applied here, set to be the controlled variable)

width = [1, 5, 25, 50];
tslide = 0 : 0.1: L;

for i = 1 : 4
    vgt_spec = [];
    for j = 1 : length(tslide)
        gabor = exp(-width(i) * (t - tslide(j)).^2);
        vg = gabor.* v;
        vgt = fft(vg);
        vgt_spec = [vgt_spec; abs(fftshift(vgt))];
    end
    figure(2)
    subplot(2, 2, i)
    pcolor(tslide, ks, vgt_spec.')
    shading interp
    colormap(hot)
    xlabel('Time in second', 'FontSize', 12)
    ylabel('Frequency(\omega)', 'FontSize', 12)
    title(['spectrogram created with fixed dt and Gaussian window width of  ', num2str(width(i))], 'FontSize', 10)
    xticks([0:1:L])
    yticks([-25000:5000:25000])
end


% explore on the idea of oversampling and undersampling, that is, applying
% different translations of gabor window (window width is fixed in this
% time as a controlled variable, an appropriate width of 25 is applied)


width = 25;
tslide = [];
dt = [0.05, 0.1, 0.5, 1.5];
tslide{1} = 0:0.05:L;
tslide{2} = 0:0.1:L;
tslide{3} = 0:0.5:L;
tslide{4} = 0:1.5:L;
for i = 1 : 4
    vgt_spec = [];
    for j = 1 : length(tslide{i})
        gabor = exp(-width * (t - tslide{i}(j)).^2);
        vg = gabor.* v;
        vgt = fft(vg);
        vgt_spec = [vgt_spec; abs(fftshift(vgt))];
    end
    figure(3)
    subplot(2, 2, i)
    pcolor(tslide{i}, ks, vgt_spec.')
    shading interp
    colormap(hot)
    xlabel('Time in second', 'FontSize', 12)
    ylabel('Frequency(\omega)', 'FontSize', 12)
    title(['spectrogram created with fixed Gaussian window width and dt of  ', num2str(tslide{i}(2) - tslide{i}(1))], 'FontSize', 10)
    xticks([0:1:L])
    yticks([-25000:5000:25000])
end



% explore on different gabor windows(Gaussian, Mexican hat, and Shannon)

% Gaussian
figure(4)
vgt_spec = [];
tslide = 0 : 0.1 : L;
for j = 1 : length(tslide)
    gabor = exp(-25 * (t - tslide(j)).^2);
    vg = gabor.* v;
    vgt = fft(vg);
    vgt_spec = [vgt_spec; abs(fftshift(vgt))];
end
subplot(3, 1, 1)
pcolor(tslide, ks, vgt_spec.')
shading interp
colormap(hot)
xlabel('Time in second', 'FontSize', 12)
ylabel('Frequency(\omega)', 'FontSize', 12)
title('spectrogram created using Gaussian filter', 'FontSize', 12)
xticks([0:1:L])
yticks([-25000:5000:25000])

% Mexican hat
vmt_spec = [];
width = 0.25;
tslide = 0:0.1:L;
for j = 1 : length(tslide)
    mexhat = 2 / (sqrt(3 * width) * pi^(1/4)) * (1 - ((t - tslide(j)) / width).^2)...
    .* exp(-((t - tslide(j)).^2) / (2 * width^2));
    vm = mexhat.* v;
    vmt = fft(vm);
    vmt_spec = [vmt_spec; abs(fftshift(vmt))];
end
subplot(3, 1, 2)
pcolor(tslide, ks, vmt_spec.')
shading interp
colormap(hot)
xlabel('Time in second', 'FontSize', 12)
ylabel('Frequency(\omega)', 'FontSize', 12)
title('spectrogram created using Mexican hat filter', 'FontSize', 12)
xticks([0:1:L])
yticks([-25000:5000:25000])

% Shannon
vst_spec = [];
width = 0.4;
for j = 1 : length(tslide)
    shannon = (abs(t - tslide(j)) < width);
    vs = shannon.* v;
    vst = fft(vs);
    vst_spec = [vst_spec; abs(fftshift(vst))];
end
subplot(3, 1, 3)
pcolor(tslide, ks, vst_spec.')
shading interp
colormap(hot)
xlabel('Time in second', 'FontSize', 12)
ylabel('Frequency(\omega)', 'FontSize', 12)
title('spectrogram created using Shannon filter', 'FontSize', 12)
xticks([0:1:L])
yticks([-25000:5000:25000])



% plot the window with original signal at the center
tau = 4.5;
gabor = exp(-25 * (t - tau).^2);
mexhat = 2 / (sqrt(3 * 0.25) * pi^(1/4)) * (1 - ((t - tau) / 0.25).^2)...
    .* exp(-((t - tau).^2) / (2 * 0.25^2));
shannon = (abs(t - tau) < 0.4);
figure(5)
subplot(3,1,1)
plot(t, v, 'k')
hold on
plot(t, gabor, 'Linewidth', 3)
xlabel('Time in second', 'FontSize', 12)
ylabel('Frequency(\omega)', 'FontSize', 12)
title('centered Gaussian window with width of 25', 'FontSize', 12)

subplot(3,1,2)
plot(t, v, 'k')
hold on
plot(t, mexhat, 'Linewidth', 3)
xlabel('Time in second', 'FontSize', 12)
ylabel('Frequency(\omega)', 'FontSize', 12)
title('centered Mexican hat window with width of 0.25', 'FontSize', 12)

subplot(3,1,3)
plot(t, v, 'k')
hold on
plot(t, shannon, 'Linewidth', 3)
xlabel('Time in second', 'FontSize', 12)
ylabel('Frequency(\omega)', 'FontSize', 12)
title('centered Shannon window with width of 0.4', 'FontSize', 12)



%%
% Part 2

% Piano
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
%{
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');
%}

L = 16;
v = y';
n = length(v);
t = (1 : length(v)) / Fs;
k= (1/L)*[0:(n/2-1), -n/2:-1];
ks = fftshift(k);
tslide = 0 : 0.1 : L;
vgt_spec = [];
frequency = [];

width = 30;
for j = 1 : length(tslide)
    gabor = exp(-width * (t - tslide(j)).^2);
    vg = v.* gabor;
    vgt = fft(vg);
    [vmax, index] = max(abs(vgt));
    frequency = [frequency; abs(k(index))];
    vgt_spec = [vgt_spec; abs(fftshift(vgt))];
end
figure(6)
plot(tslide, frequency)
xlabel('time in second', 'FontSize', 15)
ylabel('frequency in hertz', 'FontSize', 15)
title('the center frequencies of the piano', 'FontSize', 15)


figure(7)
pcolor(tslide, ks, vgt_spec.')
shading interp
colormap(hot)
xlabel('time in second', 'FontSize', 15)
ylabel('frequency(\omega)', 'FontSize', 15)
title('the spectrogram of the piano', 'FontSize', 15)
ylim([0, 1000])
xticks([0:1:L])

%%
% recorder
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
%{
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
%}

L = 14.2;
v = y';
v = v(1 : length(v));
n = length(v);
t = (1 : length(v)) / Fs;
k= (1/L)*[0:(n/2-1), -n/2:-1];
ks = fftshift(k);
tslide = 0 : 0.1 : L;
vgt_spec = [];
frequency = [];

width = 30;
for j = 1 : length(tslide)
    gabor = exp(-width * (t - tslide(j)).^2);
    vg = v.* gabor;
    vgt = fft(vg);
    [vmax, index] = max(abs(vgt));
    frequency = [frequency; abs(k(index))];
    vgt_spec = [vgt_spec; abs(fftshift(vgt))];
end
figure(8)
plot(tslide, frequency)
xlabel('time in second', 'FontSize', 15)
ylabel('frequency in hertz', 'FontSize', 15)
title('the center frequencies of the recorder', 'FontSize', 15)


figure(9)
pcolor(tslide, ks, vgt_spec.')
shading interp
colormap(hot)
xlabel('time in second', 'FontSize', 15)
ylabel('frequency(\omega)', 'FontSize', 15)
title('the spectrogram of the recorder', 'FontSize', 15)
ylim([0, 3000])
xticks([0:1:L])
