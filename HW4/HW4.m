%%
% Test 1: Band Classification (different genres)
clear; close all; clc;

% Dan Bodan
x1 = audioread('Anton1.wav');
x2 = audioread('Anton2.wav');
x3 = audioread('Anton4.wav');
x4 = audioread('City4.wav');
x5 = audioread('City2.wav');
x6 = audioread('City3.wav');
x7 = audioread('National1.wav');
x8 = audioread('National2.wav');
x9 = audioread('National3.wav');

x = [x1, x2, x3, x4, x5, x6, x7, x8, x9];
for i = 1 : 9
    x(:, i) = (x(:, 2*i - 1) + x(:, 2*i)) / 2;
end
x = x(:, 1 : end/2);



% Density & Time
y1 = audioread('Clean1.wav');
y2 = audioread('Clean2.wav');
y3 = audioread('Clean3.wav');
y4 = audioread('Dream4.wav');
y5 = audioread('Dream2.wav');
y6 = audioread('Dream3.wav');
y7 = audioread('Day4.wav');
y8 = audioread('Day1.wav');
y9 = audioread('Day3.wav');

y = [y1, y2, y3, y4, y5, y6, y7, y8, y9];
for i = 1 : 9
    y(:, i) = (y(:, 2*i - 1) + y(:, 2*i)) / 2;
end
y = y(:, 1 : end/2);

% Dan Lebowitz
z1 = audioread('Rag3.wav');
z2 = audioread('Rag1.wav');
z3 = audioread('Rag4.wav');
z4 = audioread('Opening4.wav');
z5 = audioread('Opening2.wav');
z6 = audioread('Opening1.wav');
z7 = audioread('Birds1.wav');
z8 = audioread('Birds2.wav');
z9 = audioread('Birds3.wav');

z = [z1, z2, z3, z4, z5, z6, z7, z8, z9];
for i = 1 : 9
    z(:, i) = (z(:, 2*i - 1) + z(:, 2*i)) / 2;
end
z = z(:, 1 : end/2);

Test1 = [x, y, z];
spec = doSpec(Test1);
feature = 8;

[U, S, V, threshold1, threshold2, w, M1, M2, M3] = trainer(spec, feature);

file(:, 1) = prepareTest('Anton3.wav');
file(:, 2) = prepareTest('City1.wav');
file(:, 3) = prepareTest('National4.wav');
file(:, 4) = prepareTest('Kallimachou1.wav');
file(:, 5) = prepareTest('Kallimachou2.wav');
file(:, 6) = prepareTest('Clean4.wav');
file(:, 7) = prepareTest('Dream1.wav');
file(:, 8) = prepareTest('Day2.wav');
file(:, 9) = prepareTest('Castlevania1.wav');
file(:, 10) = prepareTest('Castlevania2.wav');
file(:, 11) = prepareTest('Rag2.wav');
file(:, 12) = prepareTest('Opening3.wav');
file(:, 13) = prepareTest('Birds4.wav');
file(:, 14) = prepareTest('Cliffsides1.wav');
file(:, 15) = prepareTest('Cliffsides2.wav');


spec = doSpec(file);
[label, pval] = check(U, w, threshold1, threshold2, spec);
sucRate = checkRate(M1, M2, M3, label);


%%
% Test 2: Band Classification (same genre, Jazz & Blues)
clear all; close all; clc;

% Quincas Moreira
file(:, 1) = prepareTest('Brooklin1.wav');
file(:, 2) = prepareTest('Brooklin2.wav');
file(:, 3) = prepareTest('Brooklin3.wav');
file(:, 4) = prepareTest('Brooklin4.wav');
file(:, 5) = prepareTest('Infusion1.wav');
file(:, 6) = prepareTest('Infusion2.wav');
file(:, 7) = prepareTest('Infusion3.wav');
file(:, 8) = prepareTest('Infusion4.wav');
file(:, 9) = prepareTest('Malan1.wav');
file(:, 10) = prepareTest('Malan2.wav');
file(:, 11) = prepareTest('Malan3.wav');
file(:, 12) = prepareTest('Malan4.wav');
file(:, 13) = prepareTest('Paradox1.wav');
file(:, 14) = prepareTest('Paradox2.wav');
file(:, 15) = prepareTest('Paradox3.wav');
file(:, 16) = prepareTest('Paradox4.wav');


i = 16;
% Chris Haugen
file(:, 1+i) = prepareTest('Bleeker1.wav');
file(:, 2+i) = prepareTest('Bleeker2.wav');
file(:, 3+i) = prepareTest('Bleeker3.wav');
file(:, 4+i) = prepareTest('Bleeker4.wav');
file(:, 5+i) = prepareTest('Et1.wav');
file(:, 6+i) = prepareTest('Et2.wav');
file(:, 7+i) = prepareTest('Et3.wav');
file(:, 8+i) = prepareTest('Et4.wav');
file(:, 9+i) = prepareTest('Old1.wav');
file(:, 10+i) = prepareTest('Old2.wav');
file(:, 11+i) = prepareTest('Old3.wav');
file(:, 12+i) = prepareTest('Old4.wav');
file(:, 13+i) = prepareTest('Sunshine1.wav');
file(:, 14+i) = prepareTest('Sunshine2.wav');
file(:, 15+i) = prepareTest('Sunshine3.wav');
file(:, 16+i) = prepareTest('Sunshine4.wav');


% Freedom trauk studio
i = 32;
file(:, 1+i) = prepareTest('Crazy1.wav');
file(:, 2+i) = prepareTest('Crazy2.wav');
file(:, 3+i) = prepareTest('Crazy3.wav');
file(:, 4+i) = prepareTest('Crazy4.wav');
file(:, 5+i) = prepareTest('Current1.wav');
file(:, 6+i) = prepareTest('Current2.wav');
file(:, 7+i) = prepareTest('Current3.wav');
file(:, 8+i) = prepareTest('Current4.wav');
file(:, 9+i) = prepareTest('Elder1.wav');
file(:, 10+i) = prepareTest('Elder2.wav');
file(:, 11+i) = prepareTest('Elder3.wav');
file(:, 12+i) = prepareTest('Elder4.wav');
file(:, 13+i) = prepareTest('Mix1.wav');
file(:, 14+i) = prepareTest('Mix2.wav');
file(:, 15+i) = prepareTest('Mix3.wav');
file(:, 16+i) = prepareTest('Mix4.wav');


spec = doSpec(file);
feature = 20;

[U, S, V, threshold1, threshold2, w, M1, M2, M3] = trainer(spec, feature);

file2(:, 1) = prepareTest('Brooklin5.wav');
file2(:, 2) = prepareTest('Infusion5.wav');
file2(:, 3) = prepareTest('Malan5.wav');
file2(:, 4) = prepareTest('Paradox5.wav');
file2(:, 5) = prepareTest('Miles1.wav');
file2(:, 6) = prepareTest('Miles2.wav');

file2(:, 7) = prepareTest('Bleeker5.wav');
file2(:, 8) = prepareTest('Et5.wav');
file2(:, 9) = prepareTest('Front1.wav');
file2(:, 10) = prepareTest('Front2.wav');
file2(:, 11) = prepareTest('Old5.wav');
file2(:, 12) = prepareTest('Sunshine5.wav');

file2(:, 13) = prepareTest('Crazy5.wav');
file2(:, 14) = prepareTest('Current5.wav');
file2(:, 15) = prepareTest('Elder5.wav');
file2(:, 16) = prepareTest('Mix5.wav');
file2(:, 17) = prepareTest('Swing1.wav');
file2(:, 18) = prepareTest('Swing2.wav');

spec2 = doSpec(file2);
[label, pval] = check(U, w, threshold1, threshold2, spec2);
sucRate = checkRate(M1, M2, M3, label);





%%
% Test 3: Genre Classification
clear; close all; clc;

% Jazz & Blues
file(:, 1) = prepareTest('Et5.wav');
file(:, 2) = prepareTest('Miles2.wav');
file(:, 3) = prepareTest('Malan2.wav');
file(:, 4) = prepareTest('Malan3.wav');
file(:, 5) = prepareTest('Bleeker3.wav');
file(:, 6) = prepareTest('Swing2.wav');
file(:, 7) = prepareTest('Current3.wav');
file(:, 8) = prepareTest('Front2.wav');
file(:, 9) = prepareTest('Sunshine1.wav');
file(:, 10) = prepareTest('Swing1.wav');

i = 10;

% Electronics
file(:, 1+i) = prepareTest('Day3.wav');
file(:, 2+i) = prepareTest('Dream4.wav');
file(:, 3+i) = prepareTest('Clean2.wav');
file(:, 4+i) = prepareTest('Macaw1.wav');
file(:, 5+i) = prepareTest('Dragonfly1.wav');
file(:, 6+i) = prepareTest('Stranger2.wav');
file(:, 7+i) = prepareTest('Tea1.wav');
file(:, 8+i) = prepareTest('Tea2.wav');
file(:, 9+i) = prepareTest('Clean4.wav');
file(:, 10+i) = prepareTest('Macaw2.wav');

i = 20;

% Punk
file(:, 1+i) = prepareTest('Beach.wav');
file(:, 2+i) = prepareTest('Fail1.wav');
file(:, 3+i) = prepareTest('Fail2.wav');
file(:, 4+i) = prepareTest('Girl1.wav');
file(:, 5+i) = prepareTest('Girl2.wav');
file(:, 6+i) = prepareTest('Runner2.wav');
file(:, 7+i) = prepareTest('Pumps1.wav');
file(:, 8+i) = prepareTest('Radio1.wav');
file(:, 9+i) = prepareTest('Radio2.wav');
file(:, 10+i) = prepareTest('Runner1.wav');

spec = doSpec(file);
feature = 8;

[U, S, V, threshold1, threshold2, w, M1, M2, M3] = trainer(spec, feature);

file2(:, 1) = prepareTest('Infusion2.wav');
file2(:, 2) = prepareTest('Paradox5.wav');
file2(:, 3) = prepareTest('Elder3.wav');
file2(:, 4) = prepareTest('Old4.wav');

file2(:, 5) = prepareTest('Firefly1.wav');
file2(:, 6) = prepareTest('Cubic.wav');
file2(:, 7) = prepareTest('Operatic2.wav');
file2(:, 8) = prepareTest('Anniversary2.wav');

file2(:, 9) = prepareTest('Boy.wav');
file2(:, 10) = prepareTest('Give1.wav');
file2(:, 11) = prepareTest('Grassy.wav');
file2(:, 12) = prepareTest('Outlet1.wav');


spec2 = doSpec(file2);
[label, pval] = check(U, w, threshold1, threshold2, spec2);
[sucRate, correctlabel] = checkRate(M1, M2, M3, label);




%%
function spec = doSpec(data)
    for k = 1 : size(data, 2)
        temp = spectrogram(data(:, k));
        temp = abs(temp);
        [m, n] = size(temp);
        spec(:, k) = reshape(temp, m * n, 1);
    end
end

function file = prepareTest(clip)
    temp = audioread(clip);
    file = (temp(:, 1) + temp(:, 2))/2;
end


function [sucRate, correctLabel] = checkRate(M1, M2, M3, label)
    Mean = [M1, M2, M3];
    Mean = sort(Mean);
    for i = 1 : 3
        if Mean(i) == M1
            ind1 = i;
        end
        if Mean(i) == M2
            ind2 = i;
        end
        if Mean(i) == M3
            ind3 = i;
        end
    end
    correctLabel = [];
    siz = length(label)/3;
    for i = 1 : siz
        correctLabel(i) = ind1;
        correctLabel(i + siz) = ind2;
        correctLabel(i + 2*siz) = ind3;
    end
    
    fail = 0;
    result = correctLabel - label;
    for i = 1 : length(result)
        if result(i) ~= 0
            fail = fail + 1;
        end
    end
    
    sucRate = 1 - fail/length(result);
end

function [label, pval] = check(U, w, threshold1, threshold2, spec)
    TestMat = U' * spec;
    pval = w' * TestMat;
    
    label = (pval > threshold1) + (pval > threshold2) + 1;
end

function [U, S, V, threshold1, threshold2, w, M1, M2, M3] = trainer(sp, feature)
[U, S, V] = svd(sp, 'econ');

PCA_proj = S * V';
U = U(:, 1 : feature);

plot(diag(S).^2/sum(diag(S).^2), 'ro')
xlabel('Number of modes', 'FontSize', 12)
ylabel('Energy captured by each mode', 'FontSize', 12)
title('Check energy components for dimension reduction', 'FontSize', 12)

Art1 = PCA_proj(1 : feature, 1 : end/3);
Art2 = PCA_proj(1 : feature, end/3 + 1 : end * 2/3);
Art3 = PCA_proj(1 : feature, end * 2/3 + 1 : end);

m1 = mean(Art1, 2);
m2 = mean(Art2, 2);
m3 = mean(Art3, 2);

Sw = 0;
for k = 1 : size(Art1, 2)
    Sw = Sw + (Art1(:, k) - m1) * (Art1(:, k) - m1)';
end
for k = 1 : size(Art2, 2)
    Sw = Sw + (Art2(:, k) - m2) * (Art2(:, k) - m2)';
end
for k = 1 : size(Art3, 2)
    Sw = Sw + (Art3(:, k) - m3) * (Art3(:, k) - m3)';
end

mu = (m1 + m2 + m3) / 3;
Sb = 0;
Sb = Sb + (m1 - mu) * (m1 - mu)';
Sb = Sb + (m2 - mu) * (m2 - mu)';
Sb = Sb + (m3 - mu) * (m3 - mu)';

[V2, D] = eig(Sb, Sw);
[~, ind1] = max(diag(D));

w = V2(:, ind1);
w = w / norm(w, 2);

vArt1 = w' * Art1;
vArt2 = w' * Art2;
vArt3 = w' * Art3;

M1 = mean(vArt1);
M2 = mean(vArt2);
M3 = mean(vArt3);

sortArt1 = sort(vArt1);
sortArt2 = sort(vArt2);
sortArt3 = sort(vArt3);
sortArt = [sortArt1; sortArt2; sortArt3];

Mean = [M1, M2, M3];
Mmax = max([M1, M2, M3]);

IndMax = Mean == Mmax;
sortMax = sortArt(IndMax, :);

sortArt = setdiff(sortArt, sortMax, 'rows');


sortNext = sortArt(2, :);
sortLeast = sortArt(1, :);


sizeMax = length(sortArt1);
t1 = sizeMax;
t2 = 1;
while sortLeast(t1) > sortNext(t2)
    t1 = t1 - 1;
    t2 = t2 + 1;
end
threshold1 = (sortLeast(t1) + sortNext(t2)) / 2;

t3 = sizeMax;
t4 = 1;
while sortNext(t3) > sortMax(t4)
    t3 = t1 - 1;
    t4 = t4 + 1;
end
threshold2 = (sortNext(t3) + sortMax(t4)) / 2;

end