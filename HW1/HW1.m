clear all; close all; clc;

load('Testdata.mat')

L=15; % spatial domain
n=64; % Fourier modes

x2=linspace(-L,L,n+1); 
x=x2(1:n); 
y=x; 
z=x;


k=(2*pi/(2*L))*[0:(n/2-1), -n/2:-1]; 
ks=fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

uave = zeros(n, n, n);

% reshape the data into 64*64*64 matrix, add them up
for j = 1:20
    u(:,:,:)=reshape(Undata(j,:),n,n,n);
    ut = fftn(u);
    uave = uave + ut;
end

% take average and and normalize it
uave = abs(fftshift(uave)) / j;
uave = uave / max(uave(:));

% find the position (Px, Py, Pz) of aveMax in uave
[aveMax, Index] = max(uave(:));
[Px,Py,Pz] = ind2sub(size(uave), Index);

% find the corresponding position (xf, yf, zf) of aveMax in frequency space
xf = Kx(Px,Py,Pz);
yf = Ky(Px,Py,Pz);
zf = Kz(Px,Py,Pz);
    
% make sure we find the correct frequencies xf, yf, zf    
figure(1)
isosurface(Kx,Ky,Kz,uave,0.4)
title('Isosurface of frequency space with isovalue of 0.4','FontSize',15)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('z','FontSize',18)
axis([-10 10 -10 10 -10 10])
grid on

figure(2)
isosurface(Kx,Ky,Kz,uave,0.9)
title('Isosurface of frequency space with isovalue of 0.9','FontSize',15)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('z','FontSize',18)
axis([-2 2 -2 2 -2 2])
grid on


%for test purpose only
%{
tau = [0.1, 0.3, 0.9];
for k = 1 : 3
    filter = exp(-tau(k) * ((Kx - xf).^2 + (Ky - yf).^2 + (Kz - zf).^2));

    trajectory = zeros(20, 3);
    for i = 1 : 20
        u(:,:,:)=reshape(Undata(i,:),n,n,n);
        ut = fftshift(fftn(u));
        utf = ifftshift(filter.* ut);
        uf = ifftn(utf);
        [traceMax, Index] = max(uf(:));
        [xt, yt, zt] = ind2sub(size(uf), Index);
        trajectory(i, 1) = X(xt, yt, zt);
        trajectory(i, 2) = Y(xt, yt, zt);
        trajectory(i, 3) = Z(xt, yt, zt);
    end
    plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'Linewidth', 2)
    hold on
end
title('test on appropriate tau value', 'FontSize', 15)
xlabel('x', 'FontSize', 18)
ylabel('y', 'FontSize', 18)
zlabel('z', 'FontSize', 18)
legend('tau = 0.1', 'tau = 0.3', 'tau = 0.9', 'FontSize', 15)
%}

% build a filter and center on data around (xf, yf, zf)
filter = exp(-0.3 * ((Kx - xf).^2 + (Ky - yf).^2 + (Kz - zf).^2));

% trace the trajectory of the marble through each data spatially
trajectory = zeros(20, 3);
for i = 1 : 20
    u(:,:,:)=reshape(Undata(i,:),n,n,n);
    ut = fftshift(fftn(u));
    utf = ifftshift(filter.* ut);
    uf = ifftn(utf);
    [traceMax, Index] = max(uf(:));
    [xt, yt, zt] = ind2sub(size(uf), Index);
    trajectory(i, 1) = X(xt, yt, zt);
    trajectory(i, 2) = Y(xt, yt, zt);
    trajectory(i, 3) = Z(xt, yt, zt);
end

figure(3)
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'Linewidth', 2)
title('trajectory of the marble in spatial domain', 'FontSize', 15)
hold on
plot3(trajectory(20, 1), trajectory(20, 2), trajectory(20, 3), 'o', 'MarkerSize', 10)
xlabel('x', 'FontSize', 18)
ylabel('y', 'FontSize', 18)
zlabel('z', 'FontSize', 18)