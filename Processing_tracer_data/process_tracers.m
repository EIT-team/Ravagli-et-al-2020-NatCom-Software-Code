
clc;
clear all;
close all;

%% Load tracer image and pre-process

A = imread('20190411 R common rev 1 5x merged.jpg');

A=A(350:1150,450:1650,:);   % Cut external black frame - this is image/dependent
A = flip(A,1);    % Flip according to point of view (up/down in microscope)

%% External perimeter - mark points and fit ellipse

% Show image with superimposed edges for choice of seed input
figure(1); hold on;
himage=imshow(A ,[]);
hold on;

np=15;
p = ginput(np);
for i=1:np
    seed_x(i) = round(axes2pix(size(A(:,:,1), 2), get(himage, 'XData'), p(i,2)));
    seed_y(i) = round(axes2pix(size(A(:,:,1), 1), get(himage, 'YData'), p(i,1)));
end

plot(seed_y,seed_x,'w*');

% Find the least squares geometric estimate
[zg, ag, bg, alphag] = fitellipse([seed_x; seed_y]);
plotellipse([zg(2) zg(1)], ag, bg, pi/2-alphag, 'w--');


%% Rotate to align one ellipse axis and resize to circle

A_rot=imrotate(A,rad2deg(-pi-alphag),'nearest','crop');
figure(2);imshow(A_rot);

[h w ~] = size(A_rot);
w = w*ag/bg;
A_rescaled = imresize(A_rot, [h w], 'nearest');

figure(3);
A_rescaled=A_rescaled(150:end-150,:,:); % chop rescaled image
imshow(A_rescaled);


%% Rotate to align cuff opening (marked by thread) as mesh/paraview

% Alignment is manual according to thread position in the image
A_circle=imrotate(A_rescaled,-150,'nearest','crop');

%% Find fascicles location

figure(4); hold on;
himage=imshow(A_circle ,[]);
hold on;

% Tibial fascicle
np=15;
p = ginput(np);
for i=1:np
    seed_x(i) = round(axes2pix(size(A_circle(:,:,1), 2), get(himage, 'XData'), p(i,2)));
    seed_y(i) = round(axes2pix(size(A_circle(:,:,1), 1), get(himage, 'YData'), p(i,1)));
end
plot(seed_y,seed_x,'w*');

[c_T, ax1_T, ax2_T, alpha_T] = fitellipse([seed_x; seed_y]);
plotellipse([c_T(2) c_T(1)], ax1_T, ax2_T, pi/2-alpha_T, 'w--');

% Sural fascicle
np=15;
p = ginput(np);
for i=1:np
    seed_x(i) = round(axes2pix(size(A_circle(:,:,1), 2), get(himage, 'XData'), p(i,2)));
    seed_y(i) = round(axes2pix(size(A_circle(:,:,1), 1), get(himage, 'YData'), p(i,1)));
end
plot(seed_y,seed_x,'w*');

[c_S, ax1_S, ax2_S, alpha_S] = fitellipse([seed_x; seed_y]);
plotellipse([c_S(2) c_S(1)], ax1_S, ax2_S, pi/2-alpha_S, 'w--');

% Peroneal fascicle
np=15;
p = ginput(np);
for i=1:np
    seed_x(i) = round(axes2pix(size(A_circle(:,:,1), 2), get(himage, 'XData'), p(i,2)));
    seed_y(i) = round(axes2pix(size(A_circle(:,:,1), 1), get(himage, 'YData'), p(i,1)));
end
plot(seed_y,seed_x,'w*');

[c_P, ax1_P, ax2_P, alpha_P] = fitellipse([seed_x; seed_y]);
plotellipse([c_P(2) c_P(1)], ax1_P, ax2_P, pi/2-alpha_P, 'w--');

% External circle (boundary of nerve)
np=15;
p = ginput(np);
for i=1:np
    seed_x(i) = round(axes2pix(size(A_circle(:,:,1), 2), get(himage, 'XData'), p(i,2)));
    seed_y(i) = round(axes2pix(size(A_circle(:,:,1), 1), get(himage, 'YData'), p(i,1)));
end
plot(seed_y,seed_x,'w*');

[c_ext, ax1_ext, ax2_ext, alpha_ext] = fitellipse([seed_x; seed_y]);
plotellipse([c_ext(2) c_ext(1)], ax1_ext, ax2_ext, pi/2-alpha_ext, 'w--');

% Plot centres of images
plot(c_T(2),c_T(1),'w.','MarkerSize',12)
plot(c_S(2),c_S(1),'w.','MarkerSize',12)
plot(c_P(2),c_P(1),'w.','MarkerSize',12)


%% Convert to proper values for external use

x0=c_ext(2);
y0=-c_ext(1);

com_tracers=[
    c_T(2)-x0 -c_T(1)-y0
    c_S(2)-x0 -c_S(1)-y0
    c_P(2)-x0 -c_P(1)-y0
];

com_tracers=com_tracers*700/mean([ax1_ext ax2_ext]);

figure; hold on; grid on; xlim([-800 800]); ylim([-800 800]); axis equal;
% Plot circle representing external boundaries
R=700;
theta_plot=[0:0.001:1]*2*pi;
x2_circle=R*cos(theta_plot);
x3_circle=R*sin(theta_plot);
plot(x2_circle,x3_circle,'k.');
plot(com_tracers(1,1),com_tracers(1,2),'r*')
plot(com_tracers(2,1),com_tracers(2,2),'bs')
plot(com_tracers(3,1),com_tracers(3,2),'g^')
legend('Ext','T','S','P')


































