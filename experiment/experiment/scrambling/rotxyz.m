function [xyz] = rotz(xyz,d,wha,plotflag)

% xyz = ndots x 3 x nframes
% wha = which axis of rotation, 1 = x, 2 = y, 3 = z;
% d = angle in degree. 0-360 deg; negative doesn't seem to work

% convert to radian
d = deg2rad(d);

% transformation matrix
mtx = zeros(3,3,3);
mtx(:,:,1) = [1 0 0; 0 cos(d) -sin(d); 0 sin(d) cos(d)];
mtx(:,:,2) = [cos(d) 0 sin(d); 0 1 0; -sin(d) 0 cos(d)];
mtx(:,:,3) = [cos(d) -sin(d) 0; sin(d) cos(d) 0; 0 0 1];

dim = size(xyz);

xyzr = zeros(dim);
for i = 1:dim(3)
    for j = 1:dim(1)
        xyzr(j,:,i) = mtx(:,:,wha)*xyz(j,:,i)';
    end
end

if plotflag
    close all; figure(1);
    for n = 1:dim(3)
        subplot(1,2,1);
        plot(xyzr(:,1,n),xyzr(:,2,n),'k.');
        axis([-200 200 -200 200]);
        title('2D View');
        subplot(1,2,2);
        plot3(xyz(:,1,n),xyz(:,3,n),xyz(:,2,n),'k.');
        axis([-200 200 -200 200 -200 200]);
        title('3D View');
        pause(.1);
    end
end

xyz = xyzr;


