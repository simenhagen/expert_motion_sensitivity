clear all
clc


%%%%%%%%%%%%%%%%%%%%%%%%
% load desire mat file
%%%%%%%%%%%%%%%%%%%%%%%%
%load eagle_2_merged.mat
%load Walk.mat
%load crane_1_merged_edited_beak.mat
%load eagle_2_merged_edited_beak.mat
%load duck_1_merged_edited_beak.mat
load crane_2_merged_edited.mat

%%%%%%%%%%%%%%%%%%%%%%%%
% set parameters
%%%%%%%%%%%%%%%%%%%%%%%%
scr = 0;            % scramble target or not
ndots = 0;         % number of desired noise dots
maxdots = 3000;     % maximum number of dots; see below

%%%%%%%%%%%%%%%%%%%%%%%%
% if bird, process
%%%%%%%%%%%%%%%%%%%%%%%%
if exist('dat')
    % remove occluded dots (replace with NaN)
    szx = size(mov,2);
    for i = 1:size(mov,4)
        idx = find(dat(:,1,i) < 20 & dat(:,2,i) > (szx-20));
        if ~isempty(idx)
            dat(idx,:,i) = NaN;
        end
    end
    span = 5;
    dats = smooth_points(dat,span); % smooth dots (moving average)
    datt = remove_translation(dats,1); % remove translation of smoothed dots
    mat = datt;
    rotflag = 0;
else
    rotflag = 1;        % rotate or not
    rotdeg = 235;   % rotation degree
    whichaxis = 2;  % y-axis of rotation    CHANGE DEPENDING ON orientationindex == 2 (INVERSION)
    plotflag = 0;   % show dots
    mat = rotxyz(mat,rotdeg,whichaxis,plotflag);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% make sure matrix is N x 2 x FRAMES
%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(mat);
mattmp = zeros(dim(1),2,dim(3));
for n = 1:dim(3)
    mattmp(:,:,n) = mat(:,1:2,n);
end
mat = mattmp;
dim = size(mat);

%%%%%%%%%%%%%%%%%%%%%%%%
% distribute evenly as possible in an arbitrary 3 x 3 grid
%%%%%%%%%%%%%%%%%%%%%%%%
xx = reshape(squeeze((mat(:,1,:))),dim(1)*dim(3),1);
yy = reshape(squeeze((mat(:,2,:))),dim(1)*dim(3),1);
xmn = floor(min(xx)); xmx = floor(max(xx)); ymn = floor(min(yy)); ymx = floor(max(yy));
xs = floor((xmx-xmn)/3); ys = floor((ymx-ymn)/3);
recti = [];
for i = 1:3
    for j = 1:3
        recti = [ recti; xmn+(i-1)*xs+1 xmn+(i*xs) ymn+(j-1)*ys+1 ymn+(j*ys) ];        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% create an arbitrarily large set of dots (maxdots)
%%%%%%%%%%%%%%%%%%%%%%%%
scramblevect_signal = [];
grid_index = [];
g = 1; % keep track of grid position
for n = 1:maxdots 
    scramblevect_signal = [ scramblevect_signal; randi([recti(g,1),recti(g,2)],[1,1]) randi([recti(g,3),recti(g,4)],[1,1]) ];
    grid_index = [ grid_index; g ];
    g = g + 1;
    if g > 9, g = 1; end
end
temp_index = randi(dim(1),maxdots,1);
temp_noise = mat(temp_index,:,:);
temp_scramble = [];
for n = 1:maxdots
    temp_scramble = [ temp_scramble; scramblevect_signal(n,:) ];
end
temp_scramble = repmat(temp_scramble,[1,1,dim(3)]);
tmp = temp_noise - repmat(temp_noise(:,:,1),[1,1,dim(3)]);
matx = tmp + temp_scramble;

% throw out dots that goes beyond border by jit amount, can change jit
% amount
checkdots = ones(maxdots,1);
jit = 25;
for i = 1:maxdots
    if any(matx(i,1,:)<(xmn-jit)) | any(matx(i,1,:)>(xmx+jit)) | any(matx(i,2,:)<(ymn-jit)) | any(matx(i,2,:)>(ymx+jit))
        checkdots(i) = 0;
    end
end
idxo = find(checkdots==1);

% for checking final distribution
check_params = [
    temp_index(idxo(ndots+1:ndots+dim(1))) grid_index(idxo(ndots+1:ndots+dim(1))); ...
    temp_index(idxo(1:ndots)) grid_index(idxo(1:ndots)) ];

mats = matx(idxo(ndots+1:ndots+dim(1)),:,:);    %%% FINAL SCR MATRIX
matx = matx(idxo(1:ndots),:,:);                 %%% FINAL NOISE MATRIX

%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%
close all; figure(1);
for loop = 1:3
    for n = 1:dim(3)
        if ~scr
            plot(mat(:,1,n),mat(:,2,n),'r.'); hold on;
        else
            plot(mats(:,1,n),mats(:,2,n),'r.'); hold on;
        end
        for xx = [xmn:xs:xmx]
            for yy = [ymn:ys:ymx]
                line([xx xx],[ymn ymx],'linestyle',':','color','k');
                line([xmn xmx],[yy yy],'linestyle',':','color','k');                
            end
        end
        line([xmn-jit xmn-jit],[ymn-jit ymx+jit],'linestyle',':','color','r');
        line([xmx+jit xmx+jit],[ymn-jit ymx+jit],'linestyle',':','color','r');
        line([xmn-jit xmx+jit],[ymn-jit ymn-jit],'linestyle',':','color','r');
        line([xmn-jit xmx+jit],[ymx+jit ymx+jit],'linestyle',':','color','r');
        plot(matx(:,1,n),matx(:,2,n),'k.'); hold off;
        %axis([-250 250 -400 400]);
        axis([xmn-50 xmx+50 ymn-50 ymx+50]);        
        if exist('dat'), axis ij, end;
        if loop==1 & n==1
            title('press space to start');
            pause
        else
            title(sprintf('loop %i',loop));
        end
        pause(.1);        
    end
end

% check
figure; hist(check_params(:,1),dim(1)); title('N per target dot');
figure; hist(check_params(:,2),9); title('N per grid position');

return;


