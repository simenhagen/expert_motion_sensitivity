function dat = remove_translation(dat,ref)

% ref    0 = use center of mass, 1 to npoints = use that point


% remove translation
% 21/05/16, qcv


% replace occluded points with NaN before passing dat in

sz = size(dat);
datt = zeros(sz);
datt(:,1,1) = dat(:,1,1);
datt(:,2,1) = dat(:,2,1);

% get reference point (frame 1)
if ref==0
    xr = mean(squeeze(dat(:,1,1)));
    yr = mean(squeeze(dat(:,2,1)));
else
    xr = dat(ref,1,1);
    yr = dat(ref,2,1);
end


% loop thru frame, subtract translation
for i = 2:sz(3)
    if ref==0
        xx = mean(squeeze(dat(:,1,i))) - xr;
        yy = mean(squeeze(dat(:,2,i))) - yr;
    else
        xx = dat(ref,1,i) - xr;
        yy = dat(ref,2,i) - yr;
    end
    datt(:,1,i) = dat(:,1,i) - xx; 
    datt(:,2,i) = dat(:,2,i) - yy; 
end


dat = datt;
