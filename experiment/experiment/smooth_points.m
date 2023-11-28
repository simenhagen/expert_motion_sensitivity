function dat = smooth_points(dat,span)

% smooth point data using moving average of span size
% 21/05/16, qcv


% replace occluded points with NaN before passing dat in

sz = size(dat);

dats = zeros(sz);


for i = 1:sz(1)             % n points
    xx = squeeze(dat(i,1,:));
    yy = squeeze(dat(i,2,:));
    xx = smooth(xx,span,'moving');
    yy = smooth(yy,span,'moving');
    dats(i,1,:) = xx; dats(i,2,:) = yy;
end


dat = dats;


