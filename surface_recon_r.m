function coTemp = surface_recon_r(poses, bound_coords, frames, xoffset, yoffset, UStoCam)

% converts poses and segmentation coordinates into an nx3 matrix (coTemp)
% of points

coTemp = zeros(1000000,5); 
ind = 0;
for jj = frames 
    framePose = poses(:,:,jj);
    v = bound_coords{jj,1};
    u = bound_coords{jj,2};
    s = size(bound_coords);
    if s(2) == 3
        is = bound_coords{jj,3};
    end
    for ii = 1:numel(v) % de-densify for faster processing
        ind = ind + 1;
        xyz = framePose * UStoCam * [u(ii)-xoffset; v(ii)-yoffset; 1];
        xyz = xyz';
        coTemp(ind, 1:3) = xyz;
        if s(2) == 3
            coTemp(ind,5) = is(ii);
        end
    end
end

coTemp = coTemp(1:ind,:);
coTemp(:, 4) = 1:size(coTemp, 1);
