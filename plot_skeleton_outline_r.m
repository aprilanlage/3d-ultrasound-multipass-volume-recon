function plot_skeleton_outline_r(poses_downsampled, bound_coords, frames, xoffset, yoffset, UStoCam)

% plots segmentation outline of specified frames 

coTemp = zeros(1000000,3); 
ind = 0;
for jj = frames 
    framePose = poses_downsampled(:,:,jj);
    v = bound_coords{jj,1};
    u = bound_coords{jj,2};

    for ii = 1:numel(v) % de-densify for faster processing
        ind = ind + 1;
        xyz = framePose * UStoCam * [u(ii)-xoffset; v(ii)-yoffset; 1];
        xyz = xyz';
        coTemp(ind, :) = xyz;
    end
end

coTemp = coTemp(1:ind,:);
coTemp(:, 4) = 1:size(coTemp, 1);

figure
plot = scatter3(coTemp(:, 1), coTemp(:, 2), coTemp(:, 3), 0.5, parula(size(coTemp, 1)));
axis equal;
xlabel('x'); ylabel('y'); zlabel('z')
