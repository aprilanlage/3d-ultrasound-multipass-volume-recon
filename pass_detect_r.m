function out_frames_4 = pass_detect_r(poses_downsampled, bound_coords, xoffset, yoffset, UStoCam)

% uses x poses to finds peaks to return changes in direction
% works well for phantom scans
% may need to adjust some limits
% not sure about in vivo

x_poses = squeeze(poses_downsampled(1,4,:));
y_poses = squeeze(poses_downsampled(2,4,:));
z_poses = squeeze(poses_downsampled(3,4,:));

fra_downsampled = (1:length(poses_downsampled));

[p, l] = findpeaks((abs(x_poses) .* abs(y_poses)),"MinPeakProminence",0.0005);
%findpeaks(-x_poses,"MinPeakProminence",0.01)
[pn, ln] = findpeaks(-(abs(x_poses) .* abs(y_poses)),"MinPeakProminence",0.0005);

out_frames = sort(cat(1,l,ln));

% optional plotting
%scatter(fra_downsampled,x_poses)
%scatter(fra_downsampled,(abs(x_poses) .* abs(y_poses)));hold on
%scatter(fra_downsampled(ln),-pn,'ko','MarkerFaceColor','r');
%scatter(fra_downsampled(l),p,'ko','MarkerFaceColor','g');

%plot_skeleton_outline(poses_downsampled, bound_coords, 2930:3038, xoffset, yoffset, UStoCam)

% record segmented frames
seg_slices = zeros(length(bound_coords),1);
for o  = 1:length(bound_coords)
    if nnz(bound_coords{o,1}) > 0
        seg_slices(o) = 1;
    end
end

% check that passes are a min number of poses
out_frames_2 = out_frames;
for t = 1:(length(out_frames)-1)
    p_1 = out_frames(t);
    p_2 = out_frames(t+1);

    % check that passes are a min number of poses
    if (p_2 - p_1) < 30
        out_frames_2(t+1) = NaN;
        disp('Too few poses')

    % check number of segmented frames
    elseif sum(seg_slices(p_1:p_2)) < 20
        out_frames_2(t+1) = NaN;
        disp('Too few segmentations')
      
    else
        continue
    end

end

out_frames_3 = zeros(3*length(out_frames_2),1);
fr_index = 1;
for k = 1:(length(out_frames_2)-1)
    if isnan(out_frames_2(k+1)) || isnan(out_frames_2(k))
        continue
    else
        % check beginning and end of pass for off-center segmentations
        % calc centroid of middle of pass (middle 30 frames)
        reject_frames = zeros(100,1);
        mid_pt = out_frames_2(k) + round((out_frames_2(k+1)-out_frames_2(k))/2);
        coTemp = surface_recon(poses_downsampled, bound_coords, (mid_pt-15):(mid_pt+15), xoffset, yoffset, UStoCam);
        mid_cen = [mean(coTemp(:, 1)), mean(coTemp(:, 2)), mean(coTemp(:, 3))];
        reject_i = 1;
    
        for fr = out_frames_2(k):mid_pt
            coTemp = surface_recon(poses_downsampled, bound_coords, fr, xoffset, yoffset, UStoCam);
            cen = [mean(coTemp(:, 1)), mean(coTemp(:, 2)), mean(coTemp(:, 3))];
            if (norm(mid_cen-cen) > 0.05)
                reject_frames(reject_i) = fr;
                reject_i = reject_i + 1;
            end
        end
    
        split_i = reject_i;
    
        for fr = mid_pt:out_frames_2(k+1)
            coTemp = surface_recon(poses_downsampled, bound_coords, fr, xoffset, yoffset, UStoCam);
            cen = [mean(coTemp(:, 1)), mean(coTemp(:, 2)), mean(coTemp(:, 3))];
            if (norm(mid_cen-cen) > 0.03)
                reject_frames(reject_i) = fr;
                reject_i = reject_i + 1;
            end
        end
    
        % make sure beginning frames were flagged
        if split_i ~= 1
            new_min = reject_frames(split_i-1) + 1;
        else
            new_min = out_frames_2(k);
        end
        
        % make sure end frames were flagged
        if split_i ~= reject_i
            new_max = reject_frames(split_i) - 1;
        else
            new_max = out_frames_2(k+1);
        end
    
        out_frames_3(fr_index) = new_min;
        out_frames_3(fr_index + 1) = new_max;
        out_frames_3(fr_index + 2) = NaN;
        fr_index = fr_index + 3;
    end
end

out_frames_3 = out_frames_3(1:fr_index-2);

out_frames_4 = zeros(length(out_frames_3),1);
ind_4 = 1;
% check volumes of proposed passes are within range
for g = 1:(length(out_frames_3)-1)
    if isnan(out_frames_3(g+1)) || isnan(out_frames_3(g))
        continue
    else
        p_1 = out_frames_3(g);
        p_2 = out_frames_3(g+1);
        coTemp = surface_recon(poses_downsampled, bound_coords, p_1:p_2, xoffset, yoffset, UStoCam);
        pass_temp = pointCloud(coTemp(:,1:3));
        [k,v_convexhull] = convhull(pass_temp.Location(:,1),pass_temp.Location(:,2),pass_temp.Location(:,3));
        disp(10^6*v_convexhull)
    
        % check if min volume
        if 10^6*v_convexhull < 100
           disp('Volume too small')
    
        elseif 10^6*v_convexhull > 250
           disp('Volume too large')
          
        else
            out_frames_4(ind_4) = out_frames_3(g);
            out_frames_4(ind_4+1) = out_frames_3(g+1);
            out_frames_4(ind_4+2) = NaN;
            ind_4 = ind_4 + 3;
        end
    end
end

out_frames_4 = out_frames_4(1:ind_4-2);
