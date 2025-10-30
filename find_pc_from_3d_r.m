function [var_pc, v_convexhull] = find_pc_from_3d_r(mat_to_use, step_size, start_x, start_y, start_z)

% input voxel grid and resolution of vox grid
% outputs both point cloud and convex hull volume estimate

% first convert to point cloud
xyz_kidney_var = [];
ind = 1;
s = size(mat_to_use);

for a = 1:s(1)
    for b = 1:s(2)
        for c = 1:s(3)
            % calculate voxel coordinates from start point and vox grid resolution
            x_ = start_x + a*step_size;
            y_ = start_y + b*step_size;
            z_ = start_z + c*step_size;

            % if a voxel is non-zero, add a point to point cloud
            if mat_to_use(a,b,c) > 0
                xyz_kidney_var(ind,:) = [x_ y_ z_];
                ind = ind + 1;
            end
        end
    end
end

% flags to allow for no volume calc if point cloud empty or very small
if isempty(xyz_kidney_var)
    v_convexhull = 0;
    var_pc = pointCloud([0 0 0]);
else
    var_pc = pointCloud(xyz_kidney_var);

    % very sprase clouds often cannot calc volume
    if var_pc.Count < 20
        v_convexhull = 0;
    else
        % simple convex hull volume (fast) 
        [k,v_convexhull] = convhull(var_pc.Location(:,1),var_pc.Location(:,2),var_pc.Location(:,3));
    end
end

