function [new_voxels] = shape_variance_unnorm_r(old_voxel_grid, cloud_to_add, step_size, start_x,start_y,start_z)
% for adding point cloud points to a voxel grid 

vox_var = old_voxel_grid;

% downsample to same density
pc_temp_down = pcdownsample(cloud_to_add,'gridAverage',step_size);

for pt = 1:pc_temp_down.Count
    point = pc_temp_down.Location(pt,:);
    rounded_pt = 2*round((point/2),3); % will need to change this for other step sizes
    x_index = int16((rounded_pt(1) - start_x)/step_size);
    y_index = int16((rounded_pt(2) - start_y)/step_size);
    z_index = int16((rounded_pt(3) - start_z)/step_size);

    % check if there are specified intensities 
    if nnz(pc_temp_down.Intensity) > 0
        vox_var(x_index,y_index,z_index) = vox_var(x_index,y_index,z_index) + pc_temp_down.Intensity(pt);
    else
        vox_var(x_index,y_index,z_index) = vox_var(x_index,y_index,z_index) + 15;
    end

    % for including neighborhood blur
    for x_wiggle = (rounded_pt(1)-step_size):step_size:(rounded_pt(1)+step_size)
        for y_wiggle = (rounded_pt(2)-step_size):step_size:(rounded_pt(2)+step_size)
            for z_wiggle = (rounded_pt(3)-step_size):step_size:(rounded_pt(3)+step_size)
                x_index = int16((x_wiggle - start_x)/step_size);
                y_index = int16((y_wiggle - start_y)/step_size);
                z_index = int16((z_wiggle - start_z)/step_size);
                vox_var(x_index,y_index,z_index) = vox_var(x_index,y_index,z_index) + 5;
            end
        end
    end
end


new_voxels = vox_var;
end