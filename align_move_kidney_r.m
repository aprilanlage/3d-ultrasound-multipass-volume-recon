function new_pc = align_move_kidney_r(kidney_pt_cloud, first_axis, second_axis, center_x, center_y, center_z)

% takes in point cloud
% aligns first pca component to first specified axis and second component to second specified axis
% zeros center of kidney point cloud

kidney_points = kidney_pt_cloud.Location;

%pcshow(kidney_pt_cloud)

% first component
% calculate existing long axis
coeffs = pca(kidney_points);
pcs = kidney_points * coeffs;

pc1Min = min(pcs(:, 1));
pc1Max = max(pcs(:, 1));
pc1Mean = mean(pcs(:, 1));

pc1Limits = [pc1Min - pc1Mean, 0, 0; pc1Max - pc1Mean, 0, 0];
xyz1Limit = pc1Limits / coeffs;

first_comp = [xyz1Limit(1, 1) - xyz1Limit(2, 1), xyz1Limit(1, 2) - xyz1Limit(2, 2), xyz1Limit(1, 3) - xyz1Limit(2, 3)];
%x_axis = [1 0 0];

theta_x = acos(dot(first_comp, first_axis)/(norm(first_comp)*norm(first_axis)));
axis_x = cross(first_comp, first_axis)/norm(cross(first_comp, first_axis));

% Construct the rotation matrix using Rodrigues' formula
K_x = [0 -axis_x(3) axis_x(2); axis_x(3) 0 -axis_x(1); -axis_x(2) axis_x(1) 0];
R_x = eye(3) + sin(theta_x)*K_x + (1-cos(theta_x))*K_x*K_x;

t_1 = rigidtform3d(R_x, [0 0 0]);

% transform input point cloud longest axis 
p_1 = pctransform(kidney_pt_cloud,t_1);
%pcshow(p)

% add second componenet
coeffs_2 = pca(p_1.Location);
pcs_2 = p_1.Location * coeffs_2;

pc2Min = min(pcs_2(:, 2));
pc2Max = max(pcs_2(:, 2));
pc2Mean = mean(pcs_2(:, 2));

pc2Limits = [0, pc2Min - pc2Mean, 0; 0, pc2Max - pc2Mean, 0];
xyz2Limit = pc2Limits / coeffs_2;

second_comp = [xyz2Limit(1, 1) - xyz2Limit(2, 1), xyz2Limit(1, 2) - xyz2Limit(2, 2), xyz2Limit(1, 3) - xyz2Limit(2, 3)];
%y_axis = [0 1 0];

theta_y = acos(dot(second_comp, second_axis)/(norm(second_comp)*norm(second_axis)));
axis_y = cross(second_comp, second_axis)/norm(cross(second_comp, second_axis));

K_y = [0 -axis_y(3) axis_y(2); axis_y(3) 0 -axis_y(1); -axis_y(2) axis_y(1) 0];
R_y = eye(3) + sin(theta_y)*K_y + (1-cos(theta_y))*K_y*K_y;

t_2 = rigidtform3d(R_y, [0 0 0]);

% transform previous point cloud to second axis
p_2 = pctransform(p_1,t_2);

% place center of kidney at zero
x_kidney = mean([p_2.XLimits(1) p_2.XLimits(2)]);
y_kidney = mean([p_2.YLimits(1) p_2.YLimits(2)]);
z_kidney = mean([p_2.ZLimits(1) p_2.ZLimits(2)]);

% move kidney 
delta_x = center_x - x_kidney;
delta_y = center_y - y_kidney;
delta_z = center_z - z_kidney;

% translate kidney
tform = rigidtform3d([0 0 0],[delta_x delta_y delta_z]);
new_pc = pctransform(p_2,tform);
