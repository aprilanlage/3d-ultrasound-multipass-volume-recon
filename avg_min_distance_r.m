function avg_dist = avg_min_distance_r(pc1, pc2)

% to calculate average minimum distance (AMD) shape metric 
% input two point clouds to measure shape similarity between them

% density will down sample sampled points
density = 5;

num_pts = pc1.Count;
dis_sum = 0;

% find nearest point on point cloud 2 from every point on point cloud 1 and
% sum the distances
for p_1 = 1:density:num_pts
    [ind, dis] = findNearestNeighbors(pc2,pc1.Location(p_1,:),1);
    dis_sum = dis_sum + dis;
end

num_pts_2 = pc2.Count;
dis_sum_2 = 0;

% find nearest point on point cloud 1 from every point on point cloud 2 and
% sum the distances
for p_2 = 1:density:num_pts_2
    [ind, dis] = findNearestNeighbors(pc1,pc2.Location(p_2,:),1);
    dis_sum_2 = dis_sum_2 + dis;
end

% find averages and account for downsampling 
a_1 = dis_sum/(num_pts/density);
a_2 = dis_sum_2/(num_pts_2/density);

% AMD is mean of two minimum average distances 
avg_dist = mean([a_1 a_2]);