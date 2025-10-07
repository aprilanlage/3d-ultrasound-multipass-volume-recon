function hausdorff_dist = hausdorff_dist(pc1, pc2)

density = 1;

num_pts = pc1.Count;
max_1 = 0;

for p = 1:density:num_pts
    [ind, dis] = findNearestNeighbors(pc2,pc1.Location(p,:),1);
    if dis > max_1
        max_1 = dis;
    end
end

num_pts_2 = pc2.Count;
max_2 = 0;

for p = 1:density:num_pts_2
    [ind, dis] = findNearestNeighbors(pc1,pc2.Location(p,:),1);
    if dis > max_2
        max_2 = dis;
    end
end

hausdorff_dist = max(max_1, max_2);