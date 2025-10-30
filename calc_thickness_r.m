function avg_thick = calc_thickness_r(vox_grid)

% this function calculates the average thickness of a closed ring
% takes in a voxel grid, outputs average thickness in pixels 

total = 0;
sum = 0;
s = size(vox_grid);

% first dimension of voxel grid
for l = 1:s(1)
    % make binary image
    im = logical(squeeze(vox_grid(l,:,:)));
    stats = regionprops(im);

    % if a ring is detected
    if(length(stats) == 1)
        %disp('Ring detected correctly.')
        center = round(stats.Centroid);

        % finds 8 points to measure thickness
        [points] = be_intersection_points_r(im, [center(2), center(1)], 8);

        d = 0;
        for y = 1:8
            if isempty(points{y})
                continue
            else
                % if points exist, find distance between them
                d = d + pdist2(points{y}(1,:),points{y}(2,:));
            end
        end

        sum = sum + d;

        total = total + 1;
    else
        total = total;
    end
end

% second dimension of voxel grid
for m = 1:s(2)
    im = logical(squeeze(vox_grid(:,m,:)));
    stats = regionprops(im);

    if(length(stats) == 1)
        %disp('Ring detected correctly.')
        center = round(stats.Centroid);
        [points] = be_intersection_points_r(im, [center(2), center(1)], 8);

        d = 0;
        for y = 1:8
            if isempty(points{y})
                continue
            else
                d = d + pdist2(points{y}(1,:),points{y}(2,:));
            end
        end

        sum = sum + d;

        total = total + 1;
    else
        total = total;
    end
end

% third dimension of voxel grid 
for n = 1:s(3)
    im = logical(squeeze(vox_grid(:,:,n)));
    stats = regionprops(im);

    if(length(stats) == 1)
        %disp('Ring detected correctly.')
        center = round(stats.Centroid);
        [points] = be_intersection_points_r(im, [center(2), center(1)], 8);

        d = 0;
        for y = 1:8
            if isempty(points{y})
                continue
            else
                d = d + pdist2(points{y}(1,:),points{y}(2,:));
            end
        end

        sum = sum + d;

        total = total + 1;
    else
        total = total;
    end
end

% average thickness over 8 points on detected rings 
avg_thick = sum/total/8;