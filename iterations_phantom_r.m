%% Main code - iterations of algorithm
% phantom
% multiple sections need to be run

% read in data
% poses
% CHANGE THIS FILE PATH
bag_name = "Phantom_data\camPoses_phantom_scan1";

% US
% CHANGE THIS FILE PATH
us_file = "Phantom_data\US_phantom_scan1";

% segmentations
% CHANGE THIS FILE PATH
b = load('Phantom_data\bc_phantom_scan1_updated.mat');
bound_coords = b.bound_coords;

[frames, times, poses, m_pix, xoffset, yoffset, UStoCam] = read_in_phantom_r('', bag_name, us_file, 'GE_LOGIQE9_curvilinearProbe');

% downsample poses
poses_downsampled = zeros(3,4,frames);
for slice = 1:(frames-120)
    poses_downsampled(:,:,slice) = poses(:, :, times(slice));
end

% pass detection
passes = pass_detect_phantom_r(poses_downsampled, bound_coords, xoffset, yoffset, UStoCam);

% optional plotting to visualize segmentations
%plot_skeleton_outline(poses_downsampled, bound_coords, passes(1):10:passes(8), xoffset, yoffset, UStoCam)

%% rough initial pass alignment - phantom scans 1, 3, 4, 5
% initial align centers of segmentations
poses_centered = poses_downsampled;

% get centers
for i = 1:(length(passes)-1)
    if isnan(passes(i+1)) || isnan(passes(i))
        continue
    else
        %disp([passes(i) passes(i+1)])
        coTemp = surface_recon_r(poses_downsampled, bound_coords, (passes(i)):(passes(i+1)), xoffset, yoffset, UStoCam);
        cen = [mean(coTemp(:, 1)), mean(coTemp(:, 2)), mean(coTemp(:, 3))];
        % put everybody at (0,0,0)
        poses_centered(1:3,4,passes(i):(passes(i+1))) = poses_downsampled(1:3,4,passes(i):(passes(i+1))) - [cen(1) cen(2) cen(3)]';
    end
end

% optional visualization
%plot_skeleton_outline(poses_centered, bound_coords, (passes(1)):(passes(2)), xoffset, yoffset, UStoCam)

%% poses centered based on first pass only - phantom scan 2 only

poses_centered = poses_downsampled;

% get centers
for i = 1:(length(passes)-1)
    if isnan(passes(i+1)) || isnan(passes(i))
        continue
    else
        %disp([passes(i) passes(i+1)])
        % calculate translation for 1st pass only
        if i == 1
            coTemp = surface_recon_r(poses_downsampled, bound_coords, (passes(i)+5):(passes(i+1)-5), xoffset, yoffset, UStoCam);
            cen = [mean(coTemp(:, 1)), mean(coTemp(:, 2)), mean(coTemp(:, 3))];
        end
        % put everybody at (0,0,0)
        poses_centered(1:3,4,passes(i):(passes(i+1))) = poses_downsampled(1:3,4,passes(i):(passes(i+1))) - [cen(1) cen(2) cen(3)]';
    end
end

% optional visualization
%plot_skeleton_outline(poses_centered, bound_coords, (passes(1)):(passes(2)), xoffset, yoffset, UStoCam)

%% align shape model for shape metric evaluation later
poses = poses_centered;
kid_model = pcread('scaled_mean_pc.ply');
kidney_points = kid_model.Location;
kidney_points_small = kidney_points/1000;
kid_model_small = pointCloud(kidney_points_small);

% using first few passes as model
coTemp = surface_recon_r(poses, bound_coords, passes(1):passes(2), xoffset, yoffset, UStoCam);
pass_temp = pointCloud(coTemp(:,1:3));
%pcshowpair(kid_model_small,pass_temp)

% move shape model to align with kidney
% for phantom scans 1, 2, and 4
aligned_kidney = align_move_kidney_r(kid_model_small, [1 0 0], [0 0 -1], 0,0,0);

% for phantom scan 3
%aligned_kidney = align_move_kidney_r(kid_model_small, [0 0 1], [-1 0 0], 0,0,0);

% for phantom scan 5
%aligned_kidney = align_move_kidney_r(kid_model_small, [0 0 1], [0 -1 0], 0,0,0);

best_t = pcregistericp(aligned_kidney,pass_temp);
kidney_reg = pctransform(aligned_kidney,best_t);

% visual check
pcshowpair(kidney_reg,pass_temp)

%% method for volume
num_passes = length(passes);
vols = zeros(num_passes,1);
min_dist = zeros(num_passes,1);
% step size of voxel grid can be changed, 2mm works best with provided data
step_size = 0.002;

% iterate
% this is an upper limit for number of iterations (most scans converge much
% faster)
it = 100;
vol_avgs = zeros(1,it);
dist_avgs = zeros(1,it);

for c = 1:it
    
    if c == 1
        poses = poses_centered;
    else
        % use previous calculated point cloud to register passes again
        poses_centered_icp = poses;

        % re-center pc that we register to
        cen = [mean(pc.Location(:,1)), mean(pc.Location(:,2)), mean(pc.Location(:,3))];
        pc_cen = pctransform(pc,-cen.*ones(length(pc.Location),3));

        % get centers
        for i = 1:(length(passes)-1)
            if isnan(passes(i+1)) || isnan(passes(i))
                continue
            else
                coTemp = surface_recon_r(poses_centered_icp, bound_coords, (passes(i)+5):(passes(i+1)-5), xoffset, yoffset, UStoCam);
                % use center of pass segmentations to register

                pass_temp = pointCloud(coTemp(:,1:3));

                % pc is calculated point cloud from previous round
                best_t = pcregistericp(pass_temp,pc_cen);
                poses_to_change = poses_centered_icp(:,:,passes(i):(passes(i+1)-1));

                for v = 1:(passes(i+1)-passes(i))
                    poses_to_change(4,:,v) = [0 0 0 1];
                    poses_to_change(:,:,v) = best_t.A*poses_to_change(:,:,v);
                end
                poses_centered_icp(:,:,passes(i):(passes(i+1)-1)) = poses_to_change(1:3,:,:);
            end
        end
        poses = poses_centered_icp;
    end

    % kidney with poses centered via poses
    % these extents are hard coded but work for all provided data
    x_extent = -0.1:step_size:0.12;
    y_extent = -0.08:step_size:0.07;
    z_extent = -0.1:step_size:0.1;
    old_voxel_grid = zeros(length(x_extent),length(y_extent),length(z_extent));

    % add first pass to empty voxel grid
    coTemp = surface_recon_r(poses, bound_coords, passes(1):passes(2), xoffset, yoffset, UStoCam);
    pass_temp = pointCloud(coTemp(:,1:3),Intensity=coTemp(:,5));
    [new_voxels] = shape_variance_unnorm_r(old_voxel_grid, pass_temp, step_size, x_extent(1), y_extent(1), z_extent(1));

    for pass_num = 2:(num_passes-1)

        % flag to check for breathing spots
        if isnan(passes(pass_num+1)) || isnan(passes(pass_num))
            continue
        else
            % add another pass to the variance grid
            coTemp = surface_recon_r(poses, bound_coords, passes(pass_num):(passes(pass_num+1)-1), xoffset, yoffset, UStoCam);
            pass_temp = pointCloud(coTemp(:,1:3),Intensity=coTemp(:,5));
            [new_voxels] = shape_variance_unnorm_r(new_voxels, pass_temp, step_size, x_extent(1), y_extent(1), z_extent(1));

            %threshold
            new_voxels_thres = zeros(size(new_voxels));
            s = size(new_voxels);
            peaks = zeros(s);
            peaks_two_erode = zeros(s);
        
            thres_factor = 0.7;

            % threshold per slice to avoid skeleton volume
            for select = 1:s(1)
                im = squeeze(new_voxels(select,:,:));
                im_new = im;

                n = nnz(im);

                if n > 0
                    % flatten image and order from brightest to darkest
                    % pixel values
                    im_long = reshape(im,[],1);
                    im_sort = sort(im_long,"descend");
                    
                    % find the pixels that correspond to above the threshold
                    threshold = im_sort(round(thres_factor*n));
                    k = im<threshold;
                    % make the rest of the pixels 0
                    im_new(k) = 0;

                    new_voxels_thres(select,:,:) = im_new;
                else
                    % if slice is blank, skip
                    new_voxels_thres(select,:,:) = 0;
                end
            end

            % test erode
            se_vol = strel('sphere',2);
            erodedvol = imerode(new_voxels_thres, se_vol);

             for i = 1:s(1)
                for j = 1:s(2)
                    for l = 1:s(3)
                        if erodedvol(i,j,l) > 0
                            % convert to a binary volume
                            peaks(i,j,l) = 1;
                        end
                    end
                end
            end

            % optional visualization
            %volumeViewer(peaks)

            % calculate average ring thickness
            % if too thick, erode again
            th = calc_thickness_r(peaks);
            se_vol = strel('sphere',1);

            if th > 3
                %disp(th)
                erodedvol_2 = imerode(erodedvol, se_vol);
                peaks_two_erode = zeros(s);
                for i = 1:s(1)
                    for j = 1:s(2)
                        for l = 1:s(3)
                            if erodedvol_2(i,j,l) > 0
                                peaks_two_erode(i,j,l) = 1;
                            end
                        end
                    end
                end
                th = calc_thickness_r(peaks_two_erode);

                if th < 2
                    % if too thin, use previous peaks, with one erode
                    peaks = peaks;

                elseif th > 5
                    % erode again if still thick
                    %disp(th)
                    erodedvol_3 = imerode(erodedvol_2, se_vol);
                    peaks_three_erode = zeros(s);
                    for i = 1:s(1)
                        for j = 1:s(2)
                            for l = 1:s(3)
                                if erodedvol_3(i,j,l) > 0
                                    %if new_voxels_thres(i,j,l) > 0
                                    peaks_three_erode(i,j,l) = 1;
                                end
                            end
                        end
                    end
                    peaks = peaks_three_erode;

                else
                    % use peaks with two erode
                    peaks = peaks_two_erode;
                end
            end
            
            % calculate metrics 
            % convert to point cloud
            [pc, v_convexhull] = find_pc_from_3d_r(peaks, step_size, x_extent(1), y_extent(1), z_extent(1));
            
            % need to re-register shape model kidney every time
            best_t = pcregistericp(kidney_reg,pc);
            kidney_reg = pctransform(kidney_reg,best_t);

            %  find AMD
            d = avg_min_distance_r(pc,kidney_reg);
            min_dist(pass_num) = d;

            % using Python meshlab to calc poisson surface and volume
            pcwrite(pc,'measure_this.ply')
            vol_mesh = pyrunfile("vol_measurement.py 'measure_this.ply'","v");
            vols(pass_num) = 10^6*abs(double(vol_mesh));

        end
    end

    y = find(vols);
    vol_avgs(c) = mean(vols(y));
    z = find(min_dist);
    dist_avgs(c) = mean(min_dist(y));
    disp([vol_avgs(c) dist_avgs(c)])

    % stopping criteria
    % currently set to 3 volumes within 10mL in a row
    if (c>2) && (abs(vol_avgs(c) - vol_avgs(c - 1)) < 10) && (abs(vol_avgs(c - 1) - vol_avgs(c - 2)) < 10)
        % calculate final AMD using Poisson surface 
        [mesh] = pc2surfacemesh(pc,"Poisson");
        pc = mesh2pc(mesh);
        d = avg_min_distance_r(pc,kidney_reg);
        hd = hausdorff_dist_r(pc, kidney_reg);

        % display stats
        disp([c mean(vol_avgs((c-2):c)) 1000*d 1000*hd])
        break
    else
        continue
    end
end
