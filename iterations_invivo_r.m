%% iterations
% in vivo

% read in data
% poses
bag_name = "C:\Users\april\Desktop\In vivo kidney study 2025\PartID_34\Camera data\camPoses_partID_34_2025-03-20-12-39-05.bag";

% US
us_file = "C:\Users\april\Desktop\In vivo kidney study 2025\PartID_34\US_scans\P3KDJDGE";

b = load('C:\Users\april\Desktop\PhD research\CODE\bc_ID34_scan2.mat');
bound_coords = b.bound_coords;
[frames, times, poses, m_pix, xoffset, yoffset, UStoCam] = read_in_r('', bag_name, us_file, 'GE_LOGIQE9_curvilinearProbe');

% downsample poses
poses_downsampled = zeros(3,4,frames);
for slice = 1:(frames-120)
    poses_downsampled(:,:,slice) = poses(:, :, times(slice));
end

%passes = pass_detect_r(poses_downsampled,bound_coords, xoffset, yoffset, UStoCam);
passes = [65 140 NaN 175 260 NaN 305 430 NaN 470 560 NaN 650 800];

plot_skeleton_outline_r(poses_downsampled, bound_coords, passes(1):passes(2), xoffset, yoffset, UStoCam)

%% pass alignment
% align centers of segmentations
poses_centered = poses_downsampled;

% get centers
for i = 1:(length(passes)-1)
    if isnan(passes(i+1)) || isnan(passes(i))
        continue
    else
        %disp([passes(i) passes(i+1)])
        coTemp = surface_recon_r(poses_downsampled, bound_coords, (passes(i)+5):(passes(i+1)-5), xoffset, yoffset, UStoCam);
        cen = [mean(coTemp(:, 1)), mean(coTemp(:, 2)), mean(coTemp(:, 3))];
        % put everybody at (0,0,0)
        poses_centered(1:3,4,passes(i):(passes(i+1))) = poses_downsampled(1:3,4,passes(i):(passes(i+1))) - [cen(1) cen(2) cen(3)]';
    end
end

%plot_skeleton_outline(poses_centered, bound_coords, (passes(61)+0):(passes(62)-0), xoffset, yoffset, UStoCam)

%% method for volume
num_passes = length(passes);
vols = zeros(num_passes,1);
step_size = 0.002;

% iterations
it = 100;
vol_avgs = zeros(1,it);

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
                coTemp = surface_recon_r(poses_centered_icp, bound_coords, (passes(i)+1):(passes(i+1)-1), xoffset, yoffset, UStoCam);
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
    x_extent = -0.08:step_size:0.1;
    y_extent = -0.12:step_size:0.1;
    z_extent = -0.11:step_size:0.15;
    old_voxel_grid = zeros(length(x_extent),length(y_extent),length(z_extent));

    coTemp = surface_recon_r(poses, bound_coords, passes(1):passes(2), xoffset, yoffset, UStoCam);
    pass_temp = pointCloud(coTemp(:,1:3));
    [new_voxels] = shape_variance_unnorm_r(old_voxel_grid, pass_temp, step_size, x_extent(1), y_extent(1), z_extent(1));

    for pass_num = 2:(num_passes-1)
        % flag to check for breathing spots
        if isnan(passes(pass_num+1)) || isnan(passes(pass_num))
            continue
        else
            % add another pass to the variance grid
            coTemp = surface_recon_r(poses, bound_coords, (passes(pass_num)+1):(passes(pass_num+1)-1), xoffset, yoffset, UStoCam);
            pass_temp = pointCloud(coTemp(:,1:3),Intensity=coTemp(:,5));
            [new_voxels] = shape_variance_unnorm_r(new_voxels, pass_temp, step_size, x_extent(1), y_extent(1), z_extent(1));

            %threshold
            new_voxels_thres = zeros(size(new_voxels));
            s = size(new_voxels);
            peaks = zeros(s);

            thres_factor = 0.7;

            % threshold per slice to avoid skeleton volume
            for select = 1:s(1)
                im = squeeze(new_voxels(select,:,:));
                im_new = im;

                n = nnz(im);

                if n > 0
                    im_long = reshape(im,[],1);
                    im_sort = sort(im_long,"descend");
                    threshold = im_sort(round(thres_factor*n));

                    k = im<threshold;
                    im_new(k) = 0;

                    new_voxels_thres(select,:,:) = im_new;
                else
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
                            peaks(i,j,l) = 1;
                        end
                    end
                end
            end

            %volumeViewer(peaks)

            % calculate average ring thickness
            % if too thick, erode again
            th = calc_thickness_r(peaks);
            if th > 3
                %disp(th)
                erodedvol_2 = imerode(erodedvol, se_vol);
                peaks_new = zeros(s);
                for i = 1:s(1)
                    for j = 1:s(2)
                        for l = 1:s(3)
                            if erodedvol_2(i,j,l) > 0
                                peaks_new(i,j,l) = 1;
                            end
                        end
                    end
                end
                th = calc_thickness_r(peaks_new);
                if th < 2
                    % use previous peaks, with one erode
                    peaks = peaks;
                elseif th > 5
                    % erode again if still thick
                    %disp(th)
                    erodedvol_3 = imerode(erodedvol_2, se_vol);
                    peaks_new_new = zeros(s);
                    for i = 1:s(1)
                        for j = 1:s(2)
                            for l = 1:s(3)
                                if erodedvol_3(i,j,l) > 0
                                    peaks_new_new(i,j,l) = 1;
                                end
                            end
                        end
                    end
                    peaks = peaks_new_new;
                else
                    % use new peaks, with two erode
                    peaks = peaks_new;
                end
            end


            % find volume
            [pc, v_convexhull] = find_pc_from_3d_r(peaks, step_size, x_extent(1), y_extent(1), z_extent(1));
            pcwrite(pc,'measure_in_vivo.ply')
            vol_mesh = pyrunfile("vol_measurement.py 'measure_in_vivo.ply'","v");
            vols(pass_num) = 10^6*abs(double(vol_mesh));

        end
    end
    y = find(vols);
    vol_avgs(c) = mean(vols(y));
    disp([vol_avgs(c) std(vols(y))])

    % stopping criteria
    if (c>2) && (abs(vol_avgs(c) - vol_avgs(c - 1)) < 10) && (abs(vol_avgs(c - 1) - vol_avgs(c - 2)) < 10)
        disp(c)
        break
    else
        continue
    end
end