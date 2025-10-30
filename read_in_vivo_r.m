function [frames, poses_downsampled, bound_coords, m_pix, xoffset, yoffset, UStoCam] = read_in_vivo_r(folder, bag_name, us_name, seg_folder, probe_type, b_mode)

% in vivo scans only
% this function takes in file paths and probe type
% unwraps data and readies it for processing

% DICOM ultrasound file
info = dicominfo(strcat(folder,us_name));

% read info
frames = info.NumberOfFrames;
fr_images = info.CineRate;
US_length = info.EffectiveDuration;

% Offset to the center of the image or top of ultrasound image depending on what you measured to
xoffset = double(info.SequenceOfUltrasoundRegions.Item_1.ReferencePixelX0);
yoffset = double(info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0) - double(info.SequenceOfUltrasoundRegions.Item_1.ReferencePixelY0);

% BAG file - from the camera, assuming length is longer than US images
filepath = strcat(folder,bag_name);
reader = rosbagreader(filepath);
bag = rosbag(filepath);
fr_poses = bag.NumMessages/(bag.EndTime-bag.StartTime);
bagselect = select(reader,"Time",[(bag.EndTime - US_length) (bag.EndTime)]);
bSel = select(bagselect, 'Topic', '/odometry/filtered', 'MessageType', 'nav_msgs/Odometry');
msg = readMessages(bSel);

% Poses from BAG file
getRot = @(m) [m.Pose.Pose.Orientation.W, m.Pose.Pose.Orientation.X, m.Pose.Pose.Orientation.Y,...
    m.Pose.Pose.Orientation.Z];
getTrans = @(m)  [m.Pose.Pose.Position.X, m.Pose.Pose.Position.Y, m.Pose.Pose.Position.Z];

Q = cellfun(getRot, msg, 'UniformOutput', false);
T = cellfun(getTrans, msg, 'UniformOutput', false);

Q = cell2mat(Q);
T = cell2mat(T);

% Cycle through all the poses
n = size(T, 1); % number of frames or data points
poses = zeros(3, 4, n);

for jj = 1:n
    translation = T(jj, :);
    rotation = quat2rotm(Q(jj, :));
    poses(:, :, jj) = [rotation, translation'];
end

times = round(1:fr_poses/fr_images:n);

% downsample poses
poses_downsampled = zeros(3,4,frames);
for slice = 1:(frames-120)
    poses_downsampled(:,:,slice) = poses(:, :, times(slice));
end

% Transformation between center of camera to ultrasound probe
% units = cm/pixel, then convert to m/pixel
m_pix = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX/100.0;

% rigid transformation between camera and ultrasound iamging plane
% depends on probe
% assumes 3D printed attachment is used
% re-calibration needed infrequently
switch probe_type
    case 'GE_LOGIQE9_curvilinearProbe'
        UStoCam = [0, 0, -39.8/1000;...
            -m_pix, 0, -17.5/1000;...
            0, -m_pix, -106.5/1000;...
            0, 0, 1];
    case 'GE_LOGIQE9_linearProbe'
        UStoCam = [0, 0, -50.8/1000;...
            -m_pix, 0, -17.5/1000;...
            0, -m_pix, -90/1000;...
            0, 0, 1];
end


% load file names
fds = fileDatastore(strcat(seg_folder,"\*.nii.gz"), 'ReadFcn', @importdata);
labelNames = fds.Files;

if b_mode == 0
    % zip images into one file - bound_coords is output
    % assuming every 2 images are segmented, starting from 1
    bound_coords = cell(frames,2);
    
    for jj = 1:(frames-10)
        if mod(jj,2) == 1
            ind = 1 + round(jj/2);
            label = niftiread(labelNames{ind});
    
            if length(size(label)) == 3
                label = label(:,:,1);
            end
    
            % this area limit is in pixels, can be changed
            label_f = bwpropfilt(logical(label),'Area',[1500 Inf]);
    
            [OUT] = edge2(label_f);
            [v,u] = find(OUT); 
            bound_coords(jj,:) = {v,u};
        else
            continue
        end
    end
elseif b_mode == 1
    % zip images into one file - bound_coords is output
    % assuming every 2 images are segmented, starting from 1
    bound_coords = cell(frames,3);
    
    for jj = 1:(frames-10)
        if mod(jj,2) == 1
            ind = 1 + round(jj/2);
            label = niftiread(labelNames{ind});
    
            if length(size(label)) == 3
                label = label(:,:,1);
            end
    
            % this area limit is in pixels, can be changed
            label_f = bwpropfilt(logical(label),'Area',[15000 Inf]);
    
            [OUT] = edge2(label_f);
            [v,u] = find(OUT); 
            bound_coords(jj,1:2) = {v,u};

            % look up pixel intensity
            img = dicomread(strcat(folder,us_name),'frames',jj);
            img = double(img(:,:,1));
            is = zeros(length(u),1);
            for h = 1:length(u)
                is(h) = img(u(h),v(h));
            end
            bound_coords(jj,3) = {is};
        else
            continue
        end
    end
end
end
