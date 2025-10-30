function [frames, times, poses, m_pix, xoffset, yoffset, UStoCam] = read_in_phantom_r(folder, bag_name, us_name, probe_type)

% phantom scans only
% this function takes in file paths and probe type
% unwraps data and readies it for processing

% DICOM ultrasound file
info = dicominfo(strcat(folder,us_name));

frames = info.NumberOfFrames; % total number of frames
fr_images = info.CineRate; % frame rate
US_length = info.EffectiveDuration;

% Offset to the center of the image or top of ultrasound image depending on what you measured to
xoffset = double(info.SequenceOfUltrasoundRegions.Item_1.ReferencePixelX0);
yoffset = double(info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0) - double(info.SequenceOfUltrasoundRegions.Item_1.ReferencePixelY0);

% BAG file - from the camera, assuming length is longer than US images
filepath = strcat(folder,bag_name);
reader = rosbagreader(filepath);
bag = rosbag(filepath);
fr_poses = bag.NumMessages/(bag.EndTime-bag.StartTime);
bagselect = select(reader,"Time",[(bag.EndTime - US_length) bag.EndTime]);
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

% down sampling between poses and images
times = round(1:fr_poses/fr_images:n);

end
