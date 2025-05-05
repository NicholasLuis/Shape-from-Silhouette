%% Metadata
% SPCV Spring 24 - Project 1
% Name: Nicholas Luis
% PSU ID: NML5604 (930841391)

% Goals:
%   - Get silhouette of object
%       > adjust threshold value until it is good (alternatively, can do
%         object tracking of the statue in between images)
%   - Determine if a voxel contains our target object
%       > get the center point of each voxel
%       > Get the voxel center points in the 2D image coordinates
%       > Check if the pixel at the coordinate is white in the thresholded
%         image. If so, add an occupancy score to that voxel
%       > repeat

%       > object tracking of the statue between images? (alternatively,
%         find a silouette threshold value that only extracts the statue)


clear; clc; close all;

%% Setup

% Use these variables to enable/disable different parts of the script.

loadImages           = true;  
displayVolumeCorners = true;
computeVisualHull    = true;
displayIsoSurface    = true;


%% Task B silhouette threshold

% This should be a suitable value between 0 and 255
silhouetteThreshold = 150;% Enter Value Here

%% Task C define bounding box

% It should be as small as possible, but still contain the whole region of
% interest.
bbox = [-1 -1 -2; 4 4 4]; %[minX minY minZ; maxX maxY maxZ]  
volumeX = 10; % Select # voxels in x
volumeY = 10; % Select # voxels in y
volumeZ = 20; % Select # voxels in z
volumeThreshold = 10; % What's a good value here?

% The volume threshold is used to identify the voxels that have no
% intersection at all. For n (18) images, the volume threshold has to be a
% number less or equal to than n. Discuss in your report the significance
% of this number.

%% Load in images

numCameras = 18;

if loadImages
    % Load silhouette images and projection matrices
    for n=1:numCameras
        % Projection matrices : 3x4 matrix of [R,t]. 
        % This does not normalize the depth value
        Ps{n} = textread(sprintf('dataSFS/david_%02d.pa',n-1));
        Ps{n} = [eye(3,2) [1 1 1]']*Ps{n};  % add 1 for one-based indices
        % Images
        ims{n} = imread(sprintf('dataSFS/david_%02d.jpg',n-1));
        % Silhouettes: pixels set to one when larger brightness than
        % threshhold
        sils{n} = rgb2gray(ims{n})>silhouetteThreshold;
        
        figure(1);
        subplot(1,2,1);
        imshow(sils{n});
        subplot(1,2,2);
        imshow(double(rgb2gray(ims{n}))/255.*sils{n});
        drawnow;
    end
end

%% Define transformation from volume to world coordinates.

T = [eye(4,3) [bbox(1,:) 1]'] * ...
    diag([(bbox(2,1)-bbox(1,1))/volumeX ...
          (bbox(2,2)-bbox(1,2))/volumeY ...
          (bbox(2,3)-bbox(1,3))/volumeZ ...
          1]);
T = [1  0 0 0; ...
     0  0 1 0; ...  % flip y and z axes for better display in matlab figure (isosurface)
     0 -1 0 0; ...
     0  0 0 1] * T;
T = T*[eye(4,3) [-[1 1 1] 1]'];  % subtract 1 for one-based indices


if displayVolumeCorners
    % Draw projection of volume corners.
    for n=1:numCameras
        figure(2);
        hold off;
        imshow(ims{n});
        hold on;
        corners = [[      0       0       0 1]' ...
                   [      0       0 volumeZ 1]' ...
                   [      0 volumeY       0 1]' ...
                   [      0 volumeY volumeZ 1]' ...
                   [volumeX       0       0 1]' ...
                   [volumeX       0 volumeZ 1]' ...
                   [volumeX volumeY       0 1]' ...
                   [volumeX volumeY volumeZ 1]'];
        pcorners = Ps{n}*T*corners; % Convert to image coordinates
        pcorners = pcorners./repmat(pcorners(3,:),3,1); % Scaling so that third element is 1
        plot(pcorners(1,:),pcorners(2,:),'g*');
        drawnow;
        pause(0.1);

    end
end
%% Task D Visual Hull

if computeVisualHull
    % Define volume. This is used to store the number of observations for each voxel.
    volume = zeros(volumeX,volumeY,volumeZ);
    
    % Visual hull computation    
    %   - For each image add one to the voxel if projection is within silhouette region.
    %   - Be careful with the order of coordinates. The point is stored as
    %     (x,y,z), but matrix element access in Matlab is mat(row,col).

    % -----------------  ENTER YOUR CODE HERE  ----------------- %
    
    % calculates center point of each voxel in volume coordinates

        % Gets the corners of the voxels 
        xCorners = linspace(bbox(1,1),bbox(2,1),volumeX+1);
        yCorners = linspace(bbox(1,2),bbox(2,2),volumeY+1);
        zCorners = linspace(bbox(1,3),bbox(2,3),volumeZ+1);

        % gets the center points of each voxel
        xCenters = NaN(length(xCorners)-1,1); % pre-allocating for speed
        yCenters = NaN(length(yCorners)-1,1); 
        zCenters = NaN(length(zCorners)-1,1);
        for indx = 1 : volumeX
            xCenters(indx) = mean(xCorners(indx):xCorners(indx+1));
        end
        for indx = 1: volumeY
            yCenters(indx) = mean(yCorners(indx):yCorners(indx+1));
        end
        for indx = 1: volumeZ
            zCenters(indx) = mean(zCorners(indx):zCorners(indx+1));
        end

        % creates a matrix of all the center points. Each cell corresponds
        % to which voxel it is. And inside each cell is a 3x1 vector that
        % describe its volume coordinates
        centers = cell(volumeX, volumeY, volumeZ); % preallocate for speed
        for i = 1:volumeX
            for j = 1:volumeY
                for k = 1:volumeZ
                    centers{i,j,k} = cell(3,1);   % each d{i,j,k} is a 3×1 cell
                end
            end
        end

    % Converts the center points of each voxel from volume coords to world coordinates 

    % Converts the center points of each voxel from world coords to pixel coords

    % Checks if the pixel is white at the given voxel pixel coordinates
    for i = 1 : volumeX
        for j = 1 : volumeY
            for k = 1: volumeZ 


                % 
                
            end
        end
    end
    

    %---------------------------------------------------------------

end

%% Display the isosurface

% The volume threshold is used to identify the voxels that have no
% intersection at all. For n images, the volume threshold has to be a
% number less or equal to than n. Discuss in your report the significance
% of this number.


if displayIsoSurface
    % display result
    figure(3);
    clf;
    grid on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    hold on;
    [xMesh, yMesh, zMesh] = meshgrid(1:volumeY,1:volumeX,1:volumeZ);
    pt = patch(isosurface(yMesh, xMesh, zMesh, volume, volumeThreshold));
    set(pt,'FaceColor','red','EdgeColor','none');
    axis equal;
    daspect([volumeX/(bbox(2,1)-bbox(1,1)) volumeY/(bbox(2,2)-bbox(1,2)) volumeZ/(bbox(2,3)-bbox(1,3))]);
    camlight(0,0);
    camlight(180,0);
    camlight(0,90);
    camlight(0,-90);
    lighting phong;
    view(30,30);
end


