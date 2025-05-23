%% Metadata
% SPCV Spring 25 - Final Project
% Name: Nicholas Luis
% PSU ID: NML5604
% I have completed this with integrity

% Goals:
%   - Get silhouette of object
%       > adjust threshold value until no part of the statue is missing 
%         (alternatively can do object tracking and remove everything else)
%   - Create boundary box filled with voxels
%       > determine the min/max values of the coordinates for the box
%         (done via trial and error) 
%       > define the number of voxels in each direction
%   - Determine if a voxel contains our target object
%       > get the center point of each voxel
%       > Transform the voxel center points into 2D image coordinates
%       > Check if the pixel at the coordinate is white in the thresholded
%         image. If so, add an occupancy score to that voxel
%       > repeat for each of the 18 images

clear; clc; close all;

%% Setup

% Use these variables to enable/disable different parts of the script.

loadImages           = true;  
displayVolumeCorners = true;
computeVisualHull    = true;
displayIsoSurface    = true;


%% Task B silhouette threshold

% This should be a suitable value between 0 and 255
silhouetteThreshold = 100;% Enter Value Here

%% Task C define bounding box

% It should be as small as possible, but still contain the whole region of
% interest.
bbox = [-0.25 -2.25 -1.75; 2.25 2.75 1.75]; %[minX minY minZ; maxX maxY maxZ]  
% Creating evenly spaced voxels
volumeX = 200; % Select # voxels in x
volumeY = 600; % Select # voxels in y
volumeZ = round( volumeX*((bbox(2,3)-bbox(1,3)) / (bbox(2,1)-bbox(1,1))) ); % Create evenly distributed voxels in the xz plane
volumeThreshold = 16; 

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
        title(sprintf('Image %d',n));
        subplot(1,2,2);
        imshow(double(rgb2gray(ims{n}))/255.*sils{n});
        title(sprintf('Image %d',n));
        drawnow;
        pause(0);
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
        title(sprintf('Image %d',n));
        drawnow;
        pause(0);

    end
end
%% Task D Visual Hull
if computeVisualHull
    % Define volume. This is used to store the number of observations for each voxel.
    volume = zeros(volumeX,volumeY,volumeZ);
    
    % calculates center point of each voxel in world coords
    xCorners = linspace( bbox(1,1), bbox(2,1), volumeX+1 );
    yCorners = linspace( bbox(1,2), bbox(2,2), volumeY+1 );
    zCorners = linspace( bbox(1,3), bbox(2,3), volumeZ+1 );
    
    xCenters = ( xCorners(1:end-1) + xCorners(2:end) )/2;
    yCenters = ( yCorners(1:end-1) + yCorners(2:end) )/2;
    zCenters = ( zCorners(1:end-1) + zCorners(2:end) )/2;

    % creates a N*4 matrix to store the 3 coordinates & a column of 1's
    centers = NaN(volumeX*volumeY*volumeZ, 4); % preallocate for speed
    pixel2voxel = NaN(volumeX*volumeY*volumeZ, 3); % Maps the coordinates of each centerpoint to a corresponding voxel
    indx = 1;
    for i = 1:volumeX
        for j = 1:volumeY
            for k = 1:volumeZ
                centers(indx, :) = [xCenters(i), yCenters(j), zCenters(k), 1];
                pixel2voxel(indx, :) = [i, j, k];
                indx = indx + 1;
            end
        end
    end

    % Plots the centers of the voxels on the images
    for n = 1:numCameras
        
        % Converts the center points of each voxel from volume coords to pixel coords
        pcenters = Ps{n}*centers'; % Convert to image coordinates
        pcenters = pcenters./repmat(pcenters(3,:),3,1); % Scaling so that third element is 1
        pcenters = round(pcenters); % Round to the nearest pixel

        % Ensures that the transformed coordinates are within the bounds of the image
        [height, width] = size(sils{n});
        for i = 1 : length(pcenters)
            % Checking the u values (columns)
            if ((pcenters(1,i) <= 0) || (pcenters(1,i) > width))
                pcenters(1,i) = NaN;
            end

            % Checking the v values (rows)
            if ((pcenters(2,i) <= 0) || (pcenters(2,i) > height))
                pcenters(2,i) = NaN;
            end
        end

        figure(3); 
        imshow(sils{n}); 
        hold on;
        plot(pcenters(1,:), pcenters(2,:), 'g.');
        title(sprintf('Image %d',n));
        hold off;
        drawnow;
        pause(0);

        % Checks if the voxel center lands on a pixel that is a white or black
        voxelCounter = zeros(length(pcenters), 1);
        for i = 1 : length(pcenters)
            image = sils{n};
            u = pcenters(1,i); 
            v = pcenters(2,i);
            
            if isfinite(u) && isfinite(v)
               
                % if the voxel center is placed on a white pixel, add it to the counter
                if (image(v,u) == 1)
                    xIndex = pixel2voxel(i,1);
                    yIndex = pixel2voxel(i,2);
                    zIndex = pixel2voxel(i,3);
    
                    volume(xIndex, yIndex, zIndex) = 1 + volume(xIndex, yIndex, zIndex);
                end
            
            end

        end % loops through each of the voxels

    end % loops through each of the 18 images

end % entire visual hull code

%% Display the isosurface

% The volume threshold is used to identify the voxels that have no
% intersection at all. For n images, the volume threshold has to be a
% number less or equal to than n. Discuss in your report the significance
% of this number.


if displayIsoSurface
    % display result
    figure(4);
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


