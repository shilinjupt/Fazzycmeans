%-------------------------------------------------------------------------
% FUZZY C-MEANS - Image Segmentation
%
% Author: Moises Mendes de Assis
% Course: Fuzzy Systems
% Date: April/2018
%
% Description: Implementation of Fuzzy C-Means clustering method, which
%              uses fuzzy logic and fuzzy sets concepts. The method is
%              applied to image segmentation.
%-------------------------------------------------------------------------
clc;
close all;
clear all;

%---------------------------------------------------------------
% Load and plot data

fnumber = '003'; %number of image in folter ImagensTeste
filename = ['ImagensTeste/photo', fnumber, '.jpg'];

img = imread(filename);

red = reshape(img(:, :, 1), [], 1);
green = reshape(img(:, :, 2), [], 1);
blue = reshape(img(:, :, 3), [], 1);

red = double(red);
green = double(green);
blue = double(blue);

figure;
plot3(red, green, blue,'b.')
grid on;

figure
imshow(img)

largura = size(img, 1);
altura = size(img, 2);

X = [red, green, blue];

%------------------------------------------------------------------------
% Fuzzy C-Means algorithm

C = 4; % define number of clusters
m = 2; % exponent aplied to the membership matrix for centroids calculation
dim = size(X,2); % number of dimensions for single point
n = size(X,1); % number of points (samples)


membershipMatrix = zeros(n,C);
clusterLabels = zeros(n,1);
clusterCenters = zeros(C,dim); % centroids vector [each line is a centroid]
oldClusters = clusterCenters;

% STEP 1: Initialize Membership Matrix
for i=1:n
    membershipMatrix(i,:) = randi(C, 1, C);
    membershipMatrix(i,:) = membershipMatrix(i,:)./sum(membershipMatrix(i,:)); 
    % above row makes sum of each row equals 1
end

changes = true;
iter = 0;

while (changes)
    
    % STEP 2: Calculate Cluster Center
    for j=1:C
        u = membershipMatrix(:, j); %vector
        uraised = u.^m; %vector
        denominator = sum(uraised); %scalar
        
        numerator = zeros(size(X(1,:))); %vector
        for i=1:n
            x = double(X(i,:));
            produto = uraised(i).*x;
            numerator = numerator + produto;
        end
        clusterCenters(j, :) = numerator/denominator; %vector
    end
        
    % STEP 3: Update Membership Value
    p = 2/(m-1); % exponent
    distance = zeros(1,C);
    
    for i=1:n
        x = double(X(i,:)); %get a point in the dataset
        
        % distance between xi and each centroid
        for j=1:C
            distance(j) = norm(clusterCenters(j, :) - x);
        end
        distance = distance.^2;
        
        for j=1:C
            lower = distance(j)/sum(distance);
            lower = lower^p;
            membershipMatrix(i,j) = 1/lower;
        end
        
        
    end
        
    % sum of the row equals 1
    for i=1:n
        membershipMatrix(i,:) = membershipMatrix(i,:)./sum(membershipMatrix(i,:));
    end
    
    % STEP 4: Get Cluster Labels
    for i=1:n
        [~, id] = max(membershipMatrix(i, :));
        clusterLabels(i) = id;
    end
    
    % stop criteria: centroids do not change anymore
    epsilon = 0.1;
    if (max(max(clusterCenters - oldClusters)) < epsilon)
        changes = false;
    else
        oldClusters = clusterCenters;
    end
    
    iter = iter + 1;
    disp(['iter = ' num2str(iter)])
end

% Assign each pixel the color of his group
newimg = X;
for i=1:length(clusterLabels)
    newimg(i, :) = clusterCenters(clusterLabels(i),:);
    
end

% Format and plot new image
r = reshape(newimg(:, 1), largura, altura);
g = reshape(newimg(:, 2), largura, altura);
b = reshape(newimg(:, 3), largura, altura);

finalimg = uint8([r,g,b]);
finalimg = reshape(finalimg,largura, altura, 3);

figure;
imshow(finalimg)

