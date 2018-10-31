%-------------------------------------------------------------------------
% FUZZY C-MEANS
%
% Author: Moises Mendes de Assis
% Course: Fuzzy Systems
% Date: April/2018
%
% Description: Implementation of Fuzzy C-Means clustering method, which
%              uses fuzzy logic and fuzzy sets concepts. The method is 
%              applied to a simple data set with 4 clusters.
%-------------------------------------------------------------------------
clc;
close all;
clear all;

% Load and plot data
load fcm_dataset.mat;
X = x;

figure(1);
hold on;
plot(X(:,1), X(:,2),'b.');
grid on;

%------------------------------------------------------------------------
% Fuzzy C-Means algorithm

C = 4; % define number of clusters
m = 2; % exponent aplied to the membership matrix for centroids calculation
dim = size(X,2); % number of dimensions for single point
n = size(X,1); % number of points (samples)


membershipMatrix = zeros(n,C);
clusterLabels = zeros(n,1);
oldClusters = zeros(n,1);
clusterCenters = zeros(C,dim); % centroids vector [each line is a centroid]

% STEP 1: Initialize Membership Matrix
for i=1:n
    membershipMatrix(i,:) = randi(C, 1, C);
    membershipMatrix(i,:) = membershipMatrix(i,:)./sum(membershipMatrix(i,:)); 
    % above row makes sum of each row equals 1
end

changes = true;
iter = 1;

while (changes)
    
    % STEP 2: Calculate Cluster Center
    for j=1:C
        u = membershipMatrix(:, j); %vector
        uraised = u.^m; %vector
        denominator = sum(uraised); %scalar
        
        numerator = zeros(size(X(1,:))); %vector
        for i=1:n
            x = X(i,:);
            produto = uraised(i).*x;
            numerator = numerator + produto;
        end
        clusterCenters(j, :) = numerator/denominator; %vector
    end
    
    % STEP 3: Update Membership Value
    p = 2/(m-1); % exponent
    distance = zeros(1,C);
    
    for i=1:n
        x = X(i,:); %get a point in the dataset
        
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
    
    % ploting the new centroids
    xdata = clusterCenters(:,1);
    ydata = clusterCenters(:,2);
    pause(1)
    h = plot(xdata, ydata, 'ko', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize', 10);
    set(h,'YDataSource','ydata')
    set(h,'XDataSource','xdata')
    refreshdata
    
    
    % stop criteria: centroids do not change anymore
    if isequal(clusterLabels,oldClusters)
        changes = false;
    else
        oldClusters = clusterLabels;
    end
    
    iter = iter + 1;
end


figure(1)
hold on;

% ploting the final clustering resulting from Fuzzy C-Means
clus = unique(clusterLabels);
colors = {'b.', 'r.', 'k.', 'm.', 'y.', 'c.', 'b.', 'r.', 'k.', 'm.', 'y.', 'c.'};
for i = 1:length(clus)
    indexes = find(clusterLabels==clus(i));
    oldIndex = indexes;
    plot(X(indexes,1), X(indexes,2), colors{i});
    
end
