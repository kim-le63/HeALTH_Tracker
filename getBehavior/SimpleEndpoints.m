function [Skeleton,Endpoints] = SimpleEndpoints(Skeleton)
%This function prunes branches off a skeletal image based on size/fraction of total length
%This is the simpler version

%Skeleton is the binary skeleton of the object with 3+ endpoints
%There must be only one object in view

%Endpoints is a binary image the size of Image, with two true values

ObjectCheck = bwconncomp(Skeleton);
if ObjectCheck.NumObjects > 1
    warning(['More than one object in skeleton sent to FindEndpoints. ', ...
    char(10),'Code may behave unpredictably'])
    Skeleton = bwareaopen(Skeleton,sum(sum(Skeleton))./2);
end

counter = 0;
Endpoints = bwmorph(Skeleton,'endpoints');
NumberOfEndpoints = sum(sum(Endpoints));
[ImRows,ImCols] = size(Skeleton);

while NumberOfEndpoints > 2 && counter < 10

    counter = counter + 1;
    %Create a binary image of only the branches, uniquely labeled
    Branchpoints = bwmorph(Skeleton,'branchpoints');
    BranchpointsPrime = imdilate(Branchpoints,[1,1,1;1,1,1;1,1,1]);
    SkeletonPrime = Skeleton & ~BranchpointsPrime;
    
    %find the size of each branch
    BranchConnComp = bwconncomp(SkeletonPrime,8);
    BranchSizes = cellfun(@numel,BranchConnComp.PixelIdxList);
    LabeledBranches = bwlabel(SkeletonPrime);
    NumBranches = max(max(LabeledBranches));
    [BranchSizes,BranchIndices] = sort(BranchSizes);
    
    %Prune any branch smaller than 30% of the next-smallest branch, or 5% of body
    if NumBranches > 1
        if BranchSizes(1) <= BranchSizes(2)*.3
            Skeleton(BranchIndices(1)) = false;
        end
    end
    BranchIndices(BranchSizes < .05*sum(sum(Skeleton))) = 0;

    %Clean up one-pixel branches and non-branching branchpoints
    Endpoints = bwmorph(Skeleton,'endpoints');
    Skeleton = bwmorph(Skeleton,'spur');
    Skeleton(Endpoints) = true; %spur deletes endpoints
    Skeleton = bwmorph(Skeleton,'thin'); %may be unnecessary
    
    
    NumberOfEndpoints = sum(sum(Endpoints));
end
