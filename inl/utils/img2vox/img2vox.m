fprintf('#############################################################\n');
fprintf('#                                                           #\n');
fprintf('# IMG2VOX - a MATLAB code to convert stack images to voxels #\n');
fprintf('#                                                           #\n');
fprintf('# by Yidong Xia, Joshua Kane (Idaho National Laboratory)    #\n');
fprintf('#                                                           #\n');
fprintf('# This code does the following tasks:                       #\n');
fprintf('#   * convert stack images to voxel data                    #\n');
fprintf('#   * porosity distribution analysis                        #\n');
fprintf('#                                                           #\n');
fprintf('#############################################################\n');

%{
Open a MATLAB built-in GUI to look for image files.
    * Some specific image types are listed.
    * See help menu for all matlab recognized image types.
    * MultiSelect allows the user to select multiple images.

Variables defined in the function
     fileName: the actual name of the file, e.g. "example.jpg";
               has class cell
     pathName: the file path, e.g. "C:\Users\KaneJJ\Desktop\";
               has class string
%}
[fileName, pathName, ~] =...
    uigetfile(...
        {'*.jpg;*.tif;*.tiff;*.png;*.gif;*.bmp;','All Image Files'},...
        'Select stack images to be processed',...
        'MultiSelect','on'...
    );

%{
Determine the number of files to be processed.
    * nFiles must be >= 0 and <= 2^16-1
%}
nFiles = uint16( size(fileName, 2) );

%{
Convert file name from cell to string.
%}
fileName = char(fileName);

%{
Get the pixel resolution of each image
%}
img = imread( [pathName fileName(1,:)] );
imgSize = size(img);

%{
Initiate a 3D matrix for the stack images
%}
if ( isa(img, 'double') )
    BW3D = zeros( imgSize(1), imgSize(2), nFiles );
elseif ( isa(img, 'single') )
    BW3D = zeros( imgSize(1), imgSize(2), nFiles, 'single' );
elseif ( isa(img, 'uint16') )
    BW3D = zeros( imgSize(1), imgSize(2), nFiles, 'uint16' );
elseif ( isa(img, 'uint8') )
    BW3D = zeros( imgSize(1), imgSize(2), nFiles, 'uint8' );
else
    BW3D = false( imgSize(1), imgSize(2), nFiles );
end

%{
Serial-reading of images to assign values in BW3D matrix
%}
for iFile=1:nFiles
    BW3D(:,:,iFile) = imread( [pathName fileName(iFile,:)] );
end

%{
Get the number of all voxels and the number of pore voxels
%}
numVoxAll = numel(BW3D);
numVoxPore = sum(sum(sum(BW3D)));

%{
Delete variables not further used from Workspace
%}
clear fileName nFiles iFile img imgSize ans

%{
Calculate connectivity with the default 26-connectivity
%}
CC = bwconncomp(BW3D,26);

%{
Get the number of voxels in each connected pore network
ordered from large to small
%}
[NumVoxConnPores,Idx] = sort(cellfun(@numel,CC.PixelIdxList),'descend');

%{
Display some useful information
%}
fprintf('Resolution in the direction of Row, Column and Page:\n');
fprintf('  %d %d %d\n',size(BW3D,1),size(BW3D,2),size(BW3D,3));

fprintf('Ratio of voxels -- 5 largest pore-networks vs. all pores:\n');
for i = 1:5
    fprintf('  %5.2f%%', 100 * NumVoxConnPores(i) / numVoxPore);
end
fprintf('\n');

fprintf('Ratio of voxels -- 5 largest pore-networks vs. input domain:\n');
for i = 1:5
    fprintf('  %5.2f%%', 100 * NumVoxConnPores(i) / numVoxAll);
end
fprintf('\n');

%{
Create a 2D matrix for coloring the pore voxels.
Each row of the matrix stores (R, C, Z, Color)
%}
VOX = zeros(numVoxPore, 4);

%{
Color the #1 largest pore-network voxels with 1
%}
[R,C,Z] = ind2sub( size(BW3D), CC.PixelIdxList{Idx(1)} );
numStart = 1;
numEnd = size(R,1);
for ii = numStart : numEnd
    iii = ii - numStart + 1;
    VOX(ii,1) = R(iii);
    VOX(ii,2) = C(iii);
    VOX(ii,3) = Z(iii);
    VOX(ii,4) = 1;
end
fprintf('Voxels of #01 pore-network: %12d\n',iii);

%{
Color the #2 - #10 largest pore-network voxels from 2 to 10
%}
for i = 2 : 10
    [R,C,Z] = ind2sub( size(BW3D), CC.PixelIdxList{Idx(i)} );
    numStart = numEnd + 1;
    numEnd = numEnd + size(R,1);
    for ii = numStart : numEnd
        iii = ii - numStart + 1;
        VOX(ii,1) = R(iii);
        VOX(ii,2) = C(iii);
        VOX(ii,3) = Z(iii);
        VOX(ii,4) = i;
    end
    fprintf('Voxels of #%02d pore-network: %12d\n',i,iii);
end

%{
Color the rest of pore voxels with 11
%}
for i = 11 : CC.NumObjects
    [R,C,Z] = ind2sub( size(BW3D), CC.PixelIdxList{Idx(i)} );
    numStart = numEnd + 1;
    numEnd = numEnd + size(R,1);
    for ii = numStart : numEnd
        iii = ii - numStart + 1;
        VOX(ii,1) = R(iii);
        VOX(ii,2) = C(iii);
        VOX(ii,3) = Z(iii);
        VOX(ii,4) = 11;
    end
end

%{
Obsolete

Data = [];
for i = 1 : 10
    [R,C,Z] = ind2sub( size(BW3D), CC.PixelIdxList{Idx(i)} );
    N = ones(size(R)) * i;
    Data = [Data; [R,C,Z,N]];
    msg = sprintf('Voxels of #%02d pore-network: %8d',i,size(N,1));
    disp(msg);
end

for i = 11 : CC.NumObjects
    [R,C,Z] = ind2sub( size(BW3D), CC.PixelIdxList{Idx(i)} );
    N = ones(size(R)) * 11;
    Data = [Data; [R,C,Z,N]];
end
%}

tStart = cputime;
outFileName = "vox_pore.txt";
fprintf('Export file %s ... ', outFileName);
fileID = fopen(outFileName,'w');
for i = 1 : numVoxPore
    fprintf(fileID,'%g %g %g %g\n',VOX(i,1),VOX(i,2),VOX(i,3),VOX(i,4));
end
fclose(fileID);
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);

tStart = cputime;
fprintf('Dilate the pore voxels with direct neighbor solid voxels... ');
DIL1 = imdilate(BW3D,true(3,3,3));
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);

tStart = cputime;
fprintf('Get the surface wall voxels that encolse the pores... ');
SURF1 = DIL1 ~= BW3D;
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);
clear DIL1

tStart = cputime;
fprintf('Get the indices of all non-zero entries... ');
indSURF1 = find(SURF1);
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);
clear SURF1

tStart = cputime;
fprintf('Convert the non-zero entry indices to I,J,K indexing... ');
[R,C,Z] = ind2sub(size(BW3D),indSURF1);
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);
clear BW3D indSURF1

tStart = cputime;
fprintf('Append the boundary wall voxels into the pore voxel matrix... ');
VOX = [[R,C,Z,zeros(size(R))]; VOX];
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);

tStart = cputime;
outFileName = "vox_wall_pore.txt";
fprintf('Export file %s ... ', outFileName);
fileID = fopen(outFileName,'w');
numVoxWallPore = size(VOX,1);
for i = 1 : numVoxWallPore
    fprintf(fileID,'%g %g %g %g\n',VOX(i,1),VOX(i,2),VOX(i,3),VOX(i,4));
end
fclose(fileID);
tElapsed = cputime - tStart;
fprintf('CPU time = %.3f\n', tElapsed);

%Dil2 = imdilate(Dil1,true(3,3,3));
%Surf2=Dil2~=Dil1;
