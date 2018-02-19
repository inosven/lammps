fprintf('#############################################################\n');
fprintf('#                                                           #\n');
fprintf('# IMG2VOX v2.0                                              #\n');
fprintf('# --- By Yidong Xia & Joshua Kane (Idaho National Lab)      #\n');
fprintf('# --- This MATLAB code is designed to:                      #\n');
fprintf('#     * Find 10 largest pores and create wall layers        #\n');
fprintf('#     * Output these pores & wall as individual data files  #\n');
fprintf('#                                                           #\n');
fprintf('#############################################################\n');

%{
User defines the number of pores they want to index seperately
%}
nRegions = 10;

%{
Open a MATLAB built-in GUI to look for image files.
    * Some specific image types are listed.
    * See help menu for all matlab recognized image types.
    * MultiSelect allows the user to select multiple images.
%}
[fileName, pathName, ~] =...
    uigetfile(...
        {'*.jpg;*.tif;*.tiff;*.png;*.gif;*.bmp;','All Image Files'},...
        'Select stack images to be processed',...
        'MultiSelect','on'...
    );

%{
Determine the number of files to be processed. >= 0 and <= 2^16-1.
Convert file name from cell to string.
Get the pixel resolution of each image.
%}
nFiles = uint16( size(fileName, 2) );
fileName = char(fileName);
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
Delete variables not further used from Workspace
%}
clear fileName nFiles iFile img imgSize ans

%{
Pad volume with a 1 voxel thick layer.
Determine size of volume for later use.
Create pore index and index all pores.
%}
BW3DP = padarray(BW3D,[1,1,1]);
sizeBW3DP = size(BW3DP);
CC = bwconncomp(BW3DP,26);

%{
Get the number of voxels in each connected pore network
ordered by volume from large to small
%}
[~,Idx] = sort(cellfun(@numel,CC.PixelIdxList),'descend');

%{
For the N-th largest pores,
determine the positions of pores and solid "outer shell" individually.
Here the shell is one layer of voxels thick.
Slight modification is needed in the 2nd for loop below
as well as the pad sizes of volumes
%}

%{
We only need the information from the N-th largest pores
%}
DD = CC;
DD.NumObjects = nRegions;
DD.PixelIdxList = [];
for ii = 1 : nRegions
    DD.PixelIdxList{ii} = CC.PixelIdxList{Idx(ii)};
end

%{
No longer need this
%}
clear CC

%{
Extract bounding box and "Cropped Image" of individual pores
%}
stats = regionprops(DD,'BoundingBox','Image');

%{
Pre-allocation
%}
PORE = [];
WALL = [];

%{
Loop over N largest pores
%}
for ii = 1 : nRegions
    %{
    Determine 3D coordinates of all voxels in Pore ii
    %}
    [R_pore,C_pore,Z_pore] = ind2sub(sizeBW3DP,DD.PixelIdxList{ii});
    PORE = [R_pore,C_pore,Z_pore];

    %{
    Determine coordinates of all voxels in solid shell around PORE-ii.
    Pull image of PORE-ii from stats and pad it.
    Extract solid shell around PORE-ii.
    %}
    BWi = padarray(stats(ii).Image,[1,1,1]);
    BWiDil = imdilate(BWi,true(3,3,3));

    %{
    Pull bounding box location and dimensions for Pore ii from stats
    %}
    BB = stats(ii).BoundingBox;

    %{
    Get voxels belonging to solid shell of Pore ii,
    and convert to absolute coordinates of initial loaded volume
    %}
    Peri = BWi ~= BWiDil;
    IndPeri = find(Peri);
    [R_wall,C_wall,Z_wall] = ind2sub(size(BWi),IndPeri);
    R_wall = R_wall + floor(BB(2)) - 1;
    C_wall = C_wall + floor(BB(1)) - 1;
    Z_wall = Z_wall + floor(BB(3)) - 1;
    WALL = [R_wall,C_wall,Z_wall];

    numVoxPore = size(PORE,1);
    numVoxWall = size(WALL,1);

    outFileName = [pathName,['region.' num2str(ii,'%02d') '.info.txt']];
    fileID = fopen(outFileName,'w');
    fprintf(fileID,'numVoxWall\n%g\n',numVoxWall);
    fprintf(fileID,'Rmin Rmax Cmin Cmax Zmin Zmax\n%g %g %g %g %g %g\n',...
        min(R_wall),max(R_wall),...
        min(C_wall),max(C_wall),...
        min(Z_wall),max(Z_wall));
    fprintf(fileID,'numVoxPore\n%g\n',numVoxPore);
    fprintf(fileID,'Rmin Rmax Cmin Cmax Zmin Zmax\n%g %g %g %g %g %g\n',...
        min(R_pore),max(R_pore),...
        min(C_pore),max(C_pore),...
        min(Z_pore),max(Z_pore));
    fclose(fileID);

    outFileName = [pathName,['region.' num2str(ii,'%02d') '.pore.txt']];
    fileID = fopen(outFileName,'w');
    for i = 1 : numVoxPore
        fprintf(fileID,'%g %g %g\n',PORE(i,1),PORE(i,2),PORE(i,3));
    end
    fclose(fileID);

    outFileName = [pathName,['region.' num2str(ii,'%02d') '.wall.txt']];
    fileID = fopen(outFileName,'w');
    for i = 1 : numVoxWall
        fprintf(fileID,'%g %g %g\n',WALL(i,1),WALL(i,2),WALL(i,3));
    end
    fclose(fileID);

    %{
    Remove pore ii from the volume
    %}
    BW3DP(DD.PixelIdxList{ii}) = false;

end

clear BWi BB BWiDil DD
clear Peri IndPeri Idx
clear R_wall C_wall Z_wall
clear R_pore C_pore Z_pore
clear PORE WALL
clear stats

%Unused script
%{
IndP=find(BW3DP);
[R,C,Z]=ind2sub(sizeBW3DP,IndP);
n=numel(R);
C1=ones(n,1)*(nRegions+1);
C2=ones(n,1);
matPore=[C1,C2,R,C,Z];


BW3DPDil=imdilate(BW3DPDil,true(3,3,3));
BW3DPPer=BW3DPDil~=BW3DP;
IndS=find(BW3DPPer);
[R,C,Z]=ind2sub(sizeBW3DP,IndS);
n=numel(R);
C1=ones(n,1)*(nRegions+1);
C2=zeros(n,1);
matWall=[C1,C2,R,C,Z];

matAll=vertcat(matPore,matWall);

filename=[PathName,['Pore' num2str(nRegions+1) 'Coordinate.txt']];
    dlmwrite(filename,matAll)
%}
