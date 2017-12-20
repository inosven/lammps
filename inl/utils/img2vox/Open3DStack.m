[FileName, PathName, ~]=uigetfile({'*.jpg;*.tif;*.tiff;*.png;*.gif;*.bmp;',...
    'All Image Files' }, 'Select stack images to be processed',...
    'MultiSelect','on');
    %This says open a matlab built in gui to look for image files
    %(Some specific image types are listed) See help menu for all matlab
    %recognized image types.  MultiSelect simply allows the user to select
    %multiple images.

%Variables defined in above function
%FileName: the actual name of the file for example "example.jpg"
%-FileName has class cell
%PathName: the file path for example "C:\Users\KaneJJ\Desktop\"
%-PathName has class string

%Determine the number of files to be processed
n=uint16(size(FileName,2));% Note that n must be >=0 but <=2^16-1

%Convert FileName from cell to string
FileName=char(FileName);

%%%%%%%%%%%%%%%%%%%%%%
%comment these lines if you don't want to run everything in parallel
% isPoolOpen=isempty(gcp('nocreate'));
% 
% if(isPoolOpen)
%    parpool('local',32)%Most computers don't have 32 cpus to use might need to change this
% end
%%%%%%%%%%%%%%%%%%%%%%%%%

I0=imread([PathName FileName(1,:)]);
s=size(I0);
 
if (isa(I0,'double'))
    BW3D=zeros(s(1),s(2),n);
elseif (isa(I0,'single'))
    BW3D=zeros(s(1),s(2),n,'single');
elseif (isa(I0,'uint16'))
    BW3D=zeros(s(1),s(2),n,'uint16');
elseif (isa(I0,'uint8'))
    BW3D=zeros(s(1),s(2),n,'uint8');
else    
    BW3D=false(s(1),s(2),n);
end

%  parfor i=1:n %change to for if you don't want to run in parallel
%      FileIn=[PathName FileName(i,:)];
%      I=imread(FileIn);
%      BW3D(:,:,i)=I;
%  end

 for i=1:n %change to for if you don't want to run in parallel
     FileIn=[PathName FileName(i,:)];
     I=imread(FileIn);
     BW3D(:,:,i)=I;
 end

clear n NumWorkers isPoolOpen s FileIn I0 ans