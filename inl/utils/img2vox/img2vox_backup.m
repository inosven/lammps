%-- 12/18/2017 05:13:29 PM --%
Open3DStack
clc
Open3DStack
CC=bwconncomp(BW3D,26)
NumVox=sort(cellfun(@numel,CC.PixelIdxList));
NumVox(1:20)
NumVox=sort(cellfun(@numel,CC.PixelIdxList),'descend');
NumVox(1:20)
sum(sum(sum(BW3D)))
100*NumVox(1:20)/ans
histogram(100*NumVox/sum(sum(sum(BW3D))))
100*NumVox(end-20:end)/sum(sum(sum(BW3D)))
0.2085*sum(sum(sum(BW3D)))
0.2085*sum(sum(sum(BW3D)))/100
numel(BW3D)
clc
[~,Idx]=sort(cellfun(@numel,CC.PixelIdxList),'descend');
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(1)});
mean(R)
mean(C)
mean(Z)
LargePore=false(size(BW3D));
LargePore(CC.PixelIdxList{Idx(1)})=true;
subplot(1,2,1)
imshow(BW3D(:,:,1))
subplot(1,2,2)
imshow(LargePore(:,:,1))
subplot(1,2,1)
imshow(BW3D(:,:,end))
subplot(1,2,2)
imshow(LargePore(:,:,end))
subplot(2,2,4)
imshow(LargePore(:,:,1))
min(R)
max(R)
min(C)
max(C)
min(Z)
max(Z)
LargePoreBB=LargePore(min(R):max(R),min(C):max(C),min(Z):max(Z));
scatter3(R,C,Z,'.')
scatter3(R,C,Z,'.',0.5)
scatter3(R,C,Z,'.',1)
clc
scatter3(R,C,Z,'.')
hold on
for i=1:10
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
scatter3(R,C,Z,,'.')
end
for i=1:10
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
scatter3(R,C,Z,'.')
end
clc
memory
mem
pmem
clc
CC=bwconncomp(BW,26);
CC=bwconncomp(BW3D,26);
[NumVox,Idx]=sort(cellfun(@numel,CC.PixelIdxList),'descend');
Idx(1)
Data=[];
for i=1:10
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
Data=[Data;[R,C,Z]]
end
clc
Data=[];
for i=1:10
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
Data=[Data;[R,C,Z]];
end
Data=[];
for i=1:10
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
N=ones(size(R))*i;
Data=[Data;[R,C,Z]];
end
Data=[];
for i=1:10
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
N=ones(size(R))*i;
Data=[Data;[R,C,Z,N]];
end
size(Data)
for i=11:CC.NumObjects
[R,C,Z]=ind2sub(size(BW3D),CC.PixelIdxList{Idx(i)});
N=ones(size(R))*11;
Data=[Data;[R,C,Z,N]];
end
Dil1=imdilate(BW3D,true(3,3,3));
Surf1=Dil1~=BW3D;
Dil2=imdilate(Dil1,true(3,3,3));
Surf2=Dil2~=Dil1;
imshow(Surf1(;,:,round(end/2));
imshow(Surf1(:,:,round(end/2));
imshow(Surf1(:,:,round(end/2)))
IndS1=find(Surf1);
IndS2=find(Surf2);
[R,C,Z]=ind2sub(size(BW3D),IndS1);
DataS1=[[R,C,Z,zeros(size(R))];Data];
[R,C,Z]=ind2sub(size(BW3D),IndS2);
DataS2=[[R,C,Z,zeros(size(R))];Data];
csvwrite(DataS1.csv, DataS1);
csvwrite('DataS1.csv', DataS1);
csvwrite('Data.csv', Data);
