function [BW1,BW2,BW3,wing_reconstruction1]=UserInputReconstructionV3(wing_reconstruction1,DLt_Coef,image1,image2,image3)
%takes in the coordinates of the wings in coordinates (no voxels) and
%returns the user specified edge for extra reconstruction
%% V3: 7-2-2021 (usa time)
%overlays images instead of points
%% user instructions (to be added later)
%% find the wing coordinates in the image
[wing1_uv1] = round(dlt_inverse(DLt_Coef(:,1),wing_reconstruction1));
[wing1_uv2] = round(dlt_inverse(DLt_Coef(:,2),wing_reconstruction1));
[wing1_uv3] = round(dlt_inverse(DLt_Coef(:,3),wing_reconstruction1));
%% convert projections to images
image_projected1=zeros(size(image1));
image_projected2=zeros(size(image2));
image_projected3=zeros(size(image3));
for i=1:length(wing1_uv1)
    image_projected1(wing1_uv1(i,2),wing1_uv1(i,1))=1;
    image_projected2(wing1_uv2(i,2),wing1_uv2(i,1))=1;
    image_projected3(wing1_uv3(i,2),wing1_uv3(i,1))=1;
end
%% user selection
%figures
fig_image1=figure;
fig_image2=figure;
fig_image3=figure;
%image 1
figure(fig_image1)
C1=imfuse(image1,image_projected1);
imshow(C1)
BW1 = roipoly;
%image2
figure(fig_image2)
C2=imfuse(image2,image_projected2);
imshow(C2)
BW2 = roipoly;
%image 3
figure(fig_image3)
C3=imfuse(image3,image_projected3);
imshow(C3)
BW3 = roipoly;
%% find indicies that are outside the respective ROI
ImageSum1=BW1+image_projected1;
ImageSum2=BW2+image_projected2;
ImageSum3=BW3+image_projected3;
[u_1, v_1]=find(ImageSum1==2);
[u_2, v_2]=find(ImageSum2==2);
[u_3, v_3]=find(ImageSum3==2);

[index1,~]=ismember(wing1_uv1,[v_1 u_1],'rows');
[index2,~]=ismember(wing1_uv2,[v_2 u_2],'rows');
[index3,~]=ismember(wing1_uv3,[v_3 u_3],'rows');

index_total=index1+index2+index3;
index_total=find(index_total<3);
wing_reconstruction1(index_total,:)=[];
%% close the figures 
close(fig_image1)
close(fig_image2)
close(fig_image3)