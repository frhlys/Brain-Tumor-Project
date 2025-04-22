%%
clc;
clear all;
close all;

% Input
[I, path] = uigetfile('*.jpg', 'Select FLAIR image');
str = strcat(path, I);
s = imread(str);

[I, path] = uigetfile('*.jpg', 'Select T1 image');
str = strcat(path, I);
s2 = imread(str);

[I, path] = uigetfile('*.jpg', 'Select T1CE image');
str = strcat(path, I);
s3 = imread(str);

[I, path] = uigetfile('*.jpg', 'Select T2 image');
str = strcat(path, I);
s4 = imread(str);

% Filtering 
num_iter = 10;
delta_t = 1 / 7;
kappa = 15;
option = 2;
disp('Preprocessing image, please wait...');

inp = imgaussfilt(s, 3); % You can adjust the standard deviation value (2 in this case)
inp = uint8(inp);

inp2 = imgaussfilt(s2, 3);
inp2 = uint8(inp2);

inp3 = imgaussfilt(s3, 3);
inp3 = uint8(inp3);

inp4 = imgaussfilt(s4, 3);
inp4 = uint8(inp4);

inp = imresize(inp, [256, 256]);
if size(inp, 3) > 1
    inp = rgb2gray(inp);
end

inp2 = imresize(inp2, [256, 256]);
if size(inp2, 3) > 1
    inp2 = rgb2gray(inp2);
end

inp3 = imresize(inp3, [256, 256]);
if size(inp3, 3) > 1
    inp3 = rgb2gray(inp3);
end

inp4 = imresize(inp4, [256, 256]);
if size(inp4, 3) > 1
    inp4 = rgb2gray(inp4);
end

figure;
subplot(2, 2, 1);
imshow(inp);
title('Filtered FLAIR Image');

subplot(2, 2, 2);
imshow(inp2);
title('Filtered T1 Image');

subplot(2, 2, 3);
imshow(inp3);
title('Filtered T1CE Image');

subplot(2, 2, 4);
imshow(inp4);
title('Filtered T2 Image');
pause(2);

% Compute the gradient of the filtered images
[gx, gy] = imgradientxy(inp);
gx2 = imgradientxy(inp2);
gy2 = imgradientxy(inp2);
gx3 = imgradientxy(inp3);
gy3 = imgradientxy(inp3);
gx4 = imgradientxy(inp4);
gy4 = imgradientxy(inp4);

% Compute divergence for each image
div_result = compute_divergence(gx, gy);
div_result2 = compute_divergence(gx2, gy2);
div_result3 = compute_divergence(gx3, gy3);
div_result4 = compute_divergence(gx4, gy4);

% Update the input images using the computed divergence
inp = double(inp) + delta_t * div_result;
inp2 = double(inp2) + delta_t * div_result2;
inp3 = double(inp3) + delta_t * div_result3;
inp4 = double(inp4) + delta_t * div_result4;


% Display the updated images after divergence update
figure;
subplot(2, 2, 1);
imshow(inp);
title('Filtered FLAIR Image (After Divergence)');

subplot(2, 2, 2);
imshow(inp2);
title('Filtered T1 Image (After Divergence)');

subplot(2, 2, 3);
imshow(inp3);
title('Filtered T1CE Image (After Divergence)');

subplot(2, 2, 4);
imshow(inp4);
title('Filtered T2 Image (After Divergence)');

%Thresholding
sout = imresize(inp,[256,256]);
t0 = 70;
th = t0 + ((max(inp(:))+ min(inp(:)))./2);
for i=1:1:size(inp,1)
    for j=1:1:size(inp,2)
        if inp(i,j)>th
            sout(i,j)=1;
        else
            sout(i,j)=0;
        end
    end
end

sout2 = imresize(inp2,[256,256]);
t0=70;
th = t0 + ((max(inp2(:)) + min(inp2(:)))./2);

for i=1:1:size(inp2,1)
    for j=1:1:size(inp2,2)
        if inp2(i,j)>th
            sout2(i,j)=1;
        else
            sout2(i,j)=0;
        end
    end
end

sout3 = imresize(inp3,[256,256]);
t0=70;
th = t0 + ((max(inp3(:)) + min(inp3(:)))./2);

for i=1:1:size(inp3,1)
    for j=1:1:size(inp3,2)
        if inp3(i,j)>th
            sout3(i,j)=1;
         else
            sout3(i,j)=0;
        end
    end
end

sout4 = imresize(inp4,[256,256]);
t0 = 120;
th = t0 + ((max(inp4(:)) + min(inp4(:)))./2);

for i=1:1:size(inp4,1)
    for j=1:1:size(inp4,2)
        if inp4(i,j)>th
            sout4(i,j)=1;
        else
            sout4(i,j)=0;
        end
    end
end

%Morphological Operation
bw = im2bw(sout, 0.7); % Thresholding the image to binary
label = bwlabel(bw); % Label connected components
stats=regionprops(label,'Solidity', 'Area', 'BoundingBox');
density = (stats.Solidity);
area = (stats.Area);
high_dense_area = density > 0.5;
max_area = max(area(high_dense_area));
tumor_label = find(area==max_area);
tumor = ismember(label, tumor_label);
if max_area>100
    figure;
    imshow(tumor)
    title('Tumor Alone', 'FontSize', 20);
else
    h = msgbox('No Tumor!!', 'status');
    return;
end

bw2 = im2bw(sout2, 0.7);
label2 = bwlabel(bw2);
stats2 = regionprops(label2, 'Solidity', 'Area', 'BoundingBox');
density2 = (stats.Solidity);
area2 = (stats2.Area);
high_dense_area2 = density2 > 0.5;
max_area2 = max(area2(high_dense_area2));
tumor_label2 = find(area2==max_area2);
tumor2 = ismember(label2, tumor_label2);
if max_area>100
    figure;
    imshow(tumor2)
    title('Tumor Alone', 'FontSize', 20);
else
    h = msgbox('No Tumor!!', 'status');
    return;
end

bw3 = im2bw(sout3, 0.7);
label3 = bwlabel(bw3);
stats3 = regionprops(label3, 'Solidity','Area','BoundingBox');
density3 = [stats3.Solidity];
area3 = [stats3.Area];
high_dense_area3 = density3 > 0.5;
max_area3 = max(area3(high_dense_area3));
tumor_label3 = find(area3==max_area3);
tumor3 = ismember(label3, tumor_label3);

se = strel('square', 5);
tumor3 = imdilate(tumor3, se);


if max_area>100
    figure;
    imshow(tumor3)
    title('Tumor Alone', 'FontSize', 20);
else
    h = msgbox('No Tumor!!', 'status');
    return;
end

bw4 = im2bw(sout4, 0.7);
label4 = bwlabel(bw4);
stats4 = regionprops(label4, 'Solidity', 'Area', 'BoundingBox');
density4 = [stats4.Solidity];
area4 = [stats4.Area];
high_dense_area4 = density4 > 0.5;
max_area4 = max(area4(high_dense_area4));
tumor_label4 = find(area4==max_area4);
tumor4 = ismember(label4, tumor_label4);

output = (double(tumor) + double(tumor2) + double(tumor3) + double(tumor4))/4;
imshow(output)
title('SEGMENTED IMAGE','FontSize', 20);

display All
figure
subplot(231);imshow(s);title('FLAIR IMAGE', 'FontSize', 14);
subplot(232);imshow(s2);title('T1 IMAGE', 'FontSize', 14);
subplot(233);imshow(s3);title('T1CE IMAGE', 'FontSize', 14);
subplot(234);imshow(s4);title('T2 IMAGE', 'FontSize', 14);
subplot(235); imshow(output); title('OUTPUT IMAGE', 'FontSize', 14);

% Function to compute divergence
function div_result = compute_divergence(gx, gy)
    div_result = gradient(gx) + gradient(gy);
end