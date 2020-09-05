%------------------------------------------------------------------------
%--------------------
% Code to illustrate DNA quantification.


captionFontSize = 14;
% Read in a standard image 
ImageName = 'Image1.png';

% Generates all image information
InfoImage = imfinfo(ImageName);


% Finding out file size from image information
B= InfoImage.FileSize;

% 9085414 is the file size of the standard calibrating image which 
% in the DNA experiments was for 129 ng/ul. This standard image will
% be used in dilutions to generate a standard curve for finding out 
% unknown concentrations of DNA later on.
% Factor F will be used later to calculate filesize based on quantity 
% of DNA

F= 9085414/B; 
 

% Grayscale image and point spread function to remove gaussian noise 


DNAImage = imread(ImageName);
DNAImage2 = imread(ImageName);
I = rgb2gray(DNAImage);

%imadjust saturates the bottom 1% and the top 1% of all pixel values so 
% only the DNA is prominent in the image.
% To separate the edges we saturate the bottom 2% and top 47% of the image.
% Based on the image being analysed this value might differ
% However, for images from same instrument its mostly constant.

I3 = imadjust(DNAImage,[0.02 0.47],[]);
imshow (I3)
 
% Motion Deblurring the image for better DNA quantification 
% Typically DVD images have scratches plus motion blurring due
% to high rotational speeds. 

I4 = imread (ImageName);
PSF = fspecial('gaussian',2,300);
V = .0001;

%imnoise expects pixel values of data type double and single to be 
%in the range [0, 1]. You can use the rescale function to adjust
%pixel values to the expected range. If your image is type double or
%single with values outside the range [0,1], then imnoise clips input
%pixel values to the range [0, 1] before adding noise.

% The imfilter function computes the value of each output pixel 
%using double-precision, floating-point arithmetic. 
% If the result exceeds the range of the data type, then 
% imfilter truncates the result to the allowed range of the data type.
% If it is an integer data type, then imfilter rounds fractional values.

BlurredNoisy = imnoise(imfilter(I,PSF),'gaussian',0,V);
BlurredNoisy = double(BlurredNoisy)/255;

% Description. B = zeros(n) returns an n -by- n matrix of zeros.
% An error message appears if n is not a scalar. B = zeros(m,n) 
% or B = zeros([m n]) returns an m -by- n matrix of zeros.

WT = zeros(size(I));
WT(5:end-4,5:end-4) = 1;
INITPSF = ones(size(PSF));

% deconvolves image I using the maximum likelihood algorithm and an initial 
% estimate of the point- spread function (PSF), psfi.
% The deconvblind function returns both the deblurred image J and a
% restored PSF, psfr.

[J P] = deconvblind(BlurredNoisy,INITPSF,20,10*sqrt(V),WT);
 
figure
imshowpair(J, I3)

% Shows the deblurred image vs the original image
% Now we use the deblurred image for further processing.
% Convert the uint8 image back to double for further processing

% Y = fft(X) computes the discrete Fourier transform (DFT) of X using a
% fast Fourier transform (FFT) algorithm.

% If X is a vector, then fft(X) returns the Fourier transform of the vector.

% If X is a matrix, then fft(X) treats the columns of X as vectors and 
% returns the Fourier transform of each column.

% If X is a multidimensional array, then fft(X) treats the values along
% the first array dimension whose size does not equal 1 as vectors and 
% returns the Fourier transform of each vector.
% double is the default numeric data type (class) in MATLAB®, 
% providing sufficient precision for most computational tasks. 
% Numeric variables are automatically stored as 64-bit (8-byte) 
% double-precision floating-point values. 

DNAImage = J
fft(double(DNAImage))


% Convert the image back to the rgb space for further processing
rgbImage = cat(3, DNAImage, DNAImage, DNAImage);

% Show the color image which will be used for further processing
imshow (rgbImage)

% There are three colors: white, blue, and pink. 
% Notice how easily you can visually distinguish these colors from one another. 
% The L*a*b* color space (also known as CIELAB or CIE L*a*b*) enables you to quantify
% these visual differences.
% The L*a*b* color space is derived from the CIE XYZ tristimulus values.
% The L*a*b* space consists of a luminosity layer 'L*', chromaticity-layer 'a*' 
% indicating where color falls along the red-green axis, and chromaticity-layer 
% 'b*' indicating where the color falls along the blue-yellow axis. 
% All of the color information is in the 'a*' and 'b*' layers.
% You can measure the difference between two colors using the Euclidean distance metric.
% Convert the image to L*a*b* color space using rgb2lab.

lab_he = rgb2lab(rgbImage);
ab = lab_he(:,:,2:3);

% J = im2single(I) converts the grayscale, RGB, or binary image I to single, 
% rescaling or offsetting the data as necessary.
% If the input image is of class single, then the output image is identical. 
% If the input image is of class logical, then im2single changes
% true-valued elements to 65535.

ab = im2single(ab);
nColors = 3;

% Clustering is a way to separate groups of objects.
% K-means clustering treats each object as having a location in space. 
% It finds partitions such that objects within each cluster are as 
% close to each other as possible, and as far from objects in other 
% clusters as possible. K-means clustering requires that you specify 
% the number of clusters to be partitioned and a distance metric to 
% quantify how close two objects are to each other.
% Since the color information exists in the 'a*b*' color space, 
% your objects are pixels with 'a*' and 'b*' values. 
% Convert the data to data type single for use with imsegkmeans. 
% Use imsegkmeans to cluster the objects into three clusters.
% repeat the clustering 3 times to avoid local minima

pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);
imshow(pixel_labels,[])

% Using pixel_labels, you can separate objects in Image by color,
% which will result in three images.The edges of the DNA chamber
% can be removed with only the DNA remaining in the image.

title('Image Labeled by Cluster Index');
mask1 = pixel_labels==1;
cluster1 = DNAImage2 .* uint8(mask1);
imshow(cluster1) 
title('Objects in Cluster 1');

% Second level filtering

mask2 = pixel_labels==2;
cluster2 = DNAImage2 .* uint8(mask2);
imshow(cluster2)
title('Objects in Cluster 2');

% Third level filtering

mask3 = pixel_labels==3;
cluster3 = DNAImage2 .* uint8(mask3);
imshow(cluster3)
title('Objects in Cluster 3');

% Cluster 3 contains the DNA. Notice that there
% are edges and DNA. 
% You can separate edges from DNA precipitate using the 'L*' 
% layer in the L*a*b* color space. The DNA would be darker.

% Recall that the 'L*' layer contains the brightness values of each color.
% Extract the brightness values of the pixels in this cluster 
% and threshold them with a global threshold using imbinarize.
% The mask is_light_blue gives the indices of light blue pixels.

L = lab_he(:,:,1);
L_DNA = L .* double(mask3);
L_DNA = rescale(L_DNA);

% BW = imbinarize(I) creates a binary image from 2-D or 3-D grayscale
% image I by replacing all values above a globally determined
% threshold with 1s and setting all other values to 0s. 
% By default, imbinarize uses Otsu's method, which chooses 
% the threshold value to minimize the intraclass variance of the 
% thresholded black and white pixels [1]. 
% imbinarize uses a 256-bin image histogram to compute
% Otsu's threshold. To use a different histogram, see otsuthresh.

idx_DNA = imbinarize(nonzeros(L_DNA));
DNA_idx = find(mask3);
mask_dark_DNA = mask3;
mask_dark_DNA(DNA_idx(idx_DNA)) = 0;

DNA_nuclei = DNAImage2 .* uint8(mask_dark_DNA);
imshow(DNA_nuclei)
title('DNA part of the image');
hold on



% Step 1: Image check to distinguish between blank and precipitate
% The U shape border on the blank image has a threshold between 85 and 95 
% when generated by the Lab on DVD player
% If the code sees the precipitate on the image then only it proceeds ahead 
% to the next step. Separate threshold range for U shape border in CMOS 
% sensor images. 
 
figure
hold on
imshow(DNAImage)
DNAImage = imcrop(DNAImage)
 
subplot(2, 2, 1);
imshow(DNAImage);
% Maximize the figure window.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
drawnow;
caption = sprintf('Original Image');
title(caption, 'FontSize', captionFontSize);
axis image; 
% Histogram generation and display.
[pixelCount, grayLevels] = imhist(DNAImage);
%bar(grayLevels, pixelCount, 'FaceColor', 'r');
 
% Binary image generation, thresholdValue = 95 for DVD images 
 
 
 
thresholdValue = 95;
 
 
binaryImg =  (DNAImage < thresholdValue);

Sum2=0

for k=85:thresholdValue 
    
    Sum2=Sum2+ pixelCount(k)*((255-grayLevels(k))/255)
  
end
 
binaryImg = imfill(binaryImg, 'holes');
 
% Image to concentration conversion Factor 115 determined experimentally from 
% the images.For the highest concentration 
% 129 ng/ul this factor gives a precipitation area of 1290 a.u. This gives
% a resolution of 1 for every 0.1ng/ul change in DNA concentration.
% This factor is separate for Camera images.
 
Sum= Sum2/115
subplot(2, 2, 2);
bar(pixelCount);
%title('Histogram of DNA image', 'FontSize', captionFontSize);
xlim([0 255]); % Scale x axis manually.
grid off;
 
% Display the threshold as a line on histogram
hold on;
maxYValue = ylim;
line([thresholdValue, thresholdValue], maxYValue, 'Color', 'b');
% Place a text label on the bar chart showing the threshold.
annotationText = sprintf('Thresholded at %d gray levels', thresholdValue(end));
 
% Display the binary image.
subplot(2, 2, 3);
imshow(binaryImg); 
title('Binary Image, obtained by thresholding', 'FontSize', captionFontSize); 
 
%For blank image elimination
 
if (Sum2>20)
 
 
DNAImage2 = imadjust(DNAImage2,[0.45 0.55],[]);
figure
imshow(DNAImage2)
 

% Step 2: Second part if image is not blank
% Crop the image to keep the region containing the DNA precipitate

DNAImage2 = imcrop(DNAImage2)
 
 
% Display the grayscale image.
% Image adjustment: it works for most images from the DVD. This contrast step % makes the DNA look prominent over its background which is necessary for its
% quantification.

subplot(2, 2, 1);
imshow(DNAImage2);

% Maximize the figure window.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
drawnow;
caption = sprintf('Original Image');
title(caption, 'FontSize', captionFontSize);
axis image; 

% Histogram generation and display.
%[counts,binLocations] = imhist(I) 
% calculates the histogram for the grayscale image I. 
% The imhist function returns the histogram counts in counts 
% and the bin locations in binLocations. 
% The number of bins in the histogram is determined by the image type.

[pixelCount, grayLevels] = imhist(DNAImage);
 
% Binary image generation 
s=nonzeros(pixelCount)
A=max(s)
 
% Quantile def from MATLAB resources:
% returns quantiles of the elements in data vector or array X for the
% cumulative probability or probabilities p in the interval [0,1].
% If X is a vector, then Y is a scalar or a vector having 
% the same length as p.
% If X is a matrix, then Y is a row vector or a matrix where the
% number of rows of Y is equal to the length of p.

Y = quantile(s,[0,0.8])
 
 
thresholdValue = find(pixelCount>Y(2),1);
 
 
binaryImg =  (DNAImage < thresholdValue);
Sum2=0

% We use additional 10 pixels here due to different intensities of DNA in
% the image even after contrasting the entire DNA is not of the same 
% intensity but rather lies at a value between 1 to 10 on a scale of 0 to 255.

for k=1:thresholdValue+10
    
    Sum2=Sum2+ pixelCount(k)*((255-grayLevels(k))/255)
  
end
 
 
 
binaryImg = imfill(binaryImg, 'holes');
 
% Factor 115 (Multiply F by suitable value  for Camera images) below is to %ensure 129 ng/ul in 10 ul volume, this factor gives a precipitation area of %1290 a.u. This gives a resolution of 1 for every 0.1ng/ul change in DNA %concentration. So, in essence the value of the precipitation level by %adjusting F should be fixed as the product of known concentration and %volume of sample.
 
%Sum= (Sum2*(F)/115)/10.46 +14 
%Calculating unknown quantity of DNA

Sum= Sum2*F/115

%Calculating precipitation level
%Precipitation level is an arbitrary quantity defined in the context of
%this code to define the amount of DNA in the sample

subplot(2, 2, 2);
bar(pixelCount);

%title('Histogram of DNA image', 'FontSize', captionFontSize);

xlim([0 20]); % Scale x axis manually.
ylim([0 1800])
grid off;
 
% Display the threshold as a line on histogram
hold on;
maxYValue = ylim;
line([thresholdValue, thresholdValue], maxYValue, 'Color', 'b');
% Place a text label on the bar chart showing the threshold.
annotationText = sprintf('Thresholded at %d gray levels', thresholdValue(end));
 
% Display the binary image.
subplot(2, 2, 3);
imshow(binaryImg); 
title('Binary Image, obtained by thresholding', 'FontSize', captionFontSize); 
 

se = offsetstrel('ball',4,4);
dilatedI = imdilate(DNA_nuclei,se);

%Only uncomment for visualisation with original image
%imshowpair(DNA_nuclei,dilatedI,'montage')

%For more prominent visualisation of the thresholded DNA image

subplot(2, 2, 4);
imshow(dilatedI);

title('Nucleic acid part of the image', 'FontSize', captionFontSize); 
axis image;
hold on;

% Color filtering part of the image

he = imread('129c.png');
imshow(he), title('H&E image');
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3);
ab = im2single(ab);

nColors = 3;
% repeat the clustering 3 times to avoid local minima
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);
imshow(pixel_labels,[])
title('Image Labeled by Cluster Index');

mask1 = pixel_labels==1;
cluster1 = he .* uint8(mask1);
imshow(cluster1)
hold on;
title('Objects in Cluster 1');

mask2 = pixel_labels==2;
cluster2 = he .* uint8(mask2);
imshow(cluster2)
hold on;
title('Objects in Cluster 2');


mask3 = pixel_labels==3;
cluster3 = he .* uint8(mask3);
imshow(cluster3)
title('Objects in Cluster 3');


% Observing a frequency response

end
 
