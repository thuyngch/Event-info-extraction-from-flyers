%------------------------------------------------------------------------------
%	Date 		: Oct 25, 2017
%	Description :
%		This file is used to run the old algorithm to extract the key information
%	from a flyer image.
%------------------------------------------------------------------------------


%% Clean
close all
clear
clc
load ocr.mat
load dt_set.mat


%% User changes this variable to the number of image that it wants to examine.
im_number = 20;


%% Remove previous sub-images and text file
cd result
delete(sprintf('%d.txt', im_number))
delete *.jpg
cd ../


%% Setup paths
input_img_path = ['im/capture/', num2str(im_number), '.jpg'];
output_txt_path = ['result/', num2str(im_number), '.txt'];


%% Execute
RecognizeText(input_img_path, output_txt_path);

