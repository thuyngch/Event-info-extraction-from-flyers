%------------------------------------------------------------------------------
%	Date 		: Oct 25, 2017
%	Description :
%		This file is used to initialize the project.
%------------------------------------------------------------------------------


%% Clean
clear
clc


%% Subfolder
addpath('fnc')
addpath('im')
addpath('result')
addpath('refine')
addpath('gui')


%% OCR engine
ocr_name = 'tesseract';
save('ocr.mat', 'ocr_name')

