%------------------------------------------------------------------------------
%	Date 		: Oct 20, 2017
%	Description :
%		This file is used to extract the key information from a flyer image 
%	using the improved algorithm.
%------------------------------------------------------------------------------
function [str, time] = recognizeText(input_img_path, output_txt_path, ocr_name, dataset, flgPlot)

%% Global parameters
thr_canny = 0.1;
sz_canny = 3;
sig_canny = 1;

hz_ang = 15;
vt_ang = 15;

wnd_sz  = 3;
edg_den_low = 0.2;
edg_den_up = 0.5;

thr_area = [1/20*1/15, 1/2];
thr_oriAng = 10;
thr_solidity = 0.65;
thr_text_hei = 20;

plotBound = flgPlot(1);		% Plot the boundary of the image
plotHomo = flgPlot(2);		% Plot the retified image
plotObj = flgPlot(3);		% Plot the text-contain regions


%% Read the image
im = imread(input_img_path);
if(size(im, 3)==3)
    im_gr = rgb2gray(im);
end
im_gr = im2double(im_gr);


%% Detect horizontal lines and vertical lines
tic
% Intialize
b_found_lines = [false false; false false]; %hz_1,hz_2,vt_1,vt_2
thetas_mm = [0 0; 0 0];
rhos_mm = [0 0; 0 0];
[hei, wid] = size(im_gr);

% Hough transform
im_gr_med = medfilt2(im_gr, [7,7]);
edge_bw = edge(im_gr_med, 'canny', thr_canny, sig_canny);
edge_bw2 = imdilate(edge_bw, strel('diamond',3));
[H, theta, rho] = hough(edge_bw2, 'ThetaResolution', 0.5, 'RhoResolution', 4);

% Range of Rho and Theta in Horizontal and Vertical
t_range_vt = theta>-vt_ang & theta<vt_ang;
t_range_hz = theta > (90-hz_ang) | theta < (hz_ang-90);
r_range_vt = rho>=-0.3*wid & rho<=1.3*wid;
r_range_hz = rho>=-hei & rho<=hei;

% Limit the result of Hough transform (1-horizontal, 2-vertical)
H1 = H(r_range_hz, t_range_hz);
H2 = H(r_range_vt, t_range_vt);

T1 = theta(:, t_range_hz);
T2 = theta(:, t_range_vt);

R1 = rho(:, r_range_hz);
R2 = rho(:, r_range_vt);

% Scan for horizontal and vertical direction
line_quad = zeros(8,2);
d_rhos = zeros(2,1);
for j=1:2
    % Choose the hor or vert line
    if j==1         % horizontal line
        hou = H1;
        Tt = T1;
        Rt = R1;
    else            % vertical line
        hou = H2;
        Tt = T2;
        Rt = R2;
    end

    % Find the lines of hough transform
    P  = houghpeaks(hou, 16, 'threshold', ceil(0.3*max(hou(:))));
    lines = houghlines(edge_bw2, Tt, Rt, P, 'FillGap', 5, 'MinLength', 1/5*size(im,3-j));

    % No line detected
    max_len = 0;
    rho_max = 0;
    rho_min = 0;
    if isempty(fieldnames(lines))
        b_found_lines(j,:) = false;

    % Detected lines
    else
        % Determine the longest lines
        nline = length(lines);
        for k = 1:nline
            xy = [lines(k).point1; lines(k).point2];
            len = norm(lines(k).point1 - lines(k).point2);
            if len > max_len
                max_len = len;
            end
        end

        % Identify 4 corners, just find the one with min and max rho.
        rhos = zeros(nline,1);
        thetas = zeros(nline,1);
        center_img = [size(im,2), size(im,1)] / 2;
        
        % Horizontal line, rho can be pos or neg
        if j==1
            for k=1:nline
                rho = lines(k).rho;
                the = deg2rad(lines(k).theta); 
                % Calculate Rhos from the center of the image
                thetas(k) = the;
                rhos(k) = rho - (center_img(1)*cos(the)+center_img(2)*sin(the));
                if sin(the) < 0
                    rhos(k) = -rhos(k);
                    thetas(k) = the+pi;
                end
            end

        % Vertical line, rho is guaranted to be positive
        else 
            for k=1:nline
                rho = lines(k).rho;
                the = deg2rad(lines(k).theta);
                % Calculate Rhos from the center of the image
                thetas(k) = the;          
                rhos(k) = rho - (center_img(1)*cos(the)+center_img(2)*sin(the));
            end
        end

        % Don't take rhos that are too close to boundary
        thr_bnd_ratio_max = 0.95;
        thr_bnd_ratio_min = 0;
        [rhos_asc, ind_asc] = sort(rhos);
        hf_dim = size(im,j)/2;

        % 1st line: Min rho(-)
        for k=1:length(ind_asc)
            if rhos_asc(k) > -thr_bnd_ratio_max*hf_dim && ...
                    rhos_asc(k) < -thr_bnd_ratio_min*hf_dim
                ind_min = ind_asc(k);
                rho_min = rhos_asc(k);
                b_found_lines(j,1) = true;
                thetas_mm(j,1) = thetas(ind_min);
                break;
            end
        end

        % If still have the edge points, take the ones on the vertical edge
        if ~b_found_lines(j,1)
            for k = length(ind_asc):-1:1       
                if rhos_asc(k) < -thr_bnd_ratio_max*hf_dim
                    ind_min = ind_asc(k);
                    rho_min = rhos_asc(k);
                    b_found_lines(j,1) = true;
                    thetas_mm(j,1) = thetas(ind_min);
                    break;
                end
            end
        end
        
        % 
        if (b_found_lines(j,1))
            line_quad(1+4*(j-1),:) = lines(ind_min).point1;
            line_quad(2+4*(j-1),:) = lines(ind_min).point2;
        end

        % 2nd line: Max rho(+)
        for k=length(ind_asc):-1:1
            if rhos_asc(k) < thr_bnd_ratio_max*hf_dim && ...
                    rhos_asc(k) > thr_bnd_ratio_min*hf_dim
                ind_max = ind_asc(k);
                rho_max = rhos_asc(k);
                b_found_lines(j,2) = true;
                thetas_mm(j,2) = thetas(ind_max);
                break;
            end
        end
        
        % If still have the edge points, take the ones on the vertical edge
        if (b_found_lines(j,2)==false)
            for k = 1:length(ind_asc)       
                if rhos_asc(k) > thr_bnd_ratio_max*hf_dim
                    ind_max = ind_asc(k);
                    rho_max = rhos_asc(k);
                    b_found_lines(j,2) = true;
                    thetas_mm(j,2) = thetas(ind_max);
                    break;
                end
            end
        end
        
        % 
        if (b_found_lines(j,2) == true)
            line_quad(3+4*(j-1),:) = lines(ind_max).point1;
            line_quad(4+4*(j-1),:) = lines(ind_max).point2;
        end

        % 
        d_rhos(j) = rho_max - rho_min;
        rhos_mm(j,:) = [rho_min, rho_max];

    end %if lines struture is not empty.
end %j=1:2


%% Figure out if there are missing edges
bh1 = b_found_lines(1,1);
bh2 = b_found_lines(1,2);

bv1 = b_found_lines(2,1);
bv2 = b_found_lines(2,2);

dist = 200;
nEdgesFound = sum(b_found_lines(:));

% Need to fill the edges
if nEdgesFound < 4

    % Check horizontal line first
    if bh1+bh2 < 2
        hf_dim = hei/2;
        if ~bh1 && bh2
            rhos_mm(1,1) = min(-0.7*hf_dim, -rhos_mm(1,2));
            thetas_mm(1,1) = thetas_mm(1,2);
        elseif bh1 && ~bh2
            rhos_mm(1,2) = max(0.7*hf_dim, -rhos_mm(1,1));        
            thetas_mm(1,2) = thetas_mm(1,1);
        elseif ~bh1 && ~bh2
            rhos_mm(1,1) = -hf_dim * thr_bnd_ratio_max;
            rhos_mm(1,2) = -hf_dim * thr_bnd_ratio_max;
            thetas_mm(1,1) = pi;
            thetas_mm(1,2) = pi;
        end
        d_rhos(1) = rhos_mm(1,2) - rhos_mm(1,1);
        [pt1, pt2] = line_from_rho_theta(rhos_mm(1,1),thetas_mm(1,1),dist);
        [pt3, pt4] = line_from_rho_theta(rhos_mm(1,2),thetas_mm(1,2),dist);
        line_quad(1,:) = pt1+center_img; line_quad(3,:) = pt3+center_img;
        line_quad(2,:) = pt2+center_img; line_quad(4,:) = pt4+center_img;
    end

    % Now check the vertical edges
    if bv1+bv2 < 2
        theta_hz = 0.5*(thetas_mm(1,2)+thetas_mm(1,1));
        hf_dim = wid/2;
        if (~bv1&&bv2) %Only one edge is available
            thetas_mm(2,1) = theta_hz*2-pi-thetas_mm(2,2);
            rhos_mm(2,1) = min(-rhos_mm(2,2),-0.7*hf_dim);
        elseif (bv1&&~bv2)
            thetas_mm(2,2) = theta_hz*2-pi-thetas_mm(2,1);
            rhos_mm(2,2) = max(-rhos_mm(2,1),0.7*hf_dim);
        elseif (~bv1 && ~bv2)
            %Just make the angle perpendicular to the horizontal edges
            theta = 0.5*(thetas_mm(1,2)+thetas_mm(1,1))-0.5*pi;
            thetas_mm(2,1) = theta; thetas_mm(2,2) = theta;
            rhos_mm(2,1) = -hf_dim*thr_bnd_ratio_max;
            rhos_mm(2,2) = -rhos_mm(2,1);
        end
        d_rhos(2) = rhos_mm(2,2) - rhos_mm(2,1);
        [pt1, pt2] = line_from_rho_theta(rhos_mm(2,1),thetas_mm(2,1),dist);
        [pt3, pt4] = line_from_rho_theta(rhos_mm(2,2),thetas_mm(2,2),dist);
        line_quad(5,:) = pt1+center_img; line_quad(7,:) = pt3+center_img;
        line_quad(6,:) = pt2+center_img; line_quad(8,:) = pt4+center_img;             
    end
end

% Find the intersection of the four lines
ln_hz_1 = line_quad(1:2,:);
ln_hz_2 = line_quad(3:4,:);
ln_vt_1 = line_quad(5:6,:);
ln_vt_2 = line_quad(7:8,:);

pt1 = linlinintersect([ln_hz_1; ln_vt_1]);
pt2 = linlinintersect([ln_hz_1; ln_vt_2]);
pt3 = linlinintersect([ln_hz_2; ln_vt_2]);
pt4 = linlinintersect([ln_hz_2; ln_vt_1]);

if plotBound
    figure(1)
    pt_in = [pt1; pt2; pt3; pt4];
    imshow(edge_bw2);
    hold on
    plot(pt_in(:,1), pt_in(:,2), 'xr', 'LineWidth', 3)
    line([pt1(1); pt2(1)], [pt1(2); pt2(2)], 'LineWidth', 2, 'Color', 'Green')
    line([pt1(1); pt4(1)], [pt1(2); pt4(2)], 'LineWidth', 2, 'Color', 'Green')
    line([pt3(1); pt4(1)], [pt3(2); pt4(2)], 'LineWidth', 2, 'Color', 'Green')
    line([pt2(1); pt3(1)], [pt2(2); pt3(2)], 'LineWidth', 2, 'Color', 'Green')
    title('Boundary on edge map')
end


%% Homographic transform, but keep the aspect ratio
pt_in = [pt1; pt2; pt3; pt4];
pt_out = [0 0; d_rhos(2) 0; d_rhos(2) d_rhos(1); 0 d_rhos(1)];
Hom = estimate_homography_simp(pt_in', pt_out');
tf = maketform('projective', Hom');
im_gr = imtransform(im_gr, tf, 'bilinear','XData',[1 d_rhos(2)], 'YData', [1 d_rhos(1)]);
im_gr = medfilt2(im_gr, [3,3]);

if plotHomo
    figure(2)
    imshow(im_gr)
    title('Homographic-transformed image')
end
time.homo = toc;


%% Detect blocks
tic
% Calcualte variance image
boxflt  = fspecial('average', [wnd_sz, wnd_sz]);
img_ex  = imfilter(im_gr, boxflt, 'replicate');
img_ex2 = imfilter(im_gr .* im_gr, boxflt, 'replicate');
img_var = img_ex2 - img_ex .* img_ex;
thr_imgvar = prctile(img_var(:), 90) * 0.08; 	% Avoid extream values
imgthr_var = img_var;
imgthr_var(img_var < thr_imgvar) = 0;

% Compute the edge density, real text should have a edge density not too big nor small
edge_den = edgeCanny(im_gr, sz_canny, sig_canny);
edge_den(edge_den<edg_den_low | edge_den>edg_den_up) = 0;

% Get the mask region
msk_region = edge_den>0 & imgthr_var>0;

% Label regions and calculate region property
msk_region = imdilate(msk_region, ones(31,31));
regLabels = bwlabel(msk_region, 8);

if plotObj
    figure(3)
    imshow(regLabels)
    title('Candidate regions')
    hold on
end

props = regionprops(regLabels, im_gr, 'BoundingBox', 'Area','Solidity',...
'Orientation','Centroid');

% Filter the region based on the region property
nTotalArea = numel(im_gr);
thr_area = thr_area * nTotalArea;
nReg = size(props,1);
bTests = zeros(nReg, 4);
nTextImg = 0;
rindexs = [];

for ir = 1 : nReg
	% Label regions
    xy = props(ir).Centroid;
	if plotObj
        figure(3)
	    rectangle('Position', props(ir).BoundingBox,'LineWidth', 1, 'EdgeColor', 'y')
        text(xy(1)-5, xy(2), num2str(ir), 'FontSize', 16, 'Color', 'b');
	end

    % Test area
    if(props(ir).Area > thr_area(1) && props(ir).Area < thr_area(2))
        bTests(ir,1) = 1;
    end
    % Test orientation
    if (abs(props(ir).Orientation) < thr_oriAng)
        bTests(ir,2) = 1;
    end
    % Test solidity
    if ((props(ir).Solidity) > thr_solidity)
        bTests(ir,3) = 1;
    end
    % Test the bounding box, won't take too thinny ones
    if (props(ir).BoundingBox(4) > thr_text_hei)
        bTests(ir,4) = 1;
    end
    
    % This region passed the test, and we need to crop it out
    if (sum(bTests(ir,:)) >= 4) 
        if plotObj
            figure(3)
	    	rectangle('Position', props(ir).BoundingBox,'LineWidth', 2, 'EdgeColor', 'g')
            text(xy(1)-5,xy(2),num2str(ir),'FontSize',16,'Color','r');
			% saveas(gcf, sprintf('result/region/%s.jpg', input_img_path(4:end-4)));
        end
        nTextImg = nTextImg + 1;
        rindexs = [rindexs, ir];
    end
end
fprintf('The image has %d candidate regions.\n', nReg)


%% Crop the original image to sub-images
imgs_text = cell(nTextImg, 1);
for ig=1:nTextImg
    ir = rindexs(ig);
    
    up    = uint32(max(1, props(ir).BoundingBox(2)));
    down  = uint32(min(size(im_gr,1), props(ir).BoundingBox(4)+up));
    left  = uint32(max(1, props(ir).BoundingBox(1)));
    right = uint32(min(size(im_gr,2), props(ir).BoundingBox(3)+left));
    
    i_txt = im_gr(up:down, left:right, :);
    imgs_text{ig} = i_txt;
end
fprintf('The image has %d candidate text-blocks.\n', nTextImg)
time.region = toc;


%% Finally, use Tesseract to analyze the text
tic
% Examine each candidate sub-images
fprintf('Extracted text:\n')
txt_id = fopen(output_txt_path, 'w');
for ig = 1 : nTextImg
    % Write the sub-image
    im_path = sprintf('./result/%d.jpg', ig);
    txt_path = sprintf('./result/%d_buff', ig);
    imwrite(imgs_text{ig}, im_path);

    % Use Tesseract by system command
    prefix = 'LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:LD_LIBRARY_PATH; export LD_LIBRARY_PATH;';
    sz_cmd = sprintf('%s %s %s %s', prefix, ocr_name, im_path, txt_path);
    [status, result] = system(sz_cmd);

    % Read and display the OCR results
    sub_txt_id = fopen([txt_path, '.txt'], 'r');
    tline = fgetl(sub_txt_id);
    while ischar(tline)
        fprintf(txt_id, '%s\n', tline);
        fprintf('\t%s\n', tline)
        tline = fgetl(sub_txt_id);
    end
    fclose(sub_txt_id);
    delete([txt_path, '.txt']);
end
fclose(txt_id);

% Adapt text
str = '';
if(nTextImg)
    str = textAdapt(output_txt_path, dataset);
end
time.tesseract = toc;


end

