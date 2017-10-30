function RecognizeText(input_img_path, output_img_path)
  
%% 
set(0,'DefaultTextFontSize',16);
set(0,'DefaultAxesFontSize',16);

%%
bVis = 0;

%% Read in each image
nTile = 4; climit=0.01;
thr_canny = 0.10;

%for i=1%:nimg
  % load the img
  img_3chB = imread(input_img_path);
  imgsz_orig = size(img_3chB);
  img_3ch = imresize(img_3chB,1.0,'bicubic');
  
  % do binarization
  img_gr = im2double(rgb2gray(img_3ch));
  % And do local adaptive thresholding
  %img_eq = adapthisteq(img_gr,'NumTiles',[nTile nTile],'ClipLimit',climit);
  %imshow([img_gr,img_eq]);
  
  % Then do the thresholding
  img = img_gr;
  %img_lv = graythresh(img); img_bw = im2bw(img,img_lv*1.2);
	%imshow(img_bw);
  
  thr_canny = 0.1; sig_canny=8;
  % performe edge detection on img
  edge_bw = edge(img,'canny',thr_canny,sig_canny);
  edge_bw2 = imdilate(edge_bw,ones(5,5));
  
  if (bVis) 
    figure; imshow(img_3ch); title(['size of image:' num2str(size(img_3ch))]);
    figure; imshow([edge_bw, edge_bw2]); grid on; shg;
  end
  
  % Do geometry correction, the hard thing is to identify the four corner points.
  % Try hough transform, detect strong lines near horizontal or vertical
  thetas_hz = -10:0.5:10;
  thetas_vt = thetas_hz+90;
  % [H,T,R] = hough(edge_bw2, 'Theta', thetas_hz);
  [H,T,R] = hough(edge_bw2, 'ThetaResolution',thetas_hz(2)-thetas_hz(1),'RhoResolution',4);
  
  %% [H,T,R] = hough(edge_bw2, 'ThetaResolution',thetas_hz(2)-thetas_hz(1),'RhoResolution',2);
  if (bVis)
    figure; imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,'InitialMagnification','fit'); 
    title('Hough Transform'); xlabel('\theta'), ylabel('\rho'); axis on; axis normal; colormap(hot); colorbar;
  end
  %% Detect horizontal lines and veritcal lines respectively
  b_found_lines = [false false; false false]; %hz_1,hz_2,vt_1,vt_2
  thetas_mm = [0 0; 0 0];
  rhos_mm = [0 0; 0 0]; %hz_1,hz_2,vt_1,vt_2
  
  hz_ang = 15; vt_ang = 15;
  t_range_vt = T>-vt_ang & T<vt_ang;
  t_range_hz = T>(90-hz_ang) | T<(hz_ang-90);
  wid = size(img,2); ht = size(img,1);
  r_range_vt = R>=-0.3*wid & R<=1.3*wid;
  r_range_hz = R>=-ht & R<=ht;
  %
  H1 = []; H2 = []; T1=[]; T2=[]; R1=[]; R2=[];
  H1 = H(r_range_hz,t_range_hz); H2 = H(r_range_vt,t_range_vt);
  T1 = T(:,t_range_hz); T2 = T(:,t_range_vt);
  R1 = R(:,r_range_hz); R2 = R(:,r_range_vt);
  
  line_quad = zeros(8,2); d_rhos = zeros(2,1);
  %
  for j=1:2  %j=1 is horizontal line, and j=2 is vertical line
    Ht = []; Tt = [];
    if j==1
      Ht = H1; Tt = T1; Rt = R1;
    else
      Ht = H2; Tt = T2; Rt = R2;
    end
    
    P  = houghpeaks(Ht,16,'threshold',ceil(0.3*max(Ht(:))));
    x = Tt(P(:,2)); y = Rt(P(:,1));

    % plot the hough domain
    %figure; imshow(imadjust(mat2gray(Ht)),'XData',Tt,'YData',Rt,'InitialMagnification','fit');
    %title('Hough Transform'); xlabel('\theta'), ylabel('\rho'); axis on; axis normal; colorbar; hold on; plot(x,y,'bs');
    % Find lines and plot them
    lines = houghlines(edge_bw2,Tt,Rt,P,'FillGap',5,'MinLength',1/5*size(img,3-j));

    % 
    max_len = 0; xy_long = []; k_long = -1; rho_long = 0;
    if (bVis)
      figure, imshow(edge_bw2), hold on;
    end
    max_len = 0; rho_max = 0; rho_min = 0;
    
    if isempty(fieldnames(lines))
      %No line detected
      b_found_lines(j,:) = false;
    else 
      for k = 1:length(lines)
         xy = [lines(k).point1; lines(k).point2];
         
         if (bVis)
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
            % Plot beginnings and ends of lines
            plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
            plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
         end
         
         % Determine the endpoints of the longest line segment
         len = norm(lines(k).point1 - lines(k).point2);
         if ( len > max_len)
            max_len = len;
            xy_long = xy;
            k_long = k;
         end
      end
      % highlight the longest line
      if (bVis)
        plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
      end

      % now identify the four corner points, just find the one with min and max rho.
      nline = length(lines);
      rhos = zeros(nline,1); thetas = zeros(nline,1);
      center_img = [size(img,2) size(img,1)] / 2;

      if j==1 %horizontal line, rho can be pos or neg
        for k=1:nline
          rho = lines(k).rho; the = deg2rad(lines(k).theta); 
          %Need to calculate Rhos from the center of the image rather than the top-left corner
          thetas(k) = the;
          rhos(k) = rho - (center_img(1)*cos(the)+center_img(2)*sin(the));
          if (sin(the)<0)
            rhos(k) = -rhos(k);
            thetas(k) = the+pi;
          end
          if (k==k_long) rho_long = rhos(k); end
        end
      else %vertical line, rho is guaranted to be positive
        for k=1:nline
          rho = lines(k).rho; the = deg2rad(lines(k).theta);
          %Need to calculate Rhos from the center of the image rather than the top-left corner
          thetas(k) = the;          
          rhos(k) = rho - (center_img(1)*cos(the)+center_img(2)*sin(the));
          if (k==k_long); rho_long = rhos(k); end
        end
      end

      thr_bnd_ratio_max = 0.95; thr_bnd_ratio_min = 0;
      %don't take rhos that are too close to boundary
      [rhos_asc,ind_asc] = sort(rhos);
      hf_dim = size(img,j)/2;
      for k=1:length(ind_asc)
        if (rhos_asc(k)>-thr_bnd_ratio_max*hf_dim && rhos_asc(k)<-thr_bnd_ratio_min*hf_dim)
          ind_min = ind_asc(k); rho_min = rhos_asc(k); b_found_lines(j,1) = true; thetas_mm(j,1) = thetas(ind_min);
          break;
        end
      end
      if (b_found_lines(j,1)==false) %if still have the edge points, then take the ones on the ver edge
        for k = length(ind_asc):-1:1       
          if rhos_asc(k)<-thr_bnd_ratio_max*hf_dim
            ind_min = ind_asc(k); rho_min = rhos_asc(k); b_found_lines(j,1) = true; thetas_mm(j,1) = thetas(ind_min);
            break;
          end
        end
      end
      if (b_found_lines(j,1))
        line_quad(1+4*(j-1),:) = lines(ind_min).point1;
        line_quad(2+4*(j-1),:) = lines(ind_min).point2;
      end
      
      
      for k=length(ind_asc):-1:1
        if (rhos_asc(k)<thr_bnd_ratio_max*hf_dim && rhos_asc(k)>thr_bnd_ratio_min*hf_dim)
          ind_max = ind_asc(k); rho_max = rhos_asc(k); b_found_lines(j,2) = true; thetas_mm(j,2) = thetas(ind_max);
          break;
        end
      end
      if (b_found_lines(j,2)==false) %if still have the edge points, then take the ones on the ver edge
        for k = 1:length(ind_asc)       
          if rhos_asc(k)>thr_bnd_ratio_max*hf_dim
            ind_max = ind_asc(k); rho_max = rhos_asc(k); b_found_lines(j,2) = true; thetas_mm(j,2) = thetas(ind_max);
            break;
          end
        end
      end      
      if (b_found_lines(j,2) == true)
        line_quad(3+4*(j-1),:) = lines(ind_max).point1;
        line_quad(4+4*(j-1),:) = lines(ind_max).point2;
      end
      d_rhos(j) = (rho_max-rho_min);
      rhos_mm(j,:) = [rho_min rho_max];
    end %if lines struture is not empty.
  end %j=1:2
  
  %% figure out if there are missing edges
  bh1 = b_found_lines(1,1); bh2 = b_found_lines(1,2);
  bv1 = b_found_lines(2,1); bv2 = b_found_lines(2,2); 
  dist = 200;
  nEdgesFound = sum(b_found_lines(:));
  if (nEdgesFound<4) %need to fill the edges
    %Check horizontal line first
    if (bh1+bh2<2)
      hf_dim = ht/2;
      if ((~bh1&&bh2)) %if only one edge is found
        rhos_mm(1,1) = min(-0.7*hf_dim,-rhos_mm(1,2));
        thetas_mm(1,1) = thetas_mm(1,2);
      elseif (bh1&&~bh2) 
        rhos_mm(1,2) = max(0.7*hf_dim,-rhos_mm(1,1));        
        thetas_mm(1,2) = thetas_mm(1,1);
      elseif (~bh1&&~bh2)
        rhos_mm(1,1) = -hf_dim*thr_bnd_ratio_max;
        rhos_mm(1,2) = -rhos_mm(1,1);
        thetas_mm(1,1) = pi; thetas_mm(1,2) = thetas_mm(1,1);
      end
      d_rhos(1) = rhos_mm(1,2) - rhos_mm(1,1);
      [pt1 pt2] = line_from_rho_theta(rhos_mm(1,1),thetas_mm(1,1),dist);
      [pt3 pt4] = line_from_rho_theta(rhos_mm(1,2),thetas_mm(1,2),dist);
      line_quad(1,:) = pt1+center_img; line_quad(3,:) = pt3+center_img;
      line_quad(2,:) = pt2+center_img; line_quad(4,:) = pt4+center_img;   
    end
    
    %Now check the vertical edges
    if (bv1+bv2<2)
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
      [pt1 pt2] = line_from_rho_theta(rhos_mm(2,1),thetas_mm(2,1),dist);
      [pt3 pt4] = line_from_rho_theta(rhos_mm(2,2),thetas_mm(2,2),dist);
      line_quad(5,:) = pt1+center_img; line_quad(7,:) = pt3+center_img;
      line_quad(6,:) = pt2+center_img; line_quad(8,:) = pt4+center_img;             
    end
  end
  % Find the intersection of the four lines
  ln_hz_1 = [line_quad(1:2,:)];
  ln_hz_2 = [line_quad(3:4,:)];
  ln_vt_1 = [line_quad(5:6,:)];
  ln_vt_2 = [line_quad(7:8,:)];
  
  pt1 = linlinintersect([ln_hz_1; ln_vt_1]);
  pt2 = linlinintersect([ln_hz_1; ln_vt_2]);
  pt3 = linlinintersect([ln_hz_2; ln_vt_2]);
  pt4 = linlinintersect([ln_hz_2; ln_vt_1]);
  
  pts = [pt1; pt2; pt3; pt4];
  if (bVis)
    figure; imshow(img); hold on; scatter(pts(:,1),pts(:,2),10*(1:2:7));
    plot(ln_hz_1(:,1),ln_hz_1(:,2),'LineWidth',1.5,'Color','green');
    plot(ln_hz_2(:,1),ln_hz_2(:,2),'LineWidth',1.5,'Color','green');
    plot(ln_vt_1(:,1),ln_vt_1(:,2),'LineWidth',1.5,'Color','green');
    plot(ln_vt_2(:,1),ln_vt_2(:,2),'LineWidth',1.5,'Color','green');
  end
  
  
  %% Then do the homographic trasform, but keep the aspect ratio
  pt_in = [pt1; pt2; pt3; pt4;];
  pt_out = [0 0; d_rhos(2) 0; d_rhos(2) d_rhos(1); 0 d_rhos(1)];
  Hom = estimate_homography_simp(pt_in',pt_out');
  tf = maketform('projective',Hom');
  img_rect = imtransform(img,tf,'bilinear','XData',[1 d_rhos(2)], 'YData', [1 d_rhos(1)]);
  
  if (bVis)
    figure;imshow(img_rect);
  end
  
  % Detect text blocks using some heruistics %%, first is the pixel value variance
  imgr = img_rect;
  wnd_sz  = 25; 
  boxflt  = fspecial('average',[wnd_sz wnd_sz]);
  img_ex  = imfilter(imgr,boxflt,'replicate');
  img_ex2 = imfilter(imgr.*imgr,boxflt,'replicate');
  img_var = img_ex2 - img_ex.*img_ex;
  thr_imgvar = prctile(img_var(:),90)*0.08; %avoid extream values
  % Use local pixel variance to determine the text region
  
  imgthr_var = img_var; 
  imgthr_var(img_var<thr_imgvar) = 0; %imgthr_var(img_var>thr_imgvar) = max(imgthr_var(:));
  %figure, imshow([img_var imgthr_var],[]); title(['thres\_var=' num2str(thr_imgvar)]); colorbar;
  
  % Also compute the edge density, real text should have a edge density not too big nor small.
  edge_rect = edge(img_rect,'canny',thr_canny);
  edge_den = imfilter(double(edge_rect),boxflt);
  
  edge_den_thr = edge_den; edge_den_thr(edge_den<0.1 | edge_den>0.5) = 0;
  if (bVis)
    figure; imshow([imgthr_var edge_den_thr],[]); colorbar; title('left: brightness variance thr; right: edge density thr');
  end
  
  % Label regions and calculate region property.
  msk_region = edge_den_thr>0 & imgthr_var>0;
  if (bVis)
    figure; imshow([imgthr_var>0, edge_den_thr>0, msk_region],[0 1]); colorbar; 
  end
  msk_region = imdilate(msk_region,ones(21,31));
  regLabels = bwlabel(msk_region,8);
  props = regionprops(regLabels,imgr,'BoundingBox','MajorAxisLength','MinorAxisLength',...
  'Area','Solidity','MeanIntensity','Orientation','Centroid');
  %imshow(regLabels,[]); colormap(jet); colorbar;
  % filter the region based on the region property
  nTotalArea = prod(size(img_rect));
  nReg = size(props,1);
  nht = size(img_rect,1); nwid = size(img_rect,2);
  % Assume img_rect is of 600x900, and the text should be bigger than 30(ht)*60(wid) gives the ratio, should not be too big
  % to contain the key information
  thr_area = [1/20*1/15 1/2]*nTotalArea; thr_oriAng = 10; thr_solidity = 0.65;
  bTests = logical(zeros(nReg,7)); thr_text_ht = 20;
  nTextImg = 0; rindexs = [];
  if (bVis)
    figure; imshow(regLabels,[]); colormap(jet); colorbar; hold on;
  end
  for ir = 1:nReg
    xy = props(ir).Centroid;
    if (bVis)
      text(xy(1)-5,xy(2),num2str(ir),'FontSize',16);
    end
    if (props(ir).Area > thr_area(1) && props(ir).Area < thr_area(2)) %area should be proper
      bTests(ir,1) = 1;
    end
    if (abs(props(ir).Orientation) < thr_oriAng) %should be horizontal
      bTests(ir,2) = 1;
    end
    if ((props(ir).Solidity) > thr_solidity) %should be horizontal
      bTests(ir,3) = 1;
    end
    %Check the bounding box, won't take too thinny ones
    if (props(ir).BoundingBox(4) > thr_text_ht )
      bTests(ir,4) = 1;
    end
    if (sum(bTests(ir,1:4)) >= 4) %this region passed the test, and we need to crop it out
      xy = props(ir).Centroid;
      if (bVis)
        text(xy(1)-5,xy(2),num2str(ir),'FontSize',16,'Color','y');
      end
      nTextImg = nTextImg + 1;
      rindexs = [rindexs, ir];
    end
  end
  
  %%
  imgs_text = cell(nTextImg,1);
  imgs_bw_text = cell(nTextImg,1);
  i_txt = []; ibw_txt = [];
  
  if (bVis)
    figure; ha = tight_subplot(ceil(nTextImg/2),2,[.01 .03],[.1 .01],[.01 .01]);
  end
  for ig=1:nTextImg %crop this part of image out
    ir = rindexs(ig);
    x1 = uint32(max(1,props(ir).BoundingBox(1)));
    y1 = uint32(max(1,props(ir).BoundingBox(2)));
    x2 = uint32(min(nwid,props(ir).BoundingBox(3)+x1));
    y2 = uint32(min(nht,props(ir).BoundingBox(4)+y1));
    i_txt = img_rect(y1:y2,x1:x2);
    % Binarization
    gamma = 1.0;
    i_txt = i_txt.^gamma;
    grThr = graythresh(i_txt);
    ibw_txt = im2bw(i_txt,grThr*1.0);
    
    imgs_text{ig} = i_txt; imgs_bw_text{ig} = imerode(ibw_txt,ones(1,1));
    %subplot (double(uint32(nTextImg/2)),2,ig); imshow([i_txt imgs_bw_text{ig}]); %title(num2str(ig));
    if (bVis)
      axes(ha(ig)); imshow([i_txt imgs_bw_text{ig}]);
    end
  end
  if (bVis)
    set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','');
  end
  %% Then use Tesseract to analyze the text %%
  %write the image to file and then do the OCR
  fn_txtOCRfinal = sprintf('%s.txt',output_img_path);
  fwid = fopen([fn_txtOCRfinal '.txt'],'w');
  for ig=1:nTextImg
    fn_txtImg = sprintf('./result/%d.jpg', ig);
    fn_txtOCR = sprintf('./result/%d', ig);
    imwrite(imgs_text{ig}, fn_txtImg);
    tbin = 'tesseract';
    prefix = 'LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:LD_LIBRARY_PATH; export LD_LIBRARY_PATH;';
    sz_cmd = sprintf('%s %s %s %s', prefix, tbin, fn_txtImg, fn_txtOCR);
    [stat, res] = system(sz_cmd);
    %Read and display the OCR results
    fid = fopen([fn_txtOCR '.txt']);
    tline = fgetl(fid);
    while ischar(tline)
      fprintf(fwid,'%s\n',tline);
      disp(tline)
      tline = fgetl(fid);
    end
    fclose(fid);
  end
  fclose(fwid);
end
