function xy_shift = ImageShift_BMH(image1, image2, max_shift_rad)
% Determines how much an image shifted between two timepoints by moving a
% subimage around in a radius and minimizing the difference between the two
% images.
%
%max_shift_rad = maximum shift expected
%image1 = first image
%image2 = second image

%Pick xy shifts to sample

mm = 0; 
for jj = 0:max_shift_rad
    for kk = 0:max_shift_rad
        if sqrt(jj^2+kk^2) <= max_shift_rad
            mm = mm + 1;
            x(mm) = jj;
            y(mm) = kk;
        end
    end
end

x = [x,-x,x,-x];
y = [y,y,-y,-y];

%make base image - original image minus border the size of the max shift
%radius

xx = 0;
yy = 0; 
[Nx,Ny] = size(image1);
im_base = image1((1+max_shift_rad-xx):(Nx-max_shift_rad-xx), (1+max_shift_rad-yy):(Ny-max_shift_rad-yy));

%shift images
for jj = 1:length(x);
   xx = x(jj);
   yy = y(jj);
   im_test = image2((1+max_shift_rad-xx):(Nx-max_shift_rad-xx),(1+max_shift_rad-yy):(Ny-max_shift_rad-yy));
   dist(jj) = norm(double(im_test-im_base));
end

[~,min_ind] = min(dist);
xy_shift = [x(min_ind), y(min_ind)];

end
