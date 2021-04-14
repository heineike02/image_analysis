function [x_vec,y_vec] = plotFoundCellsBMH(imdir, imfname, celldata)
%Given the directory of an image, plots the image with the xy locations of
%the cells superimposed
im = imread([imdir,imfname]);

%x_vec = [celldata.Cxloc];
x_vec = [celldata.Cxloc];
x_vec_comb = [x_vec.Combined];
%y_vec = [celldata.Cyloc];
y_vec = [celldata.Cyloc];
y_vec_comb = [y_vec.Combined];


figure
clf 
hold on
imagesc(im)
%plot(y_vec,x_vec,'kx','MarkerSize',5,'LineWidth',2)
plot(y_vec_comb,x_vec_comb,'kx','MarkerSize',5,'LineWidth',2)