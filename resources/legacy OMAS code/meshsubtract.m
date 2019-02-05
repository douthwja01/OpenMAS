close all; clear all;

figure(1)
% [X,Y] = meshgrid(-2:2,-2:2);
% Z = ones(length(X),length(Y));
% % change the next line for more complex geometries
% Z((X.^2+Y.^2)<=1) = NaN; % for your example
% X(isnan(Z)) = NaN;
% Y(isnan(Z)) = NaN;
% plot(X,Y,'LineStyle','none','marker','*')

x=[2.5 4 6 18 9]; 
y=[12 3 7.5 1 10]; 
z=[3 15 16 8 11.5];
tri = delaunay(x, y);
trimesh(tri, x, y, z);