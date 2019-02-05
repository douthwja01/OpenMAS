%% Barycentric Coordinates
%  To work out if a point lie within a triangle
function col = Barycentric_Coordinates(p1,p2,p3,p)

    alpha = ((p2(:,2) - p3(:,2))*(p(:,1) - p3(:,1)) + (p3(:,1) - p2(:,1))*(p(:,2) - p3(:,2))) / ((p2(:,2) - p3(:,2))*(p1(:,1) - p3(:,1)) + (p3(:,1) - p2(:,1))*(p1(:,2) - p3(:,2)));
    beta  = ((p3(:,2) - p1(:,2))*(p(:,1) - p3(:,1)) + (p1(:,1) - p3(:,1))*(p(:,2) - p3(:,2))) / ((p2(:,2) - p3(:,2))*(p1(:,1) - p3(:,1)) + (p3(:,1) - p2(:,1))*(p1(:,2) - p3(:,2)));
    gamma = 1 - alpha - beta;

    if ((alpha > 0) && (beta > 0) && (gamma > 0))
        col = p;
    end
end