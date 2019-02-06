## SPHERETRI

is a set of high-performance vectorized MATLAB functions for building
a triangulation of a unit sphere based on recursive partitioning of each
of Icosahedron faces into 4 triangles with vertices in  the middles of 
original face edgeMidMat. `spheretri` function takes a number of  requested 
points as an input. An exact number of points in the returned
triangulation may not match this number, the function will choose a depth
of Icosahedron partitioning minimally sufficient to provide the requested
number of points. Alternatively you can use `spheretribydepth` function 
that expects a partitioning depth as an input.

##### Examples:
 
    [vMat, fMat] = spheretri(500);
    patch('Vertices',vMat,'Faces',fMat,'FaceColor','g','EdgeColor','k');

![Result: ](Pictures/spheretriresult.png)

    [vMat, fMat] = spheretribydepth(3);
    patch('Vertices',vMat,'Faces',fMat,'FaceColor','g','EdgeColor','k');

![Result: ](Pictures/spheretribydepthresult.png)

##### Copyright:
 Peter Gagarinov, PhD, <br>
 Moscow State University, <br>
 &ensp;&ensp;&ensp; Faculty of Computational Mathematics and Computer Science,<br> &ensp;&ensp;&ensp; System Analysis Department 2011-2017
