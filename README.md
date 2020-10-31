# voronoi

Here are some utilities for creating and displaying Voronoi Cells.

## Demo

Live demo utilizing R Shiny: https://olimay.shinyapps.io/MidtermVoronoi/

It displays the cell generation points in red and the cell boundaries in black.

_Note: The source code for the demo above reuses code, with many thanks, developed by
course staff for the Fall 2020 Harvard Extension School course,
*MATH E-151: Classic Mathematics with a Modern User Interface* that is not part of this
repository._

## Voronoi Diagram

### About Voronoi Diagrams

A Voronoi Diagram is a partitioning of a plane into separate cells containing a
particular point called the &quot;site&quot; of the cell.  The boundaries of
each cell are constructed so that all points in a particular cell are closer to
that cell's &quot;site&quot; as they are to any other cell's site. These cells
are called *Voronoi Cells* [5].

The Voronoi Diagram is named after Russian mathematician [Georgy Feodosevich
Voronoi](https://en.wikipedia.org/wiki/Georgy_Voronoy).

I am extremely new to this topic, having stumbled on the term while looking for
information on regular plane tilings.

According to Wikipedia [2]:

> They are used in many areas of science, such as the analysis of spatially
> distributed data, having become an important topic in geophysics, meteorology,
> condensed matter physics, and Lie groups. These tessellations are widely used
> in many areas of computer graphics, from architecture to film making and video
> games.

### About this repository

These source files include functions for generating Voronoi Cells for a set of
points, functions for generating a random set of points, and examples of tests.

The R code was originally developed as part of a midterm/final project for
the Fall 2020 Harvard Extension School course, E-151, *Classic Mathematics with
a Modern User Interface*.

It can generate two types of Voronoi Cells for two distance metrics:
Euclidean distance (equivalently known as L2) and Manhattan distance (known as
L1).

The algorithm used to generate the cells was adapted from [this
assignment](https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html)
from the University of Washington [1].

I adapted C++ code for computing spatial relationships of points and line
segments from [Dan Sunday's page on computational
geometry](http://geomalgorithms.com/a05-_intersect-1.html) [4].

I figured out how to generalize certain concepts used in the algorithm,
particularly that of "perpendicular bisector", to Taxicab Geometry/Manhattan
Distance with help from [Aubrey Kemp's 2018 dissertation at Georgia State
University](https://scholarworks.gsu.edu/math_diss/58) [3].

### References

1. Anderson, Ruth. "Algorithm for generation of Voronoi diagrams." CSE326: Data Structures.  https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html (accessed Oct 11, 2020).
2. "Georgy Voronoy." Wikipedia.  https://en.wikipedia.org/wiki/Georgy_Voronoy (accessed Oct 14, 2020).
3. Kemp, Aubrey, "Generalizing and Transferring Mathematical Definitions from Euclidean to Taxicab Geometry." Dissertation, Georgia State University, 2018.  https://scholarworks.gsu.edu/math_diss/58
4. Sunday, Dan. "Line and Plane Intersection (2D &amp; 3D)." Geometry Algorithms.  http://geomalgorithms.com/a05-_intersect-1.html (accessed Oct 14, 2020).
5. "Voronoi Diagram." Wikipedia.  https://en.wikipedia.org/wiki/Voronoi_diagram (accessed Oct 14, 2020).
