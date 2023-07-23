Simple implementation of segmented least squares algorithm.
\hline
There are two programs for solving the problem:
One has user messages and instructions for input and comments along the code and the other has none.
Input in the latter is done via the "input.txt" file in the following format:
-------
numberOfPoints(int)
point1x(double)	point1y(double)
.
.
.
pointnx(double) pointny(double)
-------
After the execution of the latter (directInput), solutions will be printed in stdout and also stored in "output1.txt" and "output2.txt".
These two files, along with "input.txt", are then used in "plot.ipynb" to generate plots for better understanding.

Three test examples are given, each with the same 9 points but different values of C.
The first test example, for which C=50 demonstrates only one line since lines are too expensive.
The second one, with C=1 demonstrates three lines; this one would correspond to the most realistic one w.r.t. the square errors and the scale of the range of given points.
The third one, with C=0 demonstrates plotting line through each two points since they do not cost and in this way the total cost of lines and errors is 0 (same would hold for similarly small positive costs of adding lines).

Note that, when the line is set to go through only one point, it is not plotted in the picture and its coefficients are given as a=b=0; this implies that any line through that point will suffice and it is implied in plots.

Regarding the time complexity, it is somewhat obviously quadratic in n, this is stated in the PowerPoint file and will be elaborated on during the presentation.

Finally, a simple program for calculating the square errors for arbitrary sets of points is there for convenience should one wish to check the correctness of square errors and linear fit coefficients.
