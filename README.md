# MathExtended
A library of various math functions created as a DLL. Each function is listed as such:

[type returned] [function name]([input values]): [short descriptions]

int GCD(int a, int b):
Finds the Greatest Common Divider between a and b.

int LCM(int a, int b):
Finds the Least Common Multiple between a and b.

int[] Factors(int x):
Returns all factors of the number x.

int[] PrimeFactors(int x):
Returns list of all prime factors of the number x.

double distanceBetweenPoints(double[] a1, double[] a2):
Calculates the distance between a1 and a2. a1 and a2 should each be an array of (x,y) coordinates.

double Slope(double[] a1, double[] a2):
Calculates the slope of the line between the points a1 and a2. a1 and a2 should each be an array of (x,y) coordinates.

double[] MidPoint(double[] a1, double[] a2):
Calculates the midpoint of the line between the points a1 and a2. a1 and a2 should each be an array of (x,y) coordinates. Returns as an array of (x, y) coordinates.

double[] slopeInterceptGivenPoints(double[] a1, double[] a2):
Calculates the slope and y-intercept given the two points a1 and a2. Returns slope as double[0] and y-intercept as double[1].

double standardDeviationFull(double[] m):
Calculates the standard deviation for a set of complete data.

double standardDeviationSample(double[] m):
Calculates the standard deviation for a set of sample data.

double zScore(double x, double m, double d):
Calculates the z-score of a normal distribution given a value x, the mean of m, and the standard deviation of d.

double zScoureWithValues(double x, double[] a):
Calculates the z-score of a normal distribution given the value x and the set of data a.

double Permutation(double n, double r):
Calculates P(n,r).

double Combination(double n, double r):
Calculates C(n,r).

double Probability(double s, double p):
Calculates the probability given s successes and p attempts.

double Percent(double? x, double? y, double? p):
Calculates the missing value for the percent equation x/y = p. x is the fraction, y is the whole, p is the percent.
Pass a null for the missing value. For example: what is 16% of 60 would be passed as Percent(null, 60, 16).

double areaOfTriangle(double b, double h):
Calculates the area of a triangle with a base of b and a height of h.

double areaOfRectangle(double l, double w):
Calculates the area of a rectangle with the length of l and the width of w.

double perimeterOfRectangle(double l, double w):
Calculates the perimeter of a rectangle with the length of l and the width of w.

double areaOfCircleRadius(double r):
Calculates the area of a circle given the radius r.

double areaOfCircleDiameter(double d):
Calculates the area of a circle given the diameter d.

double circumferenceOfCircleRadius(double r):
Calculates the circumference of a circle given the radius r.

double circumferenceOfCircleDiameter(double d):
Calculates the circumference of a circle given the diameter d.

double lengthOfArc(double a, double r):
Calculates the length of an arc of a circle given the angle a and the radius r.

double areaOfCircleSector(double a, double r):
Calculates the area of a sector of a circle given the angle a and the radius r.

double volumeOfBox(double l, double w, double h):
Calculates the volcume of a box given the length of l, the width of w, and the height of h.

double surfaceOfBox(double l, double w, double h):
Calculates the surface area of a box given the length of l, the width of w, and the height of h.

double volumeOfSphereRadius(double r):
Calculates the volume of a sphere given the radius r.

double valumeOfSphereDiameter(double d):
Calculates the volume of a sphere given the diameter d.

double volumeOfCylinder(double r, double h):
Calculates the volume of a Cylinder given the radius r and the height h.

double volumeOfCone(double r, double h):
Calculates the volume of a cone given the radius r and the height h.

double volumeOfPyramid(double b, double h):
Calculates the volume of a pyramid given the base b and the height h.

double areaOfTrapezoid(double b1, double b2, double h):
Calculates the area of a Trapezoid given the bases b1 and b2 and the height h.

double areaOfParallelgram(double b, double h):
Calculates the area of a parallelogram with the base of b and the height of h.

double sumOfInteriorAnglesOfPolygon(double n):
Calculates the interior sum of angles of a polygon given the number of angles n.

double pythMissingAorB(double ab, double c):
Calculates the missing length of a triangle of the formula a^2 + b^2 = c^2 when a or b is missing.

double pythMissingC(double a, double b):
Calculates the missing length of a triangle of the formula a^2 + b^2 = c^2 when c is missing.

Object[] quadraticEquation(double a, double b, double c):
Calculates the result of a quadrati Equation using the quadratic formula given the form ax^2 + bx + c = 0. 
If Object[2] is true, the answers are imaginary and Object[0] is the real portion of the answer and Object[1] is the imaginary part.
If Object[2] is false, the answers are real and Object[0] and Object[1] are each possible answers.

double Discrimiant(double a, double b, double c):
Calculates the discrimiant of a quadratic equation given the form ax^2 + bx + c = 0.

double vertexOfParabola(double a, double b, double c):
Calculates the vertex of a parabola given the form ax^2 + bx + c = 0.

double meanAverage(double[] a):
Calculates the mean average of the set of data a.

double modeAverage(double[] a):
Calculates the mode average of the set of data a.

double medianAverage(double[] a):
Calculates the median average of the set of data a.

double rangeAverage(double[] a):
Calculates the range of the data set a.
