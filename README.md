# MathExtended
A library of various math functions created as a DLL. Each function is listed as such:

**[type returned] [function name]([input values]):**

[short descriptions]

## Set Operations

**double[] ArrayToSet(double[] a):**

Calculates set from array a. 

**double[] UnionSet(double[] a, double[] b):**

Returns the union of set a and set b.

**double[] IntersectionSet(double[] a, double[] b):**

Returns the intersection of set a and set b.

**double[] ComplimentSet(double[] s, double[] a):**

Gives the compliment of a given the universal set s.

**double[] SubtractSet(double[] a, double[] b):**

Subracts set b from set a.

**bool IsDisjointSet(double[] a, double[] b):**

Checks if a and b are disjoint sets.

**bool AreEqualSets(double[] a, double[] b)**

Checks if a and b are equal sets.

## Sequencing

**double GeometricSequence(double a, double r, double n):**

Calculates the Geometric Sequence where a is is the starting term, r is the rate of growth, and n is the desired nth position in the sequence.

**double ArithmeticSequence(double a, double n, double d):**

Calculates the Arithmetic Sequence where a is the starting term, n is the desired nth position in the sequence and d is the growth of the sequence.

## Matrix Operations

**double[][] AddMatrix(double[][] x, double[][] y):**

Adds the matrices x and y together. Will return null of the matrices can't be added.

**double[][] SubtractMatrix(double[][] x, double[][] y):**

Subtracts matrix y from matrix x. Will return null if matrices can't be subtracted.

**double[][] MultiplyMatrixByNumber(double n, double[][] m):**

Multiplies the matrix m by the coefficent n. Will return null if matrix size is 0.

**double[][] MultiplyMatrices(double[][] x, double[][] y):**

Multiples the matrices x and y together. Will return null if one of the matrixes is 0 or if multiplication can't be done.

**double? MatrixDetermient(double[][] m):**

Finds the Determient of the matrix m. Will return null if matrix is 0, or matrix is not square.

**double[][] InverseMatrix(double[][] m):**

Finds the inverse matrix of the matrix m. Will return null if matrix is 0 or matrix does not have a determient.

## Equation Solving

**double[][] SolveByMatrix(double[][] m, double[][] a):**

Calculates the solution to an equation using the matrix method. Will return null if the equation can't be solved. Equation Matrix m must be horizontal and answer a matrix must be vertical. Given the equation set x + y + z = 6, 2y + 5z = -4, and 2x + 5y - z = 27, the matrix should be as followed


m = { {1, 1, 1}, {0, 2, 5}, {2, 5, -1} }


a = { {6}, {-4}, {27} }

## Greatest Common Divider and Least Common Multiple

**int GCD(int a, int b):**

Finds the Greatest Common Divider between a and b.

**int LCM(int a, int b):**

Finds the Least Common Multiple between a and b.

## Factoring

**int[] Factors(int x):**

Returns all factors of the number x.

**int[] PrimeFactors(int x):**

Returns list of all prime factors of the number x.

## Point Operations

**double distanceBetweenPoints(double[] a1, double[] a2):**

Calculates the distance between a1 and a2. a1 and a2 should each be an array of (x,y) coordinates.

**double Slope(double[] a1, double[] a2):**

Calculates the slope of the line between the points a1 and a2. a1 and a2 should each be an array of (x,y) coordinates.

**double[] MidPoint(double[] a1, double[] a2):**

Calculates the midpoint of the line between the points a1 and a2. a1 and a2 should each be an array of (x,y) coordinates. Returns as an array of (x, y) coordinates.

**double[] slopeInterceptGivenPoints(double[] a1, double[] a2):**

Calculates the slope and y-intercept given the two points a1 and a2. Returns slope as double[0] and y-intercept as double[1].

## Statistic Operations
**double standardDeviationFull(double[] m):**

Calculates the standard deviation for a set of complete data.

**double standardDeviationSample(double[] m):**

Calculates the standard deviation for a set of sample data.

**double zScore(double x, double m, double d):**

Calculates the z-score of a normal distribution given a value x, the mean of m, and the standard deviation of d.

**double zScoureWithValues(double x, double[] a):**

Calculates the z-score of a normal distribution given the value x and the set of data a.


## Factorial
**int Factorial(int n):**

Calculates n!.

## Permutation and Combinations
**int Permutation(double n, double r):**

Calculates P(n,r).

**int Combination(double n, double r):**

Calculates C(n,r).

## Probability and Percent
**double Probability(double s, double p):**

Calculates the probability given s successes and p attempts.

**double Percent(double? x, double? y, double? p):**

Calculates the missing value for the percent equation x/y = p. x is the whole, y is the fraction, p is the percent.
Pass a null for the missing value. For example: what is 16% of 60 would be passed as Percent(null, 60, 16). If all
three values are filled in, the function automatically returns y / x. if more then one space is null, it will return
a -1.

## Shape Operations

### Triangle

**double areaOfTriangle(double b, double h):**

Calculates the area of a triangle with a base of b and a height of h.


### Rectangle
**double areaOfRectangle(double l, double w):**

Calculates the area of a rectangle with the length of l and the width of w.

**double perimeterOfRectangle(double l, double w):**

Calculates the perimeter of a rectangle with the length of l and the width of w.


### Circle
**double areaOfCircleRadius(double r):**

Calculates the area of a circle given the radius r.

**double areaOfCircleDiameter(double d):**

Calculates the area of a circle given the diameter d.

**double circumferenceOfCircleRadius(double r):**

Calculates the circumference of a circle given the radius r.

**double circumferenceOfCircleDiameter(double d):**

Calculates the circumference of a circle given the diameter d.

**double lengthOfArc(double a, double r):**

Calculates the length of an arc of a circle given the angle a and the radius r.

**double areaOfCircleSector(double a, double r):**

Calculates the area of a sector of a circle given the angle a and the radius r.


### 3d shapes
**double volumeOfBox(double l, double w, double h):**

Calculates the volcume of a box given the length of l, the width of w, and the height of h.

**double surfaceOfBox(double l, double w, double h):**

Calculates the surface area of a box given the length of l, the width of w, and the height of h.

**double volumeOfSphereRadius(double r):**

Calculates the volume of a sphere given the radius r.

**double valumeOfSphereDiameter(double d):**

Calculates the volume of a sphere given the diameter d.


### Other shapes
**double volumeOfCylinder(double r, double h):**

Calculates the volume of a Cylinder given the radius r and the height h.

**double volumeOfCone(double r, double h):**

Calculates the volume of a cone given the radius r and the height h.

**double volumeOfPyramid(double b, double h):**

Calculates the volume of a pyramid given the base b and the height h.

**double areaOfTrapezoid(double b1, double b2, double h):**

Calculates the area of a Trapezoid given the bases b1 and b2 and the height h.

**double areaOfParallelgram(double b, double h):**

Calculates the area of a parallelogram with the base of b and the height of h.

**double sumOfInteriorAnglesOfPolygon(double n):**

Calculates the interior sum of angles of a polygon given the number of angles n.

## Pythagorean Theorem
**double pythMissingAorB(double ab, double c):**

Calculates the missing length of a triangle of the formula a^2 + b^2 = c^2 when a or b is missing.

**double pythMissingC(double a, double b):**

Calculates the missing length of a triangle of the formula a^2 + b^2 = c^2 when c is missing.

## Quadratic Operations
**Object[] quadraticEquation(double a, double b, double c):**

Calculates the result of a quadrati Equation using the quadratic formula given the form ax^2 + bx + c = 0. 

If Object[2] is true, the answers are imaginary and Object[0] is the real portion of the answer and Object[1] is the imaginary part.

If Object[2] is false, the answers are real and Object[0] and Object[1] are each possible answers.

## Fraction Operations
**int[] DecimalToFraction(double x):**

Returns the the decimal x in a fractional version with the numerator being int[0] and the demoninator being int[1];

**double FractionToDecimal(int[] x):**

Returns the fraction x with the numerator being x[0] and denominator being x[1] as a decimal.

**double Discrimiant(double a, double b, double c):**

Calculates the discrimiant of a quadratic equation given the form ax^2 + bx + c = 0.

**double vertexOfParabola(double a, double b, double c):**

Calculates the vertex of a parabola given the form ax^2 + bx + c = 0.

## Averages
**double meanAverage(double[] a):**

Calculates the mean average of the set of data a.

**double modeAverage(double[] a):**

Calculates the mode average of the set of data a.

**double medianAverage(double[] a):**

Calculates the median average of the set of data a.

**double rangeAverage(double[] a):**

Calculates the range of the data set a.

## Distance Rate and Time
**double DistanceRateTime(double? r, double? t, double? d):**

Calculates the missing value for the d = r x t equation. One value should be set to null as the value to be solved. If
all three values are set, the function will automatically return r x t. If more then one value is set to null, the function
will return a -1.
