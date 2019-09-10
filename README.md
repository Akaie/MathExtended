# MathExtended
A library of various math functions created as a DLL. Each function is listed as such:

**[type returned] [function name]([input values]):**

[short descriptions]

## Table of Contents
[Set Operations](#Set Operations)
[Sequencing](#Sequencing)
[Matrix Operations](#Matrix Operations)
[Greatest Common Divider and Least Common Multiple](#Greatest Common Divider and Least Common Multiple)
[Fraction Operations](#Fraction Operations)
[Factoring](#Factoring)
[Point Operations](#Point Operations)
[Statistic Operations](#Statistic Operations)
[Factorial](#Factorial)
[Permutation and Combinations](#Permutation and Combinations)
[Probability and Percent](#Probability and Percent)
[Shape Operations](#Shape Operations)
[Pythagorean Theorem](#Pythagorean Theorem)
[Quadratic Operations](#Quadratic Operations)
[Distance Rate and Time](#Distance Rate and Time)

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

**int TriangularSequence(int n):**

Gets the nth ideration of the triangular Sequence.

**int FibonacciSequence(int n):**

Gets the nth ideration of the Fibonacci Sequence.

**int CatalanSequence(int n):**

Gets the nth ideration of the Catalan Sequence.
## Matrix Operations

**double[][] AddMatrix(double[][] x, double[][] y):**

Adds the matrices x and y together.

**double[][] SubtractMatrix(double[][] x, double[][] y):**

Subtracts matrix y from matrix x.

**double[][] MultiplyMatrixByNumber(double n, double[][] m):**

Multiplies the matrix m by the coefficent n.

**double[][] MultiplyMatrices(double[][] x, double[][] y):**

Multiples the matrices x and y together.

**double? MatrixDetermient(double[][] m):**

Finds the Determient of the matrix m.

**double[][] InverseMatrix(double[][] m):**

Finds the inverse matrix of the matrix m.

**double[][] SolveByMatrix(double[][] m, double[][] a):**

Calculates the solution to an equation using the matrix method. Will return null if the equation can't be solved. Equation Matrix m must be horizontal and answer a matrix must be vertical. Given the equation set x + y + z = 6, 2y + 5z = -4, and 2x + 5y - z = 27, the matrix should be as followed


m = { {1, 1, 1}, {0, 2, 5}, {2, 5, -1} }


a = { {6}, {-4}, {27} }

## Greatest Common Divider and Least Common Multiple

**int GCD(int a, int b):**

Finds the Greatest Common Divider between a and b.

**int LCM(int a, int b):**

Finds the Least Common Multiple between a and b.

## Fraction Operations

**int[] DecimalToFraction(double x):**

Finds the fractional form of the decimal number x.

**int[] SimplifyFraction(int[] x):**

Simplifies the fraction x[0] / x[1].

**double FractionToDecimal(int[] x):**

Returns the decimal value of the fraction x[0] / x[1].

**double[] AddFrazctions(double[] x, double[] y):**

Adds the fractions x[0] / x[1] and y[0] / y[1].

**double[] SubractFractions(double[] x, double[] y):**

Subracts the fraction y[0] / y[1] from x[0] / x[1].

**double[] MultiplyFractions(double[] x, double[] y):**

Multiplies the fractions x[0] / x[1] and y[0] / y[1].

**double[] DivideFractions(double[] x, double[] y):**

Divides the fraction x[0] / x[1] by y[0] / y[1].

## Factoring

**int[] Factors(int x):**

Returns all factors of the number x.

**int[] PrimeFactors(int x):**

Returns list of all prime factors of the number x.

## Point Operations

**double DistanceBetweenPoints(double[] a1, double[] a2):**

Calculates the distance between a1 and a2. a1 and a2 should each be an array of (x,y) coordinates.

**double Slope(double[] a1, double[] a2):**

Calculates the slope of the line between the points a1 and a2. a1 and a2 should each be an array of (x,y) coordinates.

**double[] MidPoint(double[] a1, double[] a2):**

Calculates the midpoint of the line between the points a1 and a2. a1 and a2 should each be an array of (x,y) coordinates. Returns as an array of (x, y) coordinates.

**double[] SlopeInterceptGivenPoints(double[] a1, double[] a2):**

Calculates the slope and y-intercept given the two points a1 and a2. Returns slope as double[0] and y-intercept as double[1].

## Statistic Operations
**double StandardDeviation(double[] m, bool s):**

Calculates the standard deviation for a set of data. s is true if the m is a sample, false if m is a full population.

**double ZScore(double x, double[] a, bool s):**

Calculates the z-score of a normal distribution given a value x and the dataset of a. s is true if the a is a sample, false if a is a full population.

**double Variance(double[] a, bool s):**

Calculates the Variance of a. s is true if the a is a sample, false if a is a full population.

**double CoefficientOfVariation(double[] a, bool s):**

Calculates the Coefficient of the Variation for the dataset a. s is true if the a is a sample, false if a is a full population.

**double StandardError(double[] a, bool s):**

Calculates the Standard Error of the dataset a. s is true if the a is a sample, false if a is a full population.

### Classical Average Operations

**double MeanAverage(double[] a):**

Calculates the mean average of the set of data a.

**double ModeAverage(double[] a):**

Calculates the mode average of the set of data a.

**double MedianAverage(double[] a):**

Calculates the median average of the set of data a.

**double RangeAverage(double[] a):**

Calculates the range of the data set a.

**double MidrangeAverage(double[] a):**

Calculates the midrange of the dataset a.

### Quartile Operations

**double FirstQuartile(double[] a):**

Calculates the first quartile of the dataset a.

**double SecondQuartile(double[] a):**

Calculates the second quartile of the dataset a.

**double ThirdQuartile(double[] a):**

Calculates the third quartile of the dataset a.

**double InterquatileRange(double[] a):**

Calculates the interquatile range of the dataset a.

**double Midhinge(double[] a):**

Calculates the midhinge of the dataset a.

### Other Means

**double WeightedMean(double [] a, double [] w):**

Calculates the weighted mean of the dataset a where w[x] is the weight of a[x]. w and a must be the same length and all w values must be postitive. Additionally, at least one weight must be greater then 0.

**double GeometricMean(double[] a):**

Calculates the geometric mean of the dataset a. All values of a must be greater then 0.

**double WeightedGeometricMean(double[] a, double[] w):**

Calculates the weighted geometric mean for the dataset a. w and a must be the same length and all w values must be postitive. Additionally, at least one weight must be greater then 0.

**double HarmonicMean(double[] a):**

Calculates the harmonic mean of the dataset a. All values of a must be greater then 0.

**double WeightedHarmonicMean(double[] a, double[] w):**

Calulates the weighted harmonic mean of the dataset a. All values of a must be greater then 0. w and a must be the same length. w values cannot be negative. At least one value of w must be greater then 0.

**double GeneralizedMean(double[] a, double p):**

Calculates the generalized mean of the dataset a according to the p-norm p.

**double QuadraticMean(double[] a):**

Calculates the quadratic mean of the dataset a.

**double CubicMean(double[] a):**

Calculates the cubic mean of the dataset a.

**double ModifiedMean(double[] a):**

Calculates the modified mean, discarding the first and last values and calculating the mean of the remaining values.

**double RoundedTruncatedMean(double[] a, int tp):**

Calculates the truncated mean by taking tp% total from the dataset a, or (tp/2)% from each side. Rounds the number of items taken to the nearest whole number. For example, if the number to be taken from each side is 2.7, it rounds to 3 and takes 3 from each side.

**double TruncatedMean(double[] a, int tp):**
Calculates the truncated mean by taking tp% total from the dataset a, or (tp/2)% from each side. Takes reduced weight from the innermost removed items. For example, if the number to be taken from each side is 2.7, it would remove 2 and take 30% of the 3rd values in.

**double InterquartileMean(double[] a):**

Calculates the interquartile mean of the dataset a.

**double WinsorizedMean(double[] a, double tp):**

Calculates the winsorized mean of the dataset at, replaing tp% of the dataset with the upper and lower values, or (tp/2)% of each side.
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

**double AreaOfTriangle(double b, double h):**

Calculates the area of a triangle with a base of b and a height of h.


### Rectangle
**double AreaOfRectangle(double l, double w):**

Calculates the area of a rectangle with the length of l and the width of w.

**double PerimeterOfRectangle(double l, double w):**

Calculates the perimeter of a rectangle with the length of l and the width of w.


### Circle
**double AreaOfCircleRadius(double r):**

Calculates the area of a circle given the radius r.

**double AreaOfCircleDiameter(double d):**

Calculates the area of a circle given the diameter d.

**double CircumferenceOfCircleRadius(double r):**

Calculates the circumference of a circle given the radius r.

**double CircumferenceOfCircleDiameter(double d):**

Calculates the circumference of a circle given the diameter d.

**double LengthOfArc(double a, double r):**

Calculates the length of an arc of a circle given the angle a and the radius r.

**double AreaOfCircleSector(double a, double r):**

Calculates the area of a sector of a circle given the angle a and the radius r.


### 3d shapes
**double VolumeOfBox(double l, double w, double h):**

Calculates the volcume of a box given the length of l, the width of w, and the height of h.

**double SurfaceOfBox(double l, double w, double h):**

Calculates the surface area of a box given the length of l, the width of w, and the height of h.

**double VolumeOfSphereRadius(double r):**

Calculates the volume of a sphere given the radius r.

**double ValumeOfSphereDiameter(double d):**

Calculates the volume of a sphere given the diameter d.


### Other shapes
**double VolumeOfCylinder(double r, double h):**

Calculates the volume of a Cylinder given the radius r and the height h.

**double VolumeOfCone(double r, double h):**

Calculates the volume of a cone given the radius r and the height h.

**double VolumeOfPyramid(double b, double h):**

Calculates the volume of a pyramid given the base b and the height h.

**double AreaOfTrapezoid(double b1, double b2, double h):**

Calculates the area of a Trapezoid given the bases b1 and b2 and the height h.

**double AreaOfParallelgram(double b, double h):**

Calculates the area of a parallelogram with the base of b and the height of h.

**double SumOfInteriorAnglesOfPolygon(double n):**

Calculates the interior sum of angles of a polygon given the number of angles n.

## Pythagorean Theorem
**double PythMissingAorB(double ab, double c):**

Calculates the missing length of a triangle of the formula a^2 + b^2 = c^2 when a or b is missing.

**double PythMissingC(double a, double b):**

Calculates the missing length of a triangle of the formula a^2 + b^2 = c^2 when c is missing.

## Quadratic Operations
**Object[] QuadraticEquation(double a, double b, double c):**

Calculates the result of a quadrati Equation using the quadratic formula given the form ax^2 + bx + c = 0. 

If Object[2] is true, the answers are imaginary and Object[0] is the real portion of the answer and Object[1] is the imaginary part.

If Object[2] is false, the answers are real and Object[0] and Object[1] are each possible answers.

**double Discrimiant(double a, double b, double c):**

Calculates the discrimiant of a quadratic equation given the form ax^2 + bx + c = 0.

**double VertexOfParabola(double a, double b, double c):**

Calculates the vertex of a parabola given the form ax^2 + bx + c = 0.

## Distance Rate and Time
**double DistanceRateTime(double? r, double? t, double? d):**

Calculates the missing value for the d = r x t equation. One value should be set to null as the value to be solved. If
all three values are set, the function will automatically return r x t. If more then one value is set to null, the function
will return a -1.
