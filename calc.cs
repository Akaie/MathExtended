using System;
using System.Collections.Generic;

namespace MathExpanded
{
    public static class Calc
    {   //-----------------------------------------Set Operations--------------------------------------
        public static double[] ArrayToSet(double[] a)
        {
            List<double> l = new List<double>();
            for(int i = 0; i < a.Length; i++)
            {
                if (l.Contains(a[i]))
                    continue;
                l.Add(a[i]);
            }
            return l.ToArray();
        }
        public static double[] UnionSet(double[] a, double[] b)
        {
            List<double> l = new List<double>(a);
            for(int i = 0; i < b.Length; i++)
            {
                if (l.Contains(b[i]))
                    continue;
                l.Add(b[i]);
            }
            return l.ToArray();
        }

        public static double[] IntersectionSet(double[] a, double[] b)
        {
            List<double> bl = new List<double>(b);
            List<double> l = new List<double>();
            for (int i = 0; i < a.Length; i++) {
                if (bl.Contains(a[i]))
                    l.Add(a[i]);
            }
            return l.ToArray();
        }

        public static double[] ComplimentSet(double[] s, double[] a)
        {
            List<double> al = new List<double>(a);
            List<double> l = new List<double>();
            for (int i=0; i < s.Length; i++)
            {
                if(!al.Contains(s[i])) {
                    l.Add(s[i]);
                }
            }
            return l.ToArray();
        }

        public static double[] SubtractSet(double[] a, double[] b)
        {
            List<double> l = new List<double>(a);
            for(int i = 0; i<b.Length; i++)
            {
                if(l.Contains(b[i]))
                {
                    l.Remove(b[i]);
                }
            }
            return l.ToArray();
        }

        public static bool IsDisjointSet(double[] a, double[] b)
        {
            List<double> l = new List<double>(a);
            for(int i = 0; i < b.Length; i++)
            {
                if (l.Contains(b[i]))
                    return false;
            }
            return true;
        }

        public static bool AreEqualSets(double[] a, double[] b)
        {
            if (a.Length != b.Length)
                return false;
            List<double> l = new List<double>(a);
            for (int i = 0; i < b.Length; i++)
            {
                if (!l.Contains(b[i]))
                    return false;
            }
            return true;
        }
        //-------------------------------------------Geometric and arithmetic----------------------------
        public static double GeometricSequence(double a, double r, double n)
        {
            return a * (1 - Math.Pow(r, n)) / (1 - r);
        }
        public static double ArithmeticSequence(double a, double n, double d)
        {
            return a + ((n - 1) * d);
        }
        //------------------------------------------Matrix Operations------------------------------------
        public static double[][] AddMatrix(double[][] x, double[][] y)
        {
            if (x.Length != y.Length)
            {
                throw new FormatException("Matrices must be the same size.");
            }
            if (x.Length == 0 || y.Length == 0 || x[0].Length == 0 || y[0].Length == 0)
            {
                throw new FormatException("Matrices cannot be 0.");
            }
            for(int i=0; i<x.Length-1; i++)
            {
                if (x[i].Length != x[i + 1].Length)
                    throw new FormatException("Matrices rows must be the same length.");
                if (y[i].Length != y[i + 1].Length)
                    throw new FormatException("Matrices rows must be the same length.");
                if (x[i].Length != y[i].Length)
                    throw new FormatException("Matrices must be the same size.");
            }
            double[][] ans = new double[x.Length][];
            for (int i = 0; i < x.Length; i++)
            {
                ans[i] = new double[i];
                for (int ii = 0; ii < x[i].Length; ii++)
                {
                    ans[i][ii] = x[i][ii] + y[i][ii];
                }
            }
            return ans;
        }

        public static double[][] SubtractMatrix(double[][] x, double[][] y)
        {
            if (x.Length != y.Length)
            {
                throw new FormatException("Matrices must be the same size.");
            }
            if(x.Length == 0 || y.Length == 0 || x[0].Length == 0 || y[0].Length==0)
            {
                throw new FormatException("Matrices cannot be 0.");
            }
            for (int i = 0; i < x.Length - 1; i++)
            {
                if (x[i].Length != x[i + 1].Length)
                    throw new FormatException("Matrices must be the same length.");
                if (y[i].Length != y[i + 1].Length)
                    throw new FormatException("Matrices must be the same length.");
                if (x[i].Length != y[i].Length)
                    throw new FormatException("Matrices must be the same size.");
            }
            double[][] ans = new double[x.Length][];
            for (int i = 0; i < x.Length; i++)
            {
                ans[i] = new double[i];
                for (int ii = 0; ii < x[i].Length; ii++)
                {
                    ans[i][ii] = x[i][ii] - y[i][ii];
                }
            }
            return ans;
        }
        public static double[][] MultiplyMatrixByNumber(double n, double[][] m)
        {
            if (m.Length == 0)
                throw new FormatException("Matrices cannot be 0.");
            for (int i = 0; i < m.Length - 1; i++)
            {
                if (m[i].Length != m[i + 1].Length)
                    throw new FormatException("Matrix rows must be the same length.");
            }
            double[][] ans = new double[m.Length][];

            for(int i=0; i<m.Length; i++ )
            {
                ans[i] = new double[m[0].Length];
                for(int j=0; j<m[0].Length; j++ )
                {
                    ans[i][j] = m[i][j] * n;
                }
            }
            return ans;
        }
        public static double[][] MultiplyMatrices(double[][] x, double[][] y)
        {
            if (x.Length == 0 || y.Length == 0 || x[0].Length == 0 || y[0].Length == 0)
            {
                throw new FormatException("Matrix cannot be 0.");
            }
            if (x[0].Length != y.Length)
            {
                throw new FormatException("Matrices not multipliable.");
            }
            for (int i = 0; i < x.Length - 1; i++)
            {
                if (x[i].Length != x[i + 1].Length)
                    throw new FormatException("Matrix rows must be the same length.");
                if (y[i].Length != y[i + 1].Length)
                    throw new FormatException("Matrix rows must be the same length.");
            }
            double[][] ans = new double[x.Length][];
            for (int i = 0; i < x.Length; i++)
            {
                ans[i] = new double[y[0].Length];
                for (int j = 0; j < y[0].Length; j++)
                {
                    ans[i][j] = 0;
                    for (int k = 0; k < x[0].Length; k++)
                    {
                        ans[i][j] += x[i][k] * y[k][j];
                    }
                }
            }
            return ans;
        }

        public static double MatrixDeterminant(double[][] m)
        {
            if (m.Length == 0)
                throw new FormatException("Matrix cannot be 0.");
            if (m.Length != m[0].Length)
                throw new FormatException("Matrix must be square.");
            for (int i = 0; i < m.Length - 1; i++)
            {
                if (m[i].Length != m[i + 1].Length)
                    throw new FormatException("Matrix rows must be the same length.");
            }
            if (m.Length == 2)
            {
                return m[0][0] * m[1][1] - m[0][1] * m[1][0];
            }
            else
            {
                double det = 0;
                for (int i = 0; i < m.Length; i++)
                {
                    double[][] m1 = new double[m.Length - 1][];
                    for (int j = 0; j < m.Length - 1; j++)
                    {
                        m1[j] = new double[m.Length - 1];
                    }
                    for (int j = 1; j < m.Length; j++)
                    {
                        int track = 0;
                        for (int k = 0; k < m.Length; k++)
                        {
                            if (k == i)
                                continue;
                            m1[j - 1][track] = m[j][k];
                            track++;

                        }
                    }
                    det += Math.Pow(-1, i) * m[0][i] * MatrixDeterminant(m1);
                }
                return det;
            }
        }

        public static double[][] InverseMatrix(double[][] m)
        {
            if (m.Length == 0)
                throw new FormatException("Matrix cannot be 0.");
            if (m.Length != m[0].Length)
                throw new FormatException("Matrix must be square.");
            for (int i = 0; i < m.Length - 1; i++)
            {
                if (m[i].Length != m[i + 1].Length)
                    throw new FormatException("Matrix rows must be the same length.");
            }
            double[][] m1 = new double[m.Length][];
            for (int i = 0; i < m.Length; i++)
            {
                m1[i] = new double[m[i].Length];
            }
            for (int i = 0; i < m.Length; i++)
            {
                for (int j = 0; j < m.Length; j++)
                {
                    int track1 = 0;
                    double[][] mhol = new double[m.Length - 1][];
                    for (int k = 0; k < mhol.Length; k++)
                        mhol[k] = new double[m[0].Length - 1];
                    for (int k = 0; k < m.Length; k++)
                    {
                        int track2 = 0;
                        if (k == i)
                            continue;
                        for (int l = 0; l < m[k].Length; l++)
                        {
                            if (l == j)
                            {
                                continue;
                            }
                            mhol[track1][track2] = m[k][l];
                            track2++;
                        }
                        track1++;
                    }
                    m1[i][j] = (double)MatrixDeterminant(mhol);
                }
            }
            double[][] cofactor = new double[m.Length][];
            int a = 1;
            for (int i = 0; i < m.Length; i++)
            {
                cofactor[i] = new double[m[i].Length];
                if (i % 2 == 1)
                    a = 1;
                else
                    a = 0;
                for (int j = 0; j < m[i].Length; j++)
                {
                    cofactor[i][j] = 1 * Math.Pow(-1, j + a);
                }
            }
            for (int x = 0; x < m1.Length; x++)
            {
                for (int y = 0; y < m1[x].Length; y++)
                {
                    m1[x][y] = m1[x][y] * cofactor[x][y];
                }
            }
            double[][] m2 = new double[m.Length][];
            for (int i = 0; i < m.Length; i++)
            {
                m2[i] = new double[m[i].Length];
            }
            for (int i = 0; i < m.Length; i++)
            {
                for (int j = 0; j < m[i].Length; j++)
                {
                    m2[j][i] = m1[i][j];
                }
            }
            return MultiplyMatrixByNumber(1 / (double)MatrixDeterminant(m), m2);
        }
        //-------------------------------------------Equation Solving-----------------------------------------
        public static double[][] SolveByMatrix(double[][] m, double[][] a)
        {
            if (m.Length == 0 || a.Length == 0)
                throw new FormatException("Matrices must both have values.");
            if (m.Length != a.Length)
                throw new FormatException("Matrices must have the same number of columns.");
            for (int i = 0; i < m.Length - 1; i++)
            {
                if (m[i].Length != m[i + 1].Length)
                    throw new FormatException("Matrix rows must be same length.");
            }
            for (int i = 0; i < a.Length; i++)
            {
                if (a[i].Length != 1)
                    throw new FormatException("Answer Matrix must only have one row.");
            }
            if (m.Length != m[0].Length)
                throw new FormatException("Equation Matrix must be square.");
            double[][] minverse = InverseMatrix(m);
            double[][] ans = MultiplyMatrices(minverse, a);
            return ans;
        }
    //-------------------------------------------GCD AND LCM-----------------------------------------
    public static int GCD(int a, int b)
        {
            if (a == 0)
                return b;
            if (b == 0)
                return a;
            if(b>a)
            {
                var hol = b;
                b = a;
                a = hol;
            }
            int newb = a % b;
            return GCD(b, newb);
        }

        public static int LCM(int a, int b)
        {
            return (a * b) / GCD(a, b);
        }
        //------------------------------------------Fractions---------------------------------------
        public static int[] DecimalToFraction(double x)
        {
            int count = 0;
            while (x % 1 != 0)
            {
                x *= 10;
                count++;
            }
            int b = (int)Math.Pow(10, count);
            return SimplifyFraction(new int[] { (int)x, b });
        }

        public static int[] SimplifyFraction(int[] x)
        {
            if (x.Length != 2)
                throw new FormatException("Array must contain two and only two values.");
            int[] copy = new int[] { x[0], x[1] };
            int g = GCD(copy[0], copy[1]);
            copy[0] /= g;
            copy[1] /= g;
            return copy;
        }
        public static double FractionToDecimal(int[] x)
        {
            if (x.Length != 2)
                throw new FormatException("Array must contain two and only two values.");
            return x[0] / x[1];
        }

        public static double[] AddFractions(double[] x, double[] y)
        {
            if (x.Length != 2 || y.Length != 2)
                throw new FormatException("Arrays must contain two and only two values.");
            if (x[1] == y[1])
            {
                return new double[] {x[0] + y[0], x[1] };
            }
            else
            {
                double bottom = x[1] * y[1];
                double top = x[0] * y[1] + y[0] * x[1];
                if (bottom % 1 == 0 && top % 1 == 0)
                {
                    int[] hold = SimplifyFraction(new int[] {(int) top, (int) bottom });
                    return new double[] {hold[0], hold[1]};
                }
                else return new double[] {top, bottom};
                    
            }
        }

        public static double[] SubtractFractions(double[] x, double[] y)
        {
            if (x.Length != 2 || y.Length != 2)
                throw new FormatException("Arrays must contain two and only two values.");
            if (x[1] == y[1])
            {
                return new double[] { x[0] - y[0], x[1] };
            }
            else
            {
                double bottom = x[1] * y[1];
                double top = x[0] * y[1] - y[0] * x[1];
                if (bottom % 1 == 0 && top % 1 == 0)
                {
                    int[] hold = SimplifyFraction(new int[] { (int)top, (int)bottom });
                    return new double[] { hold[0], hold[1] };
                }
                else return new double[] { top, bottom };
            }
        }

        public static double[] MultiplyFractions(double[] x, double[] y)
        {
            if (x.Length != 2 || y.Length != 2)
                throw new FormatException("Arrays must contain two and only two values.");
            double top = x[0] * y[0];
            double bottom = x[1] * y[1];
            if (bottom % 1 == 0 && top % 1 == 0)
            {
                int[] hold = SimplifyFraction(new int[] { (int)top, (int)bottom });
                return new double[] { hold[0], hold[1]};
            }
            else return new double[] { top, bottom };
        }

        public static double[] DivideFractions(double[] x, double[] y)
        {
            if (x.Length != 2 || y.Length != 2)
                throw new FormatException("Arrays must contain two and only two values.");
            double top = x[0] * y[1];
            double bottom = x[1] * y[0];
            if (bottom % 1 == 0 && top % 1 == 0)
            {
                int[] hold = SimplifyFraction(new int[] { (int)top, (int)bottom });
                return new double[] { hold[0], hold[1] };
            }
            else return new double[] { top, bottom };
        }
        //------------------------------------------Factoring-----------------------------------------
        public static int[] Factors(int x)
        {
            List<int> l = new List<int>();
            for (int i = 1; i < Math.Sqrt(x); i++)
            {
                if (x % i == 0)
                    l.Add(i);
            }
            l.Sort();
            return l.ToArray();
        }

        public static int[] PrimeFactors(int x)
        {
            List<int> l = new List<int>();
            while (x % 2 == 0)
            {
                l.Add(2);
                x /= 2;
            }
            for (int i = 3; i <= Math.Sqrt(x); i += 2)
            {
                while (x % i == 0)
                {
                    l.Add(i);
                    x /= i;
                }
            }
            if (x > 2)
            {
                l.Add(x);
            }
            l.Sort();
            return l.ToArray();
        }

        //-------------------------------------------Point calculations-----------------------------------------
        public static double distanceBetweenPoints(double[] a1, double[] a2)
        {
            return Math.Sqrt(Math.Pow(a2[0] - a1[0], 2) + Math.Pow(a2[1] - a1[1], 2));
        }
        public static double Slope(double[] a1, double[] a2)
        {
            return (a2[1] - a1[1]) / (a2[0] - a1[0]);
        }
        public static double[] MidPoint(double[] a1, double[] a2)
        {
            double[] ans = new double[2];
            ans[0] = (a1[0] + a2[0]) / 2;
            ans[1] = (a1[1] + a2[1]) / 2;
            return ans;
        }
        public static double[] slopeInterceptGivenPoints(double[] a1, double[] a2)
        {
            double slope = Slope(a1, a2);
            double b = a1[1] - (slope * a1[0]);
            return new double[] { slope, b };
        }
        //-------------------------------------------Statistics operations-----------------------------------------
        public static double standardDeviationFull(double[] m)
        {
            double mean = meanAverage(m);
            double value = 0;
            for(int i = 0; i < m.Length; i++)
            {
                value += Math.Pow((m[i] - mean),2);
            }
            return Math.Sqrt(value / m.Length);
        }

        public static double standardDeviationSample(double[] m)
        {
            double mean = meanAverage(m);
            double value = 0;
            for (int i = 0; i < m.Length; i++)
            {
                value += Math.Pow((m[i] - mean), 2);
            }
            return Math.Sqrt(value / (m.Length-1));
        }

        public static double zScore(double x, double m, double d)
        {
            return (x - m) / d;
        }

        public static double zScoreWithValues(double x, double[] a)
        {
            return (x - meanAverage(a)) / standardDeviationFull(a);
        }
        //-------------------------------------------Factorials-----------------------------------------
        public static int Factorial(int n)
        {
            int factn = 1;
            for (int i = 2; i <= n; i++)
            {
                factn *= i;
            }
            return factn;
        }
        //-------------------------------------------Perms and Combos-----------------------------------------
        public static int Permutation(int n, int r)
        {
            int ans = 1;
            for(int i = n; i > (n-r); i--)
            {
                ans *= i;
            }
            return ans;
        }

        public static int Combination(int n, int r)
        {
            int ans = 1;
            for (int i = n; i > (n - r); i--)
            {
                ans *= i;
            }
            int factr = Factorial(r);
            return ans / factr;
        }
        //-------------------------------------------Probability and Percent-----------------------------------------
        public static double Probability(double s, double p)
        {
            return s / p;
        }

        public static double Percent(double? x, double? y, double? p)
        {
            if(x is null)
            {
                if (y is null || p is null)
                    throw new FormatException("Only one value can be null.");
                return (double)y / ((double)p/100);
            }
            if(y is null)
            {
                if (p is null || x is null)
                    throw new FormatException("Only one value can be null.");
                return ((double)p / 100) * (double)x;
            }
            if (y is null || x is null)
                throw new FormatException("Only one value can be null.");
            return (double)y / (double)x;
        }
        //-------------------------------------------Shape operations-----------------------------------------
        // Triangle
        public static double areaOfTriangle(double b, double h)
        {
            return 0.5 * b * h;
        }
        //Rectangle
        public static double areaOfRectangle(double l, double w)
        {
            return l * w;
        }
        public static double perimeterOfRectangle(double l, double w)
        {
            return (2*l) + (2*w);
        }
        //Circle
        public static double areaOfCircleRadius(double r)
        {
            return Math.PI * Math.Pow(r, 2);
        }
        public static double areaOfCircleDiameter(double d)
        {
            return Math.PI * Math.Pow(d / 2, 2);
        }
        public static double circumferenceOfCircleRadius(double r)
        {
            return 2 * Math.PI * r;
        }
        public static double circumferenceOfCircleDiameter(double d)
        {
            return 2 * Math.PI * (d / 2);
        }
        public static double lengthOfArc(double a, double r)
        {
            return (a/360)*2*Math.PI*r;
        }
        public static double areaOfCircleSector(double a, double r)
        {
            return Math.PI * Math.Pow(r, 2) * (a / 360);
        }
        //3d shapes
        public static double volumeOfBox(double l, double w, double h)
        {
            return l * w * h;
        }
        public static double surfaceOfBox(double l, double w, double h)
        {
            return (2 * l * w) + (2 * w * h) + (2 * l * h);
        }
        public static double volumeOfSphereRadius(double r)
        {
            return (4 / 3) * Math.PI * Math.Pow(r, 3);
        }
        public static double volumeOfSphereDiameter(double d)
        {
            return (4 / 3) * Math.PI * Math.Pow(d/2, 3);
        }
        public static double volumeOfCylinder(double r, double h)
        {
            return Math.PI * Math.Pow(r, 2) * h;
        }
        public static double volumeOfCone(double r, double h)
        {
            return (1 / 3) * Math.PI * Math.Pow(r, 2) * h;
        }
        public static double volumeOfPyramid(double b, double h)
        {
            return (1 / 3) * b * h;
        }
        //Other shapes
        public static double areaOfTrapezoid(double b1, double b2, double h)
        {
            return ((b1 + b2) / 2) * h;
        }
        public static double areaOfParallelogram(double b, double h)
        {
            return b * h;
        }
        public static double sumOfInteriorAnglesOfPolygon(double n)
        {
            return (n - 2) * 180;
        }
        //-------------------------------------------Pythagorean Theorem-----------------------------------------
        public static double pythMissingAorB(double ab, double c)
        {
            return Math.Sqrt(Math.Pow(c, 2) - Math.Pow(ab, 2));
        }
        public static double pythMissingC(double a, double b)
        {
            return Math.Pow(a, 2) + Math.Pow(b, 2);
        }

        //-------------------------------------------Quadratic Formula and equation-----------------------------------------
        public static Object[] quadraticEquation(double a, double b, double c)
        {
            double sq = Discriminant(a, b, c);
            double bottom = 2 * a;
            Object[] ans = new Object[3];
            if (sq < 0)
            {
                ans[0] = (-b) / bottom;
                ans[1] = Math.Sqrt(-sq) / bottom;
                ans[2] = true;
            }
            else
            {
                ans[0] = ((-b + Math.Sqrt(sq)) / bottom);
                ans[1] = ((-b - Math.Sqrt(sq)) / bottom);
                ans[2] = false;
            }
            return ans;
        }

        public static double Discriminant(double a, double b, double c)
        {
            return Math.Pow(b, 2) - (4 * a * c);
        }

        public static double[] vertexOfParabola(double a, double b, double c)
        {
            double[] ans = new double[2];
            ans[0] = (-b) / (2 * a);
            ans[1] = (-Discriminant(a,b,c)) / (4*a);
            return ans;
        }
        //------------------------------------------Averages-----------------------------------------
        public static double meanAverage(double[] a)
        {
            double total = 0;
            for (int i = 0; i < a.Length; i++)
            {
                total = total + a[i];
            }
            return total / a.Length;
        }

        public static double[] modeAverage(double[] a)
        {
            var map = new Dictionary<double, int>();
            for (int i = 0; i < a.Length; i++)
            {
                if (map.ContainsKey(a[i]))
                {
                    map[a[i]] += 1;
                }
                else
                {
                    map.Add(a[i], 1);
                }
            }
            if (map.Count == a.Length)
            {
                return a;
            }
            int highest = 0;
            List<double> arr = new List<double>();
            foreach (var i in map) {
                if(i.Value > highest)
                {
                    highest = i.Value;
                    arr.Clear();
                }
                arr.Add(i.Value);
            }
            return arr.ToArray();
        }

        public static double medianAverage(double[] a)
        {
            double[] acopy = new double[a.Length];
            Array.Copy(a, acopy, a.Length);
            Array.Sort(acopy);
            if (a.Length % 2 == 0)
            {
                return (acopy[a.Length / 2] + acopy[a.Length / 2 + 1]) / 2;
            }
            else
            {
                return acopy[(int)(a.Length / 2)];
            }
        }

        public static double rangeAverage(double[] a)
        {
            double lowest = Double.MaxValue;
            double highest = Double.MinValue;
            for(int i=0; i< a.Length; i++)
            {
                if (a[i] < lowest)
                    lowest = a[i];
                if (a[i] > highest)
                    highest = a[i];
            }
            return highest - lowest;
        }
        //------------------------------------------Distance/Rate/Time-----------------------------------------
        public static double DistanceRateTime(double? r, double? t, double? d)
        {
            if (r is null)
            {
                if (d is null || t is null)
                    throw new FormatException("Only one value can be null.");
                return (double)d / (double)t;
            }
            if(t is null)
            {
                if (d is null || r is null)
                    throw new FormatException("Only one value can be null.");
                return (double)d / (double)r;
            }
            if (r is null || t is null)
                throw new FormatException("Only one value can be null.");
            return (double)r * (double)t;
        }
    }
}
