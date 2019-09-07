using System;
using System.Numerics;
using System.Collections.Generic;

namespace MathExpanded
{
    public static class Calc
    {
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
        //-------------------------------------------Perms and Combos-----------------------------------------
        public static double Permutation(double n, double r)
        {
            double factn = 1;
            for(int i = 2; i <= n; i++)
            {
                factn = factn * i;
            }
            double factrn = 1;
            for(int i=2; i<= (n-r); i++)
            {
                factrn = factrn * i;
            }
            return factn / factrn;
        }

        public static double Combination(double n, double r)
        {
            double factn = 1;
            for(int i=2; i<=n; i++)
            {
                factn = factn * i;
            }
            double factr = 1;
            for(int i = 2; i<=r; i++)
            {
                factr = factr * i;
            }
            double factrn = 1;
            for (int i = 2; i <= (n - r); i++)
            {
                factrn = factrn * i;
            }
            return factn / (factr * factrn);
        }
        //-------------------------------------------Shape operations-----------------------------------------
        public static double areaOfTriangle(double b, double h)
        {
            return 0.5 * b * h;
        }
        public static double areaOfRectangle(double l, double w)
        {
            return l * w;
        }
        public static double perimeterOfRectangle(double l, double w)
        {
            return (2*l) + (2*w);
        }
        public static double areaOfParallelogram(double b, double h)
        {
            return b * h;
        }
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
        public static double areaOfCircleSector(double c, double r)
        {
            return Math.PI * Math.Pow(r, 2) * (c / 360);
        }
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
        public static double volumeOfCylinder(double h, double r)
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
        public static double areaOfTrapezoid(double b1, double b2, double h)
        {
            return ((b1 + b2) / 2) * h;
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
    }
}
