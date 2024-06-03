// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com
using MathNet.Numerics.Statistics;
namespace SimChA.Optimization;

public class StatisticMeasures
{
    // The first Wasserstein Distance is obtained by calculating CDFs for both lists
    // and then the integral of the absolute difference between the CDFs.
    // https://arxiv.org/abs/1509.02237
    public static double WassersteinDistance(List<double> A, List<double> B)
    {
        // Calculate CDFs
        var cdfA = GetCDF(A);
        var cdfB = GetCDF(B);

        // Calculate 1-Wasserstein distance
        double distance = 0.0;
        for (int i = 0; i < A.Count; i++)
        {
            distance += Math.Abs(cdfA[i] - cdfB[i]);
        }
        return distance / A.Count;
    }

    public static double WassersteinDistance(Histogram histA, Histogram histB)
    {
        // Calculate CDFs
        var cdfA = GetCDF(histA);
        var cdfB = GetCDF(histB);

        // Calculate 1-Wasserstein distance
        double distance = 0.0;
        for (int i = 0; i < histA.BucketCount; i++)
        {
            distance += Math.Abs(cdfA[i] - cdfB[i]) * histA[i].Width;
        }
        return distance / histA.BucketCount;
    }


    public static List<double> GetCDF(List<double> list)
    {
        var cdf = new List<double>();
        double sum = list.Sum();
        double cumulative = 0;
        foreach (var item in list)
        {
            cumulative += item / sum;
            cdf.Add(cumulative);
        }
        return cdf;
    }

    public static List<double> GetCDF(Histogram hist)
    {
        var cdf = new List<double>();
        double cumulative = 0;
        for (int i = 0; i < hist.BucketCount; i++)
        {
            cumulative += hist[i].Count / hist.DataCount;
            cdf.Add(cumulative);
        }
        return cdf;
    }
}