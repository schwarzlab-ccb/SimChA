// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.Optimization;

public class StatisticMeasures<T> where T:IConvertible
{
    // The first Wasserstein Distance is obtained by calculating CDFs for both lists
    // and then the integral of the absolute difference between the CDFs.
    // https://arxiv.org/abs/1509.02237
    public static double WassersteinDistance(List<T> A, List<T> B)
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


    public static List<double> GetCDF(List<T> list)
    {
        var cdf = new List<double>();
        double sum = list.Sum(x => Convert.ToDouble(x));
        double cumulative = 0;
        foreach (var item in list)
        {
            cumulative += Convert.ToDouble(item) / sum;
            cdf.Add(cumulative);
        }
        return cdf;
    }
}