namespace SimChA;

public static class Extensions
{
    public static IEnumerable<T> Shuffle<T>(this IEnumerable<T> source, Random rnd)
        => source.OrderBy(_ => rnd.Next());

    public static bool CoinFlip(this Random rnd)
        => rnd.Next(0, 2) == 0;
}