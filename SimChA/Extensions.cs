namespace SimChA;

public static class Extensions
{
    public static IEnumerable<T> Shuffle<T>(this IEnumerable<T> source)
    {
        var rnd = new Random();
        return source.OrderBy(_ => rnd.Next());
    }
    
    public static bool CoinFlip(this Random random)
        => random.Next(0, 2) == 0;
}