namespace SimChA.Misc;

public static class Extensions
{
    public static IEnumerable<T> Shuffle<T>(this IEnumerable<T> source, Random rnd)
        => source.OrderBy(_ => rnd.Next());

    public static bool CoinFlip(this Random rnd)
        => rnd.Next(0, 2) == 0;
    
    public static int GetRandIndex<T>(this IEnumerable<T> source, Random rnd)
        => rnd.Next(source.Count());

    public static IEnumerable<int> GetRandIndices<T>(this IEnumerable<T> source, int count, Random rnd)
        => Enumerable.Range(0, source.Count()).Shuffle(rnd).Take(count);

    public static IEnumerable<T> GetValues<T>() where T : Enum 
        => Enum.GetValues(typeof(T)).Cast<T>();
    
    public static void ForEach<T>(this IEnumerable<T> sequence, Action<T> action)
    {
        foreach (var item in sequence)
        {
            action(item);
        }
    }
}