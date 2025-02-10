using SimChA.EventData;

namespace SimChA.Data;

public struct Clone
{
    public string CloneId { get; }
    public Karyotype Karyotype { get; }
    public List<CNEventDesc> Events { get; }

    public Clone(string cloneId, Karyotype kar, List<CNEventDesc> events)
    {
        CloneId = cloneId;
        Karyotype = kar;
        Events = events;
    }
}