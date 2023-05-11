// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.EventData;

[Serializable]
public record CNEventPars(CNEventType Type, double Prob = 1, Dictionary<string, double>? Params = null) : IHasProb
{
    public double Get(string name, double defaultValue)
    {
        if (Params is null || !Params.ContainsKey(name))
        {
            return defaultValue;
        }
        return Params[name];
    }
    
    public int Get(string name, int defaultValue)
    {
        if (Params is null || !Params.ContainsKey(name))
        {
            return defaultValue;
        }
        return (int) Math.Round(Params[name]);
    }

    public long Get(string name, long defaultValue)
    {
        if (Params is null || !Params.ContainsKey(name))
        {
            return defaultValue;
        }
        return (long) Math.Round(Params[name]);
    }
}
