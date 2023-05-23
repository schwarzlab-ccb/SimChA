// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.EventData;

[Serializable]
public record CNEventPars(CNEventType Type, double Prob = 1, Dictionary<string, double>? Pars = null) : IHasProb
{
    public double Get(string name, double defaultValue)
    {
        if (Pars is null || !Pars.ContainsKey(name))
        {
            return defaultValue;
        }
        return Pars[name];
    }
    
    public int Get(string name, int defaultValue)
    {
        if (Pars is null || !Pars.ContainsKey(name))
        {
            return defaultValue;
        }
        return (int) Math.Round(Pars[name]);
    }

    public long Get(string name, long defaultValue)
    {
        if (Pars is null || !Pars.ContainsKey(name))
        {
            return defaultValue;
        }
        return (long) Math.Round(Pars[name]);
    }
}
