// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.EventData;

[Serializable]
public record CNEventPars(CNEventType Type, double Prob = 1, Dictionary<string, double>? Pars = null) : IHasProb
{
    public double GetDouble(string name)
    {
        if (Pars is null || !Pars.ContainsKey(name))
        {
            throw new Exception($"Double parameter {name} not found in {Type} event.");
        }
        return Pars[name];
    }
    
    public int GetInt(string name)
    {
        if (Pars is null || !Pars.ContainsKey(name))
        {
            throw new Exception($"Int parameter {name} not found in {Type} event.");
        }
        return (int) Math.Round(Pars[name]);
    }

    public long GetLong(string name)
    {
        if (Pars is null || !Pars.ContainsKey(name))
        {
            throw new Exception($"Long parameter {name} not found in {Type} event.");
        }
        return (long) Math.Round(Pars[name]);
    }
}
