// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de

namespace SimChA.IO;


public class Abberation{
    public string CloneName;
    public string AbberationEnum;
    public string? Region;
    public float? DeltaFitness;
    public Abberation(string cloneName, string abberation, string region = "", float deltaFitness = 0){
        CloneName = cloneName;
        AbberationEnum = abberation;
        Region = region;
        DeltaFitness = deltaFitness;
        AbberationList.ListAbberation.Add(this);
    }
}

public static class AbberationList
{    
    public static List<Abberation> ListAbberation = new List<Abberation>();
}
