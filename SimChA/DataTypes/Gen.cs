// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de

namespace SimChA.DataTypes;


public class Gen {
    public string name { get; set;}

    // TODO change below to a region
    public ChromNum chr { get; set;}
    public int start { get; set;}
    public int stop { get; set;}
    public float deltaFitness { get; set;}
}