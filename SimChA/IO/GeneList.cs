// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de
using SimChA.DataTypes;

namespace SimChA.IO;

public static class GeneList
{
    public static readonly Dictionary<ChromNum, List<Gene>> TsgOgList = new();
    public static readonly Dictionary<ChromNum, List<Gene>> EssentialList = new();
}