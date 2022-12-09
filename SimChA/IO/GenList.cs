// Created by Felix Schifferdecker, 2022, felix.schifferdecker@gmx.de
using SimChA.DataTypes;

namespace SimChA.IO;

public static class GenList
{
    public static Dictionary<ChromNum, List<Gen>> TsgOgList = new Dictionary<ChromNum, List<Gen>>();

    public static Dictionary<ChromNum, List<Gen>> EssentialList = new Dictionary<ChromNum, List<Gen>>();
}