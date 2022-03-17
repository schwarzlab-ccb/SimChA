// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct ResultSummary
{
    public long aliveCount;
    public long totalCount;
    public int generations;
    public int treeDepth;
    public int nodeCount;
    public int leafCount;
    public float branching;
    public int subcloneTotal;
    public int subcloneSelect;

    public string ToString()
        => $"{aliveCount},{totalCount},{generations},{treeDepth},{nodeCount},{leafCount},{branching},{subcloneTotal},{subcloneSelect}";

    public static string Header()
        => "aliveCount,totalCount,generations,treeDepth,nodeCount,leafCount,branching,subcloneTotal,subcloneSelect";

    public string ToLine()
        => string.Join(", ", Header().Split(",").Zip(ToString().Split(","), (label, val) => $"{label}:{val}"));
}