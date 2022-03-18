// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct ResultSummary
{
    public long AliveCount;
    public long TotalCount;
    public int Generations;
    public int TreeDepth;
    public int NodeCount;
    public int LeafCount;
    public float Branching;
    public int SubcloneTotal;
    public int SubcloneSelect;
    public float clonalDiversityFiltered;
    public float clonalDiversity;

    public override string ToString()
        => $"{AliveCount},{TotalCount},{Generations},{TreeDepth},{NodeCount},{LeafCount},{Branching},{SubcloneTotal},{SubcloneSelect},{clonalDiversityFiltered},{clonalDiversity}";

    public static string Header()
        => "aliveCount,totalCount,generations,treeDepth,nodeCount,leafCount,branching,subcloneTotal,subcloneSelect,clonalDiversityFiltered,clonalDiversity";

    public string ToLine()
        => string.Join(", ", Header().Split(",").Zip(ToString().Split(","), (label, val) => $"{label}:{val}"));
}