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
    public float treeBalance;
    public float treeBalanceConnected;
    public float clonalDiversity;
    public float clonalDiversityFiltered;

    public override string ToString()
        => $"{AliveCount},{TotalCount},{Generations},{TreeDepth},{NodeCount},{LeafCount},{Branching}," +
           $"{SubcloneTotal},{SubcloneSelect},{clonalDiversityFiltered},{clonalDiversity},{treeBalance},{treeBalanceConnected}";

    public static string Header()
        => "aliveCount,totalCount,generations,treeDepth,nodeCount,leafCount,branching,subcloneTotal," +
           "subcloneSelect,clonalDiversityFiltered,clonalDiversity,treeBalance,treeBalanceConnected";

    public string ToLine()
        => string.Join(", ", Header().Split(",").Zip(ToString().Split(","), (label, val) => $"{label}:{val}"));
}