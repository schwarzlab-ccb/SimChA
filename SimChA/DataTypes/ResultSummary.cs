// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct ResultSummary
{
    public int Generations;
    public int ClusterCount;
    public long AliveCount;
    public long TotalCount;
    public int TreeDepth;
    public int NodeCount;
    public int LeafCount;
    public float Branching;
    public int SubcloneTotal;
    public int SubcloneSelect;
    // public double treeBalance;
    public double TreeBalance;
    // public double clonalDiversity;
    public double ClonalDiversity;
    // public double meanDriversPerCell;
    public double MeanDriversPerCell;

    public override string ToString()
        => $"{AliveCount},{TotalCount},{ClusterCount},{Generations},{TreeDepth},{NodeCount},{LeafCount},{Branching}," +
           $"{SubcloneTotal},{SubcloneSelect},{ClonalDiversity}," +
           $"{TreeBalance},{MeanDriversPerCell}";

    public static string Header()
        => "aliveCount,totalCount,clusterCount,generations,treeDepth,nodeCount,leafCount,branching," +
           "subcloneTotal,subcloneSelect,clonalDiversity,treeBalance,meanDriversPerCell";

    public string ToText()
        => "\t" + string.Join(",\n\t", Header().Split(",").Zip(ToString().Split(","), (label, val) => $"{label}: {val}"));
}
