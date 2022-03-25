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
    public double treeBalance;
    public double treeBalanceFiltered;
    public double clonalDiversity;
    public double clonalDiversityFiltered;
    public double meanDriversPerCell;
    public double meanDriversPerCellFiltered;

    public override string ToString()
        => $"{AliveCount},{TotalCount},{Generations},{TreeDepth},{NodeCount},{LeafCount},{Branching}," +
           $"{SubcloneTotal},{SubcloneSelect},{clonalDiversity},{clonalDiversityFiltered},{treeBalance}," +
           $"{treeBalanceFiltered},{meanDriversPerCell},{meanDriversPerCellFiltered}";

    public static string Header()
        => "aliveCount,totalCount,generations,treeDepth,nodeCount,leafCount,branching," +
           "subcloneTotal,subcloneSelect,clonalDiversity,clonalDiversityFiltered,treeBalance," +
           "treeBalanceFiltered,meanDriversPerCell,meanDriversPerCellFiltered";

    public string ToLine()
        => string.Join(", ", Header().Split(",").Zip(ToString().Split(","), (label, val) => $"{label}:{val}"));
}
