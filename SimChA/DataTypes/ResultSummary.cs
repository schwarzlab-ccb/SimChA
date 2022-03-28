// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.Computation;

namespace SimChA.DataTypes;

public struct ResultSummary
{
    public int RepeatId;
    public int PopCount;
    public int GenerationId;
    public int Generations;
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
        =>
            $"{RepeatId},{GenerationId},{PopCount},{AliveCount},{TotalCount},{Generations},{TreeDepth}," +
            $"{NodeCount},{LeafCount},{Branching},{SubcloneTotal},{SubcloneSelect},{ClonalDiversity}," + 
            $"{TreeBalance},{MeanDriversPerCell}";

    public static string Header()
        => "RepeatId,GenerationId,popCount,aliveCount,totalCount,generations,treeDepth,nodeCount,leafCount,branching," +
           "subcloneTotal,subcloneSelect,clonalDiversity,treeBalance,meanDriversPerCell";

    public string ToText()
        => "\t" + string.Join(",\n\t",
            Header().Split(",").Zip(ToString().Split(","), (label, val) => $"{label}: {val}"));


    public ResultSummary(int repeatId, int generationId, ParentTree connectedTree, List<SubClone> aboveCutOff,
        int cloneCount, int sampleCount, int stepNo, List<(long total, long alive)> popSizes, int popCount)
    {
        RepeatId = repeatId;
        GenerationId = generationId;
        PopCount = popCount;
        (NodeCount, LeafCount, TreeDepth, Branching)
            = TreeAnalysis.ComputeTreeSize(connectedTree);
        TreeBalance = TreeAnalysis.ComputeTreeBalance(LeafCount, connectedTree);
        ClonalDiversity = TreeAnalysis.ComputeClonalDiversity(aboveCutOff);
        MeanDriversPerCell = TreeAnalysis.ComputeMeanDriversPerCell(aboveCutOff);
        SubcloneTotal = cloneCount;
        SubcloneSelect = sampleCount;
        Generations = stepNo;
        AliveCount = popSizes.Last().alive;
        TotalCount = popSizes.Last().total;
    }
}