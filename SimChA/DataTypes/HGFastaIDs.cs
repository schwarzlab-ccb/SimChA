// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public static class HGFastaIDs
{
    public const int CHR_COUNT = 46;
    public const int AUTOSOME_COUNT = 22;
    
    public static string Sex(bool isXX) => isXX ? "XX" : "XY";

    // https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
    public static readonly Dictionary<string, ChrNo> HG38IDs = new Dictionary<string, ChrNo>()
    {
        {"CM000663.2", ChrNo.chr1},
        {"CM000664.2", ChrNo.chr2},
        {"CM000665.2", ChrNo.chr3},
        {"CM000666.2", ChrNo.chr4},
        {"CM000667.2", ChrNo.chr5},
        {"CM000668.2", ChrNo.chr6},
        {"CM000669.2", ChrNo.chr7},
        {"CM000670.2", ChrNo.chr8},
        {"CM000671.2", ChrNo.chr9},
        {"CM000672.2", ChrNo.chr10},
        {"CM000673.2", ChrNo.chr11},
        {"CM000674.2", ChrNo.chr12},
        {"CM000675.2", ChrNo.chr13},
        {"CM000676.2", ChrNo.chr14},
        {"CM000677.2", ChrNo.chr15},
        {"CM000678.2", ChrNo.chr16},
        {"CM000679.2", ChrNo.chr17},
        {"CM000680.2", ChrNo.chr18},
        {"CM000681.2", ChrNo.chr19},
        {"CM000682.2", ChrNo.chr20},
        {"CM000683.2", ChrNo.chr21},
        {"CM000684.2", ChrNo.chr22},
        {"CM000685.2", ChrNo.chrX},
        {"CM000686.2", ChrNo.chrY}
    };

}