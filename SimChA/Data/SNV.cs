namespace SimChA.Data;

// TODO @CODY IMO SNV should hold both REF and ALT
public record SNV(long Pos, string Chrom, Nucleotide Alt)
{
    public string Header() 
        => "chrom\tpos\talt";
    
    public override string ToString() 
        => $"{Chrom}\t{Pos}\t{Alt}";
}