namespace SimChA.Data;

// @CODY IMO SNV should hold either both REF and ALT or neither of them (i.e. sampling should be either base aware, or done during output)
public record SNV(long Pos, string Chrom, Nucleotide Alt);