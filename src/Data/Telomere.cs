namespace SimChA.Data;

// A telomere is attached to the region that contains the complete reference-terminal
// interval. Region edits discard it when any part is lost, so only intact telomeres
// can make a contig end eligible for a telomere-bound event.
public class Telomere(long start, long end, string chrom) : GenRange(start, end, chrom);
