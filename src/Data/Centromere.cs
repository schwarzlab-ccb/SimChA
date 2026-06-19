namespace SimChA.Data;

// A centromere is attached to the region it falls within and rides along with that region
// through structural operations, in the same way genes do.
public class Centromere(long start, long end, string chrom) : GenRange(start, end, chrom);
