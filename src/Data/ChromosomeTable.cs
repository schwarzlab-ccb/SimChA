namespace SimChA.Data;

// The contents of an assembly's chromosomes.tsv: one row per chromosome carrying its length, the
// sex it belongs to, and the reference coordinates of its centromere and telomeres.
//
// Telomeres used to be synthesised from a hard-coded 50 kb constant rather than read, and
// centromeres lived in a separate centromeres.tsv. Both are now data, so an assembly is described
// by a single file.
public record ChromosomeTable(
    Dictionary<string, int> Lengths,
    Dictionary<string, SexType> Sex,
    Dictionary<string, GenRange> Centromeres,
    Dictionary<string, List<GenRange>> Telomeres);
