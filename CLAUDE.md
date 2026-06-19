# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Build
dotnet build

# Run (default config, outputs to ./out)
dotnet run

# Run with options
dotnet run -- --config ./configs/main_config.json -m evolution -O ./out

# Run all tests
dotnet test

# Run tests in the test project explicitly
dotnet test test
```

The `.csproj` is at the repo root; sources live in `src/`. `dotnet run` runs the simulator directly without specifying a project.

## Architecture

### Data model (bottom-up)

**`Region`** (`src/Data/Region.cs`) — an oriented genomic segment: zero-indexed, start-inclusive, end-exclusive (`[start, end)`). Tracks chromosome name, haplotype (`Hap1`), SNVs, gene annotations, and centromeres. `AbsStart`/`AbsEnd` are absolute genome coordinates. Genes and centromeres are attached at creation (`RefGen.GetRegion`) and ride along through structural operations — copied with the region and dropped by `UpdateRegion` when no longer fully `IsInsideOf` the (only-ever-shrinking) region. They surface as `Contig.Genes` / `Contig.Centromeres` (cached, invalidated by the `Regions` setter).

**`Contig`** (`src/Data/Contig.cs`) — an ordered list of `Region`s representing a derived chromosome. Contigs may span regions from multiple original chromosomes (e.g. after a translocation). All structural operations (`DeleteRange`, `DuplicateRange`, `InvertRange`, `Split`, `Join`, `Scatter`, `Bridge`, etc.) are delegated to `RegionOps` (`src/Computation/RegionOps.cs`).

**`Karyotype`** (`src/Data/Karyotype.cs`) — the full genome state: a list of `Contig`s plus a `GeneCounts[geneType][geneId]` matrix maintained incrementally. Every `Apply*` method follows the pattern: `RemoveGenes(contig)` → mutate → `AddGenes(contig)`. Empty contigs are kept in the list to preserve stable integer indices.

**`RefGen`** (`src/Data/RefGen.cs`) — immutable reference data loaded at startup: chromosome lengths, centromere positions, and gene lists (TSG, OG, Essentials) per sex.

### Simulation flow

`Program.cs` parses CLI options → reads config → calls `Factory.GetSimulator()` → calls `simulator.Simulate()` → scores samples → writes output.

**Simulator hierarchy** (`src/Simulation/`):
- `Simulator` — base class; MonteCarlo/basic mode: picks events uniformly at random, applies them without checking fitness.
- `EvoSimulator` — overrides `SampleEvents`; uses Metropolis-like acceptance: `exp(ΔFitness − Acceptance) > U(0,1)`.
- `MatchSimulator` — overrides `SampleEvents`; minimizes distance to a per-node target fitness, with a `Decay` parameter that tightens acceptance as events progress.

`Factory.GetSimulator()` returns the right subclass based on `SelectionMode` (`MonteCarlo` / `Evolution` / `FitnessMatching`).

### Event system

`CNEventType` (enum) → `CNEventPars` (record: type + Prob + Frac + Frag + Signature) → concrete `BaseEventData` subclasses in `src/EventData/` generate the specific genomic coordinates → `eventData.ApplyEvent(karyotype)` mutates the karyotype in place.

`Sampling.GenerateCNEventData()` dispatches to the right `EventData` constructor based on event type.

### Signatures and mixing

A config's `Signatures` list is flattened into a single `List<CNEventPars>` by `Factory.MixSignatures()`, scaling each event's probability by its parent signature's probability. The `MixtureType` (`Single` / `Constant` / `Dirichlet`) controls per-sample variation in signature weights.

### Fitness

`Fitness.Calculate()` (`src/Computation/Fitness.cs`) combines three weighted terms:
- **Stress**: penalizes genome length above diploid reference.
- **TsgOg**: log-scaled contribution of oncogenes (positive) and tumor suppressors (negative).
- **Essentiality**: penalizes homozygous loss of essential genes.

Gene counts are maintained incrementally in `Karyotype.GeneCounts` rather than recomputed from scratch.

### Execution modes

- **Repeats** (`-R N`): N independent samples, each rooted at the same diploid/tetraploid karyotype.
- **Tree** (`-T file`): samples share ancestry; simulation recurses depth-first through `CTreeNode` parents/children.
- **Profiles** (`-P file`): skips simulation; reads existing CN profiles and scores them.

### Config structure

`SimChAConfig` (JSON) has four sections:
- `SimParams` — seed, assembly, sex, mutation rate distribution, mixture type.
- `FitParams` — weights for stress/TsgOg/essentiality, gene set folder name.
- `EvoParams` — acceptance threshold, max tries, decay (required for evolution/matching modes).
- `Signatures` — array of named signature objects, each with a `Prob` and `Events` array.

Pre-built configs are in `configs/` (hg19 optimized), `configs/hg38/`, and `configs/basic/`. Cancer type abbreviations follow TCGA conventions.
