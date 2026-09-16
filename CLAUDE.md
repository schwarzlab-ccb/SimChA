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

## Code style

**Post-fix notation.** Prefer post-fix (fluent, extension-method) calls over prefix calls that wrap a
value in parentheses: write `value.ToCount()`, not `ToCount(value)`, and
`SampleContigWeighted(...).ToList()`, not `AsList(SampleContigWeighted(...))`. Chains read
left-to-right in the order they execute, and a new step appends to the end rather than nesting the
whole expression one level deeper. Helper conversions therefore go in `Extensions`
(`src/Computation/Extensions.cs`) — or as `private static ... (this T value)` members of the static
class that uses them — and are named for what they return (`ToList`, `ToCount`). Static factory and
sampling entry points that do not act *on* a receiver (`Sampling.GetExpSeg`, `Factory.GetSimulator`,
MathNet's `Gamma.Sample`) stay prefix.

## Architecture

### Data model (bottom-up)

**`Region`** (`src/Data/Region.cs`) — an oriented genomic segment: zero-indexed, start-inclusive, end-exclusive (`[start, end)`). Tracks chromosome name, haplotype (`Hap1`), SNVs, gene annotations, and centromeres. `AbsStart`/`AbsEnd` are absolute genome coordinates. Genes and centromeres are attached at creation (`RefGen.GetRegion`) and ride along through structural operations — copied with the region and dropped by `UpdateRegion` when no longer fully `IsInsideOf` the (only-ever-shrinking) region. They surface as `Contig.Genes` / `Contig.Centromeres` (cached, invalidated by the `Regions` setter).

**`Contig`** (`src/Data/Contig.cs`) — an ordered list of `Region`s representing a derived chromosome. Contigs may span regions from multiple original chromosomes (e.g. after a translocation). All structural operations (`DeleteRange`, `DuplicateRange`, `InvertRange`, `Split`, `Join`, `Scatter`, `Bridge`, etc.) are delegated to `RegionOps` (`src/Computation/RegionOps.cs`).

**`Karyotype`** (`src/Data/Karyotype.cs`) — the full genome state: a list of `Contig`s plus a `GeneCounts[geneType][geneId]` matrix maintained incrementally. Every `Apply*` method follows the pattern: `RemoveGenes(contig)` → mutate → `AddGenes(contig)`. Empty contigs are kept in the list to preserve stable integer indices.

**`RefGen`** (`src/Data/RefGen.cs`) — immutable reference data loaded at startup: chromosome lengths, centromere and telomere positions, and gene lists (TSG, OG, Essentials) per sex. Lengths, sex, centromeres and telomeres all come from a single `chromosomes.tsv` per assembly (7 columns: `chrom`, `length`, `sex`, `tel_p_end`, `cen_start`, `cen_end`, `tel_q_start`; `#` comments allowed, `.` marks an absent feature) parsed into a `ChromosomeTable`. The separate `centromeres.tsv` is gone, and telomere extents are data rather than the former hard-coded `RefGen.TelomereLength = 50_000`.

### Simulation flow

`Program.cs` parses CLI options → `SimChAConfig.Load()` (reads + resolves config) → calls `Factory.GetSimulator()` → calls `simulator.Simulate()` → scores samples → writes output. `Program.cs` itself is kept thin; configuration setup lives in `SimChAConfig.Load()` (see *Config structure*).

**Simulator hierarchy** (`src/Simulation/`):
- `Simulator` — base class; MonteCarlo/basic mode: picks events uniformly at random, applies them without checking fitness.
- `EvoSimulator` — overrides `SampleEvents`; accepts a proposal with the Glauber (Fermi) probability
  `1 / (1 + exp(−ΔFitness))` (`Fitness.AcceptProb`). This replaced a Metropolis
  `min(1, exp(ΔFitness))`, which clipped to exactly 1 for every `ΔFitness ≥ 0` and so had a zero
  derivative — with respect to every fitness weight — over a third of all accepted events at the
  fitted optimum. The two agree in the deleterious tail; the new form is strictly monotone in
  `ΔFitness` everywhere, so no event is selection-free (99.3% of accepted events land in the graded
  band `0.02 < p < 0.98` on the 2026-09-10 spice cohort). **There is no acceptance offset**: the
  rule is symmetric, a neutral proposal accepted half the time. `EvoParams.Acceptance` (delta) used
  to shift that midpoint and was removed — see `Fitness.AcceptProb` for why nothing could identify
  it, and what would have to change to make it fittable again.
- `MatchSimulator` — overrides `SampleEvents`; minimizes distance to a per-node target fitness, with a `Decay` parameter that tightens acceptance as events progress.

`Factory.GetSimulator()` returns the right subclass based on `SelectionMode` (`MonteCarlo` / `Evolution` / `FitnessMatching`).

The base `Simulator` wraps each sample's `SampleEvents` in `SampleEventsLimited`, which enforces `SimParams.MaxWGD` (default `-1` = no limit): the target event count is drawn once, and if a generated sample contains more whole-genome doublings than allowed it is re-simulated from the parent karyotype **with the same event count** (so the cap does not bias the event-count distribution). The count is only redrawn — and the per-count try budget reset — after `SimParams.MaxWgdTries` consecutive failures; a fixed tree distance (`Distance >= 0`) cannot be redrawn, so it aborts instead. This applies uniformly to all three modes.

### Event system

`CNEventType` (enum) → `CNEventPars` (record: type + Prob + Frac + Frag + Signature) → concrete `BaseEventData` subclasses in `src/EventData/` generate the specific genomic coordinates → `eventData.ApplyEvent(karyotype)` mutates the karyotype in place.

`Sampling.GenerateCNEventData()` dispatches to the right `EventData` constructor based on event type.

Internal, tail and telomere events take a separate path: `SampleTerminalArmWeighted` first picks one
`TerminalArm` (`src/Data/TerminalArm.cs`) across the whole karyotype, weighted by arm length, and the
constructor then draws the event against that arm. `TerminalArm.ArmLength` runs to the **middle** of
the nearest centromere, because that is the arm definition the empirical proportions are fitted
against upstream (`project-simcha/src_empirical_dist/lib_create_config.py`, `FIXED_BETA_ALPHA` and
`beta_shape_for_mean` mirror `Sampling.FixedBetaAlpha` / `Sampling.GetBetaShape`, which are the
*defaults* a config without a `Shape` gets); `UsableLength`
stops at the near edge of the centromere and bounds the event, so an over-long draw collapses onto a
whole-arm event instead of entering the centromere. Keep the two sides in step: changing the arm
denominator here invalidates every `Frac` in `configs/`.

### Event lengths

`Frac` sets the **mean** of an event's length distribution and `Shape` its dispersion; together they
determine it. Internal events use a bounded Pareto on [0,1], telomere-bound ones a Beta, everything
else an exponential that ignores `Shape`. A config without `Shape` gets `FixedParetoShape` (0.5) and
`FixedBetaAlpha` (0.6), so nothing written before the field existed changes behaviour.

The Beta side was always exact: `beta = alpha*(1-mean)/mean` gives `Beta(alpha, beta)` the requested
mean for any alpha, so alpha only ever moved the spread. The Pareto side was not, and the repair is
worth recording because two separate things were wrong.

`GetParetoScale` used to be `1/2*(-1 + sqrt(1 + 4*mean^2))`, which **contains no shape term** — it is
not the inverse of any mean, and at the configured 0.5 it left the realized proportion 10-13% below
the configured `Frac` (told 0.1058 it delivered 0.0948). Worse, the error was shape-dependent enough
to compress the classes together: interior gain and loss, configured 24% apart at 0.1058 and 0.1312,
both came out at 0.110 in a real cohort, so the per-class mean was not being controlled at all.
`GetParetoScale` now bisects `BoundedParetoMean` — which has a closed form — for the scale that
actually delivers the requested mean, cached per `(mean, shape)` because the pair is fixed for a run.

`SampleParetoLim` used to draw from the *unbounded* Pareto and reject anything above 1, with a
fallback after a thousand tries. That is the same distribution, but it costs 1.2 draws per sample at
shape 0.5 and 2.0 at 0.12, and it needs a fallback that silently piles mass at exactly 1. The bounded
Pareto has a closed-form quantile, `L*(1 - u*(1 - (L/(L+1))^a))^(-1/a) - L`, so one draw always
suffices and no fallback is needed.

Two numerical traps live in `BoundedParetoMean`, both of which produced wrong answers rather than
exceptions, and both of which have regression tests:

- **Do not use `double.LogP1` / `double.ExpM1`.** .NET's are the naive `Math.Log(1 + x)` and
  `Math.Exp(x) - 1`; at x = 1e-12 they are already wrong in the fifth digit, which made the mean
  return 1.9e9 for a quantity bounded by one half. `Sampling` carries Kahan's formulations instead.
- **The integral needs two branches.** Above scale 1, `(L+1)^(1-a)` and `L^(1-a)` agree to many
  digits and their difference cancels to noise (828 at L = 1e9), so the `expm1` identity is used;
  below 1 that identity would overflow, since `log((L+1)/L)` reaches 691 at the smallest scale
  searched, so the direct form is used — with `L^a` distributed into the bracket, because evaluating
  the two powers separately gives `0 * infinity` at 1e-300 once the shape passes one.

The bounded Pareto tends to uniform on [0,1] as its scale grows, so **no shape can reach a mean of
one half**; `GetParetoScale` throws rather than approximate. `GetParetoSeg` saturates at `Frac >= 1`
to the whole arm, mirroring `GetBetaSeg` and matching what the rejection loop did by exhaustion.

### Signatures and mixing

A config's `Signatures` list is flattened into a single `List<CNEventPars>` by `Factory.MixSignatures()`, scaling each event's probability by its parent signature's probability. The `MixtureType` (`Single` / `Constant` / `Dirichlet`) controls per-sample variation in signature weights.

### Fitness

`Fitness.Calculate()` (`src/Computation/Fitness.cs`) combines three weighted terms:
- **Stress**: penalizes genome length above diploid reference.
- **TsgOg**: log-scaled contribution of oncogenes (positive) and tumor suppressors (negative).
- **Essentiality**: the price of **haploinsufficiency**, and only that — the summed score of every
  essential gene sitting below its reference copy number (`Fitness.EssTerm`). Outright loss is not
  priced here at all: `FitParams.ProhibitEssentialLoss` (default true) makes `EvoSimulator` reject
  any proposal that would leave an essential gene at zero copies, before the fitness is computed
  (`Fitness.AnyEssentialLost`), because that is inviability rather than unfitness. The reference is
  sex-aware (`RefGen.SexGeneRefCNs`): two on the autosomes, one for a male's X and Y, so a male is
  not charged for being hemizygous on X.

  This replaced a graded `DosageLoss(CN)^k * score` whose exponent `HaploExponent` (k) set the CN==1
  cost *relative to* CN==0. With CN==0 unreachable that ratio has no referent — k was exactly
  degenerate with `Essentiality`, so only their product was determined — and it was retired. Before
  that the term was `min(CN-1, 0) * score`, a homozygous-loss-only step that was inactive in ~89% of
  simulated samples, which had made `Essentiality` the least determined parameter of every fit.

Gene counts are maintained incrementally in `Karyotype.GeneCounts` rather than recomputed from scratch.

### Execution modes

- **Repeats** (`-R N`): N independent samples, each rooted at the same diploid/tetraploid karyotype.
- **Tree** (`-T file`): samples share ancestry; simulation recurses depth-first through `CTreeNode` parents/children.
- **Profiles** (`-P file`): skips simulation; reads existing CN profiles and scores them.

### Config structure

`SimChAConfig` (JSON) has four sections plus two optional top-level fields:
- `SimParams` — seed, assembly, sex, mutation rate distribution, mixture type, `MaxWGD`, `MaxWgdTries` (re-simulations at a fixed event count before the count is redrawn; default 100).
- `FitParams` — weights for stress/TsgOg/essentiality, `ProhibitEssentialLoss` (reject an event that would zero an essential gene rather than pricing it; default true, and configs predating the field deserialise to that), gene set folder name. The retired `HaploExponent` is ignored if a config still carries it.
- `EvoParams` — max tries, decay (required for evolution/matching modes). Configs still carrying the removed `Acceptance` key deserialise fine; it is ignored.
- `Signatures` — array of named signature objects, each with a `Prob` and `Events` array.
- `Root` (top level, optional) — base directory for resolving relative paths.
- `Version` (top level) — ignored on input; stamped with the running version on output.

`SimChAConfig.Load(CmdOptions)` reads the config file and resolves `root`, `assembly`, `gene_set`, and
the RNG `seed` with precedence **command line > config > default** (`-r`/`--root`, `-a`/`--assembly`,
`-g`/`--gene_set`, `--seed`; a negative seed — from either source — draws a random one).
It applies the effective root as the working directory and stamps the effective `SimParams.Assembly`,
`SimParams.Seed`, `FitParams.GeneSet`, `Root`, and `Version` back into the returned config, so the
output `sim_params.json` records exactly what was used (and can be fed back in as input).

Path resolution in `FileIO.ReadGenRef`: an **absolute** assembly/gene-set value is used directly as the
folder; a **relative** value is resolved under the data folder (`data/<assembly>`) and assembly folder
(`<assembly>/<gene_set>`) respectively.

Pre-built configs are in `configs/` (hg19 optimized), `configs/hg38/`, and `configs/basic/`. Cancer type abbreviations follow TCGA conventions.
