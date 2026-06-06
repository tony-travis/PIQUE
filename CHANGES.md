# Changes

## SNP-level heterozygosity filter (`-y`)

A new option, `-y het_max`, has been added to `pique-input`. It removes SNPs where the
observed proportion of heterozygous individuals exceeds the supplied threshold. The value
must be a number between 0 and 1 (e.g. `-y 0.005` removes SNPs where more than 0.5% of
individuals are heterozygous).

The filter runs after the missingness filter (`-x`) but before the MAF filter (`-m`).
This ordering ensures that heterozygosity is estimated on the post-missingness SNP set,
and that the MAF threshold is then applied to the already-filtered data. When `-y` is
not supplied, the pipeline behaves exactly as before.

The filter works by running `plink2 --hardy` on the current dataset and parsing the
`O(HET_A1)` column directly in Perl. Two files are written alongside the existing output:

- `[prefix].hardy` -- per-SNP Hardy-Weinberg statistics from plink2
- `[prefix].sid` -- variant IDs removed by the filter (empty if none were removed)

These follow the same convention as the `.het` and `.hid` files produced by the
sample-level heterozygosity filter (`-z`).

## Changes to `qc()`

The internal `qc()` function was updated to accept an optional mode argument (`'geno'`
or `'maf'`), so that the missingness and MAF filters can be called separately when `-y`
is active. All existing call sites pass no mode argument and are unaffected.

## Tests

Two new test targets have been added to `test/Makefile`:

- `input_snphet` / `run_snphet` -- runs `pique-input` with `-y 0.005` and then EMMAX
- `input_snpzhet` / `run_snpzhet` -- runs both `-y 0.005` and `-z 0.999` together, to
  confirm that the SNP-level and sample-level filters coexist correctly

The threshold of 0.005 was chosen from inspection of the sativas413 test data. Most
SNPs in this highly inbred population have zero observed heterozygosity, but 649 SNPs
exceed 0.005, giving a verifiable subset to filter without disrupting downstream analysis.

## Known issue (not introduced here)

There is a pre-existing bug in the `het()` function where the sample exclusion file is
silently empty when `pique-input` is run without `-v`. This does not affect the new
`-y` filter. It will be fixed in a separate follow-up.
