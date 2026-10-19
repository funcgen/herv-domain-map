# T2T-CHM13 ERV-L Gag scope correction within v2

Baseline: published T2T v2, https://doi.org/10.5281/zenodo.23474864.
GitHub baseline: ed541ec4fc8f256392860dc04ce1673539ce896e.
This minor correction updates the T2T v2 annotation; it uses the existing human scope configuration and adds InterProScan JSON validation.

## Selection and evidence

The human scope already includes HERVL-int, ERVL-int, ERVL-B4-int and ERVL-E-int.
The configuration applies to both GRCh38 and T2T. This correction regenerates
T2T-specific data with that policy, without changing selector logic or thresholds.

All 18,628 published v2 selections are retained. The 67 additions comprise 42
ERVL-B4-int, 20 ERVL-int and 5 ERVL-E-int Gag-homologous records; 66 distinct peptide
sequences occur across those 67 loci. Two loci encode the same 186-aa peptide and
retain separate identifiers. The corrected catalog contains 18,695 domains in
7,275 loci; new-model Gag increases from 686 to 753. All scope winners use
GAG_ERVL_Hs. Coverage ranges from 11.1% to 34.2%; none reaches 40%.

Both sequence/domain scores must be >=95 bits for GAG_ERVL_Hs or >=130 for
GAG_MERVL, with independent-domain E <=1e-5. Existing sequence/conditional-domain
E <=1e-5 filters and deterministic best-per-locus/domain-type ranking are retained.
No minimum coverage is introduced. These are homology annotations, not demonstrated
complete or functional Gag proteins. All 66,000 exact-composition shuffles fail both
new-model score guards (1,000 per distinct peptide, seed 20261020); this is not a
biological FDR estimate. Reference similarity is not independent of model development.

The recovered native HMMER 3.4 rescan covers 130,923 ORFs and 237,349 reported-domain
rows, matching the complete published raw file. It is reverified, not rerun in this
correction. Current public selection reproduces initial v2 under the restricted
scope and 18,695 records under the expanded scope. All ORF/domain DNA translations
and native/UCSC coordinate mappings are validated. All earlier records are preserved.

## InterProScan and export validation

The full XML, TSV, JSON and GFF3 archive is rebuilt for the final catalog using
pinned InterProScan 5.75-106.0. All merged matches are checked against exact-peptide
source evidence. The 18,628 old identifiers reuse their published XML; 56 scope
records reuse completed GRCh38 member analyses and the remaining 11 represent ten
unique peptides with verified completed results in a recovered 13-peptide shared
batch. The historical batch is frozen for input/output hash verification; it is
not a new database scan claimed in this task. All 17 applications are accounted for.

PANTHER Fv1-family PTHR48195 recognizes 32 of 67 added records; none has a CDD match.
Full-catalog CDD regeneration yields the same 2,944 rows, with all 45,050 explicit
residue tokens matching. All 2,123 Env peptides and their Phobius results are unchanged.

The published v2 JSON and the pinned CONVERT export omit the final top-level array
and object closing tokens. The corrected export restores only those tokens after
validating every full record. XML and JSON include all 18,695 identifiers; TSV and
GFF3 report member matches for 16,444 identifiers. Missing TSV/GFF matches do not
mean missing XML/JSON sequence records.

`scripts/16_validate_interproscan_json.py` validates JSON sequences/MD5s/identifiers
against the complete source FASTA and optional XML. With `--repair-omitted-root`,
it can append the missing final results-array/root-object closing tokens only when
the resulting JSON is valid and represents the entire expected sequence set. Other
corruption is rejected without file mutation. The InterProScan wrapper calls this
validator before reporting success; it does not change search/model settings.

Example for an existing export:

```bash
python3 scripts/16_validate_interproscan_json.py domains.fasta.json \
  --fasta domains.fasta --xml domains.fasta.xml \
  --repair-omitted-root --report domains.fasta.json.validation.json
```

The validator was tested on the actual malformed published v2 JSON, the corrected
full export, and an internally corrupted JSON that must remain untouched. The shell
wrapper is syntax-checked; the full member-database searches are not rerun solely
to exercise its added post-processing validation call.

## Release contents and downstream effects

Replace six files: domain TSV, native domain BED, UCSC domain BED, domain AA FASTA,
domain coding-DNA FASTA and the complete InterProScan ZIP. Retain the other ten
published inputs/results, including CDD, Phobius, raw HMMER, models and classification.
Add scope/threshold files, individual candidate review, English report and evidence ZIP.
The original six affected v2 files are frozen for same-version auditability.

At >=40%, >=80%, >=90%, >=95% and >95%, selections and Gag/Pol/Env trio counts are
unchanged. Whole-catalog counts, family summaries, locus presence and complete HMM
alignments are regenerated. Python equivalents of affected article analyses were
executed; R itself was not executed. The current v2 coordinates and coding-DNA
orientation are already correct; no additional one-base correction is introduced.

This correction changes T2T documentation and InterProScan export validation. It does not
change the human or mouse scope configuration, profile models or score thresholds.
