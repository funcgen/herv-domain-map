# Mouse ERV-L Gag scope correction within v2

Baseline: mouse Zenodo record 23469639 (initial v2), https://doi.org/10.5281/zenodo.23469639.
GitHub baseline: ed541ec4fc8f256392860dc04ce1673539ce896e.
This minor correction updates the v2 mouse annotation; the scope change below is the only pipeline change.

## Scope and thresholds

The mouse list in `config/ERVL_primary_scope.json` includes `MERVL-int`,
`MERVL_2A-int`, `ERVL-int` and `ERVL-B4-int`. The human list is preserved.
The existing selector, deterministic ranking, 316 frozen HMM profiles and
score/E-value guards are unchanged. Both sequence and domain scores must
be >=130 bits for GAG_MERVL or >=95 bits for GAG_ERVL_Hs, with independent
domain E <=1e-5; legacy sequence and conditional-domain E <=1e-5 also apply.
No minimum model-coverage filter is added. Other families retain the review route.

## Result

The expanded selection adds three Gag-homologous regions (three unique peptides):
two ERVL-int and one ERVL-B4-int. All 48,798 original v2 selections are preserved.
The corrected catalog has 48,801 domains in 17,946 loci; the new-model Gag count
increases from 958 to 961. Profile spans are 21.0%, 65.7% and 75.4%.
Two additions pass >=40% coverage; none reaches >=80%, >=90% or >=95%.
Counts at >=40% change from 27,342 to 27,344 domains. Gag/Pol/Env trio counts
and all higher-coverage selections are unchanged.

Recovered HMMER 3.4 rescan output contains all 166,616 released ORFs and reproduces
all 613,939 published reported-domain rows. That completed scan was reverified,
not rerun for this correction. Running the current public selector on the complete
published raw data reproduces initial v2 with the restricted list and the corrected
48,801 records with the expanded list. All ORFs and selected-domain coding DNA
translate exactly; all previous BED and sequence records are preserved.

Each unique added peptide was reshuffled 1,000 times with exact composition retained
(seed 20261019). None of the 3,000 controls passes either new model's bit guards.
These controls do not estimate biological false-discovery rate. Similarity to
reference and reconstructed ancestral MuERV-L Gag is supporting evidence and is
not independent of earlier model development/validation.

Completed InterProScan 5.75-106.0 member runs were recovered and validated by
input/output hashes and exact peptide sequence. All 17 member applications are
accounted for. Two additions have PANTHER PTHR48195 Fv1-family matches; none has a
CDD match. Homology and model span do not establish functional or intact Gag.
The raw completed runs use an earlier 13-peptide batch; the mouse-specific exported
XML, TSV, JSON and GFF3 contain only the three additions.

## Release files and analyses

Replace the domain TSV, domain BED, domain peptide FASTA and domain coding-DNA
FASTA with the complete regenerated catalogs. Retain the other nine existing
mouse files, including all 16,975 CDD-site rows and the Phobius output for the
same 3,160 Env peptides. Add the reviewed scope JSON, frozen thresholds,
per-candidate review, English report and evidence ZIP. Original affected v2 bytes
are frozen in the evidence archive for the same-version correction.

Historical complete mouse InterProScan raw XML is unavailable. The audit provides
the new three-peptide subset, not an invented full-catalog raw archive. Historical
carried-forward CDD annotation versions remain unknown where the original evidence
recorded no version. Downstream summaries and complete HMM alignments are regenerated.
CDD joins use exact domain peptide identifiers, with >=95% and >95% distinguished.
Equivalent Python analyses were executed; R itself was not executed.

This correction changes only mouse scope and its documentation. It does not regenerate
T2T data or modify human scope. See the accompanying English mouse release report
for coordinates, per-candidate evidence, exact file manifest and reproduction steps.
