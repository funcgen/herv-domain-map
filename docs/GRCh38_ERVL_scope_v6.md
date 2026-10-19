# GRCh38 ERV-L Gag scope correction (v6 minor correction)

Baseline: the initial Zenodo v6 release, https://doi.org/10.5281/zenodo.23469632,
and GitHub `cbcdf518ee88ae214cacfbd920f3915d720103b2`. This minor correction
updates the v6 human annotation; the scope change below is the only pipeline change.

## Selection policy

The prior scan already searched all 128,427 released GRCh38 ORFs against all
316 profiles. Native HMMER 3.4 rescanning (`hmmscan --cpu 0 --noali -Z 316`;
33 query chunks) reproduced all 235,630 reported domain rows. The combined
HMM SHA256 is `27eb9f5f41cf6c0d054d879669a34e0ffdc9720950a8b4dd13a2dcf5bc0942b5`.
No model, profile threshold, or winner-ranking change is required.

Only the human subfamily allowlist is expanded from `HERVL-int` to
`HERVL-int`, `ERVL-int`, `ERVL-B4-int`, and `ERVL-E-int`. Mouse is unchanged.
Other families still go to `outside_scope_review` when appropriate. This is
an explicitly limited scope expansion based on reviewed candidates, not a
claim that every ERVL family or unrelated ERV group has been validated.

Both sequence and domain bit scores must be >=95 for `GAG_ERVL_Hs` or >=130
for `GAG_MERVL`; independent domain E <=1e-5. Legacy sequence and conditional
domain E <=1e-5 are retained. Coverage is rounded to three decimals, with no
new minimum-coverage filter. The best annotation per locus/domain type is
reselected from the complete raw scan, using the existing deterministic ranking.

## Result and evidence

The expanded selection contains 18,540 domains in 7,220 loci. All 18,473
published selected domains, peptide sequences, DNA sequences and BED rows
are preserved. Added Gag fragments: 42 `ERVL-B4-int`, 20 `ERVL-int`, 5
`ERVL-E-int`; 67 loci, 66 distinct peptides. There are 41 newly annotated
loci and 26 previously annotated loci gaining Gag; 25 of the latter have RT.
The new-model Gag count increases from 690 to 757; the total Gag count from
3,374 to 3,441. None of the 67 loci previously had a selected Gag annotation.

All candidates pass the frozen model guards (sequence bits 96.0-182.1,
domain bits 95.8-181.9; profile coverage 11.1%-32.2%). All are fragments:
45 have <20% coverage and 22 have >=20%; none reaches 40%. Composition
controls shuffled each of 66 unique peptides 1,000 times (seed 20261019),
preserving length and exact amino-acid composition. None of 66,000 controls
reached both 95-bit guards. This is a composition control, not an empirical
biological false-discovery estimate or evidence of protein function.

Local BLOSUM62 alignments to the published MuERV-L Gag reference
Y12713.1/CAA73250.1 support sequence homology. The reference backbone was
also used in human model development, so this comparison is not an independent
training/test validation. Additional member-database annotations and per-hit
evidence are provided in the accompanying release audit. Negative member
database results are recorded, not converted into positive functional evidence.

All 128,427 ORFs and all 18,540 selected domains translate from their released
coding-oriented DNA sequences. Genomic coordinates were regenerated through
the zero-based internal BED, with both strands verified. The public BED,
sequence-extraction, summary and alignment scripts reproduce the outputs.

## Affected resources

Replace the primary TSV, domain BED, domain AA FASTA, domain NT FASTA and
complete InterProScan ZIP. The fresh peptides are analyzed with the same pinned
InterProScan 5.75-106.0 member databases; unchanged peptides reuse their exact
published XML matches. XML, TSV, JSON and GFF3 are rebuilt for the full catalog.
The CDD-site table is retained if its regenerated rows remain identical.

All 2,114 Env identifiers and peptides are identical, so the published Phobius
output is retained. All coverage >=40%, >=80%, >=95% and >95% summary counts
are unchanged, and all three published XLSX workbook selections are verified
unchanged. Whole-catalog counts, family summaries and coverage distributions
must use the corrected selection. The 316 HMMs, classification table, raw
scan, internal regions, original ORFs and RepeatMasker input are retained.

This change does not regenerate or publish mouse or T2T datasets. Because the
configuration key is `human`, future T2T runs will use the expanded human
allowlist; their outputs require a separate release audit.

## Completed member-database audit

All 17 pinned InterProScan 5.75-106.0 applications completed for 66 new unique peptides. PANTHER Fv1 PTHR48195 recognizes 33 of the 67 added locus fragments (32 distinct peptides). CDD matches in the additions: 0. The complete CDD-site table was regenerated and is identical to the initial v6 release, so its original bytes are retained. These computational matches support homology; they do not demonstrate Gag or Fv1 function.
