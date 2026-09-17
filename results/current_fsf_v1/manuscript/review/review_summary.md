# FSF biological interpretation review package

This directory contains evidence views and historical-rule candidates for author review. It contains no accepted biological theme authority and no final manuscript prose.

- Current main-class gene sets reviewed: 40
- Historical interpretation rules inventoried: 52
- GO historical mappings with current matching evidence: 11
- KEGG historical mappings with current matching evidence: 10
- Profile status counts: AMBIGUOUS_REQUIRES_AUTHOR_REVIEW=13; NO_CURRENT_EVIDENCE=3; PARTIALLY_SUPPORTED=6; SUPPORTED_BUT_LABEL_UPDATE=4; SUPPORTED_UNCHANGED=14

## Historical mapping evidence

- GO mappings with matching current terms: Translation / ribosome; Protein targeting / secretion; Genome maintenance / DNA repair; Chromosome / nuclear organization; Cell cycle / cell division; Development / morphogenesis; Neural / projection organization; Movement / behavior / taxis; Cell adhesion / junction; Metabolism / oxidation-reduction; Broad regulatory / cellular process
- KEGG mappings with matching current terms: Translation / ribosome; Genome maintenance / DNA repair; Cell cycle regulation; Protein turnover / proteostasis; Energy / core metabolism; Stress signaling / signal transduction; Cell adhesion / extracellular interaction; Neural signaling; Endocrine / physiological regulation; Broad proliferation / regulatory pathway
- Historical mappings without current matches: Cytoskeleton / cellular remodeling
- Broad/default and combined-profile mappings remain manual biological judgments even where their source terms are present.

## Condition review

- DES, GAM, and LT contain Low Stability rows requiring the locked label migration; their non-UV memberships are unchanged.
- HT and OSM are entirely Highly Stable in the region authority; biological wording still requires author review.
- Statistical count differences versus historical curated tables are recorded per gene set in combined_profile_review.tsv.

## UV review gate

UV is Transitional-dominated (6,697 of 15,187); its historical profiles remain REQUIRES_AUTHOR_REVIEW because current class membership differs materially.

## Figure 5 candidates

Candidate mappings are the GO/KEGG rules marked SUPPORTED_UNCHANGED in the review tables. Presence of matching enriched terms does not constitute author acceptance of a theme.

## Figure 7 candidate structure

Current condition architecture -> current signal class -> GO/KEGG statistical evidence -> author-reviewed theme/profile. Causal or adaptive language is outside the automatic evidence layer.

Author decisions and notes are intentionally blank.
