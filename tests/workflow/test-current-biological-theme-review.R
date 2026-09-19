#!/usr/bin/env Rscript
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
out<-file.path(root,"results/current_fsf_v1/manuscript/biological_themes")
sem<-file.path(root,"results/current_fsf_v1/manuscript/semantic_authority")
builder<-file.path(root,"scripts/current/biological_themes/build_current_biological_theme_review.R")
rd<-function(p)read.delim(p,check.names=FALSE,quote="",stringsAsFactors=FALSE)
sha<-function(p)sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
blank<-function(x)all(is.na(x)|!nzchar(x))
expected<-c(go_term_authority.tsv="608020c5f4c95263a34749b54b23a797c621184a27979dd59652506fcee1a218",
 kegg_term_authority.tsv="9652ad9bf7a89cf6cbb7f73da5cf99e8808241be284531bd47f85d21eef50671",
 current_significant_go_annotated.tsv="6bcf373e531b9dbfd8e1e89dc25ba77e2b2407c9c5dc3741408da281571786dd",
 current_significant_kegg_annotated.tsv="6bf9666965f84c7f183e75331f90cb878f2eabe0b8caf7adff6d520dfcb39735",
 biological_semantic_review_packet.tsv="c3b720079f405a8f9b20a4704435ae09e3a4e0887b559106bd4b85f996c42403")
stopifnot(all(vapply(file.path(sem,names(expected)),sha,"")==expected))

source_packet<-rd(file.path(sem,"biological_semantic_review_packet.tsv"))
trace<-rd(file.path(out,"biological_theme_evidence_trace.tsv"))
stopifnot(nrow(trace)==nrow(source_packet))
stopifnot(isTRUE(all.equal(trace[names(source_packet)],source_packet,check.attributes=FALSE)))
stopifnot(!anyDuplicated(trace$evidence_row_id))

mapping<-rd(file.path(out,"candidate_term_theme_mapping.tsv")
)
ga<-rd(file.path(sem,"go_term_authority.tsv"));ka<-rd(file.path(sem,"kegg_term_authority.tsv"))
stopifnot(all(mapping$term_id[mapping$source=="GO"]%in%ga$go_id))
km<-mapping[mapping$source=="KEGG",]
stopifnot(all(paste(km$term_id,ifelse(grepl("^K[0-9]",km$term_id),"KO","PATHWAY"))%in%
 paste(ka$kegg_id,ka$kegg_identifier_type)))
allowed_basis<-c("ONTOLOGY_RELATION","OFFICIAL_NAME_MATCH","KEGG_PATHWAY_IDENTITY",
 "HISTORICAL_HYPOTHESIS_SUPPORTED","CURRENT_EVIDENCE_DISCOVERY")
stopifnot(all(mapping$mapping_basis%in%allowed_basis))

files_with_author<-c("candidate_term_theme_mapping.tsv","candidate_theme_evidence.tsv",
 "within_condition_region_contrast.tsv","historical_theme_reconciliation.tsv",
 "biological_theme_author_review_packet.tsv")
for(f in files_with_author){x<-rd(file.path(out,f));stopifnot(blank(x$author_decision),blank(x$author_notes))}
ev<-rd(file.path(out,"candidate_theme_evidence.tsv"))
allowed_status<-c("STRONG_MULTI_SOURCE_SUPPORT","GO_SUPPORTED","KEGG_SUPPORTED","PARTIAL_SUPPORT",
 "WEAK_OR_BROAD_SUPPORT","NO_CURRENT_SUPPORT","REQUIRES_AUTHOR_REVIEW")
stopifnot(all(ev$evidence_status%in%allowed_status))
review<-rd(file.path(out,"biological_theme_author_review_packet.tsv"))
stopifnot(all(review$recommended_action_options==
 "RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER"))

regions<-c("Low Stability","Transitional","Stable","Highly Stable")
conditions<-c("DES","GAM","HT","LT","OSM","UV")
stopifnot(setequal(unique(trace$condition),conditions))
stopifnot(all(unique(trace$FSF_region)%in%regions))
for(f in c("candidate_term_theme_mapping.tsv","candidate_theme_evidence.tsv",
 "condition_region_biological_profile.tsv","within_condition_region_contrast.tsv")){
 x<-rd(file.path(out,f));cols<-intersect(names(x),c("FSF_region","region_A","region_B"))
 stopifnot(all(unlist(x[cols])%in%regions))
}
stopifnot(!any(tolower(unlist(ev[c("FSF_region")]))%in%c("instable","weakly stable","weakly_stable")))

defs<-rd(file.path(out,"candidate_theme_definitions.tsv"))
profile<-rd(file.path(out,"condition_region_biological_profile.tsv"))
summary_text<-readLines(file.path(out,"biological_theme_review_summary.md"),warn=FALSE)
authored<-tolower(c(defs$candidate_theme,profile$candidate_themes,profile$strongest_themes,
 profile$GO_only_themes,profile$KEGG_only_themes,profile$cross_source_themes,summary_text))
prohibited<-c("adaptive program","protective","fitness","causal","mechanistic","functionally required","robust program")
stopifnot(!any(vapply(prohibited,function(z)any(grepl(z,authored,fixed=TRUE)),logical(1))))

groups<-rd(file.path(out,"go_redundancy_groups.tsv"))
go<-rd(file.path(sem,"current_significant_go_annotated.tsv"))
members<-unlist(strsplit(groups$member_term_ids,";",fixed=TRUE))
stopifnot(setequal(members,unique(go$go_term)))

figure_expected<-c(
 "scripts/current/figures/build_current_manuscript_figures.R"="345a6077d2c23f8ba24c6a890422521a6359bba852120c56ff99e52c1cbb834c",
 "results/current_fsf_v1/manuscript/figures/Figure5.pdf"="ca8dbc79c55081cc4dc84bba87dcec54cb5e4e854f716aff674a77edebe4f8a7",
 "results/current_fsf_v1/manuscript/figures/Figure5.png"="e8405a4f8d9c421fe5a56143a0328b26087b187c7b2b1ecd994ecc8d544efa7f",
 "results/current_fsf_v1/manuscript/figures/Figure7.pdf"="c11babb4deed9bc75bcd09bcdc832d18481e309930907dce1db9d5e6a7ef6c17",
 "results/current_fsf_v1/manuscript/figures/Figure7.png"="df442be3cdfce8a5f23107afd6efeca9981d5fdb9d7e018478d914db632b0ef0")
stopifnot(all(vapply(names(figure_expected),sha,"")==figure_expected))
stopifnot(system2("git",c("diff","--quiet","HEAD","--","manuscript","docs"))==0L)

tmp<-tempfile("biological-theme-review-");dir.create(tmp)
log<-system2("Rscript",c("--vanilla",builder,tmp),stdout=TRUE,stderr=TRUE)
rebuilt<-rd(file.path(tmp,"biological_theme_output_hashes.tsv"))
frozen<-rd(file.path(out,"biological_theme_output_hashes.tsv"))
stopifnot(identical(frozen,rebuilt))
unlink(tmp,recursive=TRUE)
cat("current biological theme review tests: PASS\n")

