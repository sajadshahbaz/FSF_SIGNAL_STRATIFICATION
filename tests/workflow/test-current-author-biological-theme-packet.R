#!/usr/bin/env Rscript
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
src<-file.path(root,"results/current_fsf_v1/manuscript/biological_themes")
out<-file.path(root,"results/current_fsf_v1/manuscript/biological_theme_decisions")
builder<-file.path(root,"scripts/current/biological_theme_decisions/build_author_biological_theme_decision_packet.R")
rd<-function(p)read.delim(p,check.names=FALSE,quote="",stringsAsFactors=FALSE)
sha<-function(p)sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
blank<-function(x)all(is.na(x)|!nzchar(x))
bhash<-c(candidate_theme_definitions.tsv="9dc65d51b39bc45d5885e000ae0e92ad7371e21732c17376527a5e1a2a450dfa",
 current_biological_evidence_inventory.tsv="bef472609c2790f1c6985493aa24149a65b7cbeb2b5f66ee7328bffefc600c58",
 biological_theme_evidence_trace.tsv="51be21df2d9484f83a4428ffe02a0d62e484415667eb102f5b785d9b5fe9ceee",
 candidate_term_theme_mapping.tsv="abab8afc65f4189d4c1b7220c495ba26ef793f84d57b1c67e67dfc76958e967b",
 go_redundancy_groups.tsv="6059c1e761f42c213222c5ea9540dc5a1ddee2155ba0c37905eee72b6faca639",
 candidate_theme_evidence.tsv="9c0a354ba0f98956ad844b35ee746e3dc54159cf39fb8d23c9d707816787c2a6",
 condition_region_biological_profile.tsv="62900e152b6f4b45529a5babe40396687dfa8a9e1fd6bcdc9531a9301d82a76c",
 within_condition_region_contrast.tsv="dad4cc1e90f58d7d731dac1d0535dd0eea6099776fe5e95429cfc7833be42f08",
 historical_theme_reconciliation.tsv="4a70702ddb57ae231247b4c0facd8f0ea31ee9bfc315410e0aaa1e7f09321958",
 biological_theme_author_review_packet.tsv="cb15514e28a556781f2a681b14f2c6f2ad56d48fd1a60fb0ead52b745c4e3f5f",
 biological_theme_review_summary.md="42335d4503e5d08e73f236fae208d941b6505403daa177fb6e76fb44090e0c9e")
stopifnot(all(vapply(file.path(src,names(bhash)),sha,"")==bhash))

defs<-rd(file.path(src,"candidate_theme_definitions.tsv"))
ev<-rd(file.path(src,"candidate_theme_evidence.tsv"))
mapping<-rd(file.path(src,"candidate_term_theme_mapping.tsv"))
trace<-rd(file.path(src,"biological_theme_evidence_trace.tsv"))
packet<-rd(file.path(out,"author_biological_theme_decision_packet.tsv"))
within<-rd(file.path(out,"within_condition_biological_review.tsv"))
global<-rd(file.path(out,"candidate_theme_global_review.tsv"))
stopifnot(nrow(global)==18L,setequal(global$candidate_theme,defs$candidate_theme))
required<-ev[ev$evidence_status=="REQUIRES_AUTHOR_REVIEW",]
required_key<-paste(required$condition,required$FSF_region,required$candidate_theme,sep="|")
stopifnot(all(required_key%in%packet$trace_key))
unit_key<-paste(ev$condition,ev$FSF_region,ev$candidate_theme,sep="|")
stopifnot(all(packet$trace_key%in%unit_key))
map_key<-paste(mapping$condition,mapping$FSF_region,mapping$candidate_theme,sep="|")
stopifnot(all(packet$trace_key%in%map_key))
trace_theme_key<-unique(unlist(lapply(seq_len(nrow(trace)),function(i){
 if(is.na(trace$candidate_themes[i])||!nzchar(trace$candidate_themes[i]))return(character())
 paste(trace$condition[i],trace$FSF_region[i],strsplit(trace$candidate_themes[i],"; ",fixed=TRUE)[[1]],sep="|")
})))
stopifnot(all(packet$trace_key%in%trace_theme_key))

for(x in list(packet,within,global))stopifnot(blank(x$author_decision),blank(x$author_notes))
options<-"RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER"
stopifnot(all(packet$recommended_action_options==options),all(global$recommended_action_options==options))
stopifnot(all(packet$review_priority%in%c("HIGH","MEDIUM","LOW")))
regions<-c("Low Stability","Transitional","Stable","Highly Stable")
stopifnot(all(packet$FSF_region%in%regions),all(within$region_A%in%regions),all(within$region_B%in%regions))

summary<-readLines(file.path(out,"author_biological_theme_review_summary.md"),warn=FALSE)
authored<-tolower(c(packet$candidate_theme,packet$within_condition_distinctiveness,
 packet$conflicting_or_ambiguous_evidence,within$review_interpretation_status,
 global$candidate_theme,global$region_specificity_summary,global$ambiguity_summary,summary))
prohibited<-c("adaptive program","protective","fitness","causal","mechanistic","functionally required","robust program")
stopifnot(!any(vapply(prohibited,function(z)any(grepl(z,authored,fixed=TRUE)),logical(1))))

fig<-c("scripts/current/figures/build_current_manuscript_figures.R"="345a6077d2c23f8ba24c6a890422521a6359bba852120c56ff99e52c1cbb834c",
 "results/current_fsf_v1/manuscript/figures/Figure5.pdf"="ca8dbc79c55081cc4dc84bba87dcec54cb5e4e854f716aff674a77edebe4f8a7",
 "results/current_fsf_v1/manuscript/figures/Figure5.png"="e8405a4f8d9c421fe5a56143a0328b26087b187c7b2b1ecd994ecc8d544efa7f",
 "results/current_fsf_v1/manuscript/figures/Figure7.pdf"="c11babb4deed9bc75bcd09bcdc832d18481e309930907dce1db9d5e6a7ef6c17",
 "results/current_fsf_v1/manuscript/figures/Figure7.png"="df442be3cdfce8a5f23107afd6efeca9981d5fdb9d7e018478d914db632b0ef0")
stopifnot(all(vapply(names(fig),sha,"")==fig))
protected<-c("results/current_fsf_v1/manuscript/semantic_authority",
 "results/current_fsf_v1/manuscript/biological_themes",
 "scripts/current/biological_themes","tests/workflow/test-current-biological-theme-review.R")
stopifnot(system2("git",c("diff","--quiet","HEAD","--",protected))==0L)

tmp<-tempfile("author-biological-theme-packet-");dir.create(tmp)
log<-system2("Rscript",c("--vanilla",builder,tmp),stdout=TRUE,stderr=TRUE)
rebuilt<-rd(file.path(tmp,"author_biological_theme_decision_hashes.tsv"))
frozen<-rd(file.path(out,"author_biological_theme_decision_hashes.tsv"))
stopifnot(identical(frozen,rebuilt))
unlink(tmp,recursive=TRUE)
cat("current author biological theme packet tests: PASS\n")

