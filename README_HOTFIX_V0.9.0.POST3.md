# CCID-Epigenome v0.9.0.post3 hotfix

The ENCODE rE2G files can contain millions of element-gene-link rows. v0.9.0.post2 loaded both files completely into pandas, which can trigger the Linux/WSL out-of-memory killer and terminate the process with only `Killed`.

v0.9.0.post3 streams each gzip file line by line and retains only intervals overlapping the frozen 16-CpG panel. It also reads the actual ENCODE header fields (`TargetGene` and `Score`) rather than relying on generic BED column positions.

After applying the patch, rerun `04k_assess_regulatory_links_and_editability.py` with the same command. Delete any incomplete `fibroblast_re2g_candidate_links_v090.csv` first. Guide-design outputs are unchanged.
