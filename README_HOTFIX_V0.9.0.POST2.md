# CCID-Epigenome v0.9.0.post2 hotfix

This hotfix repairs parsing of the full GPL6883 platform table downloaded from
NCBI GEO using `view=data&form=text`.

The downloaded file is a GEO SOFT document.  Metadata lines precede a tabular
section delimited by:

```text
!platform_table_begin
...
!platform_table_end
```

Earlier code attempted to read the entire document as a TSV, causing:

```text
pandas.errors.ParserError: Expected 1 fields ..., saw 27
```

Version `0.9.0.post2` extracts only the delimited platform table before loading
it with pandas.  No candidate, score, threshold, or biological analysis is
changed.
