pandoc --filter pandoc-citeproc paper.md --bibliography=references.bib -s -o paper.pdf --verbose -M date="$(date "+%B %e, %Y")"
