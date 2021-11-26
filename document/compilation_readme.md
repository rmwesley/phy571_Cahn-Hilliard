To compile the document just run `pdflatex report.tex`.
As for dependencies, installing the package `texlive-latex-extra` should be enough.

Remember, for `pdflatex`, order is important.
I lost some time to compile the report from the root directory of the project.

`pdflatex -halt-on-error -output-directory document document/report.tex`


Some minimalist evangelism: use `zathura` as a pdf reader!

`zathura document/report.pdf`

Also, you can compile this `.md` file with pandoc:

`pandoc compilation_readme.md -o temp.pdf`
