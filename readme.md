
[//]: # (To compile this file as a pdf run the following command:
pandoc readme.md -o readme.pdf)

# Project desciption
Subject: Phase separation of the components of a binary fluid (AC)

Instructions/advices:

> The main paper to start with for the possible algorithms leading to a
practical scheme for solving the Cahn-Hilliard equation is Lee_CMS2014.pdf

> Jamet.pdf gives details about the physics of phase separation for a
binary fluid and the derivation of the Cahn-Hilliard equation you will
solve numerically.

> DJ.Eyre_1998.pdf is useful if you are interested in the way to
perform an operator splitting that will stabilize some of the schemes
you may write naively.

> TCGPhysica1995.pdf will possibly give you some ideas of diagnostics
that you can implement, once your solver will be running, to
characterize the patterns obtained in phase separation.
\

## Compiling report.tex 

To compile the LaTeX report in `document/` just run the command `pdflatex`&nbsp;`report.tex` from the `document/` directory.
As for dependencies, installing the package `texlive-latex-extra` should be enough.

Remember, for the `pdflatex` command the order of the options matters.
If you want to set an output directory, the `-output-directory` option should come *before* the name of the `.tex` file to be compiled.
This way you can compile the report directly from the root directory of the project into `document/`.

`pdflatex -output-directory document document/report.tex`

`-halt-on-error` is also a good option to avoid unnecessary error prompts. 
\
\

**Project done in collaboration by students Rayen Mahjoub, Gaspard Daumas and Wesley Rodrigues Machado. Supervised by professor Arnaud Couairon.**
