
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

## Compiling report.tex 

To compile the LaTeX report in document/ just run the command `pdflatex report.tex` from the document/ directory after getting all correct dependencies.
Specifically, you can install texlive-latex-extra from the texlive distribution to get the basic packages.
Then, to get Polytechnique's styling and graphical packages you need to also install the polytechnique package of the group typographix.
It isn't enough anymore for the newer versions of report.tex to compile it using just packages readily available on latex distributions like texlive.
You need to be able to use Polytechnique's styling packages.
Install them system-wide or decompress the package's files on the document/ folder.

Remember, for the `pdflatex` command the order of the options matters.
If you want to set an output directory, the `-output-directory` option should come *before* the name of the `.tex` file to be compiled.
This way you can compile the report directly from the root directory of the project into `document/`.

`pdflatex -output-directory document document/report.tex`

`-halt-on-error` is also a good option to avoid unnecessary error prompts.

## Example animation:

https://user-images.githubusercontent.com/64560431/159530799-3d379950-dbc6-48ae-a0bc-7ca31ae3f2c5.mp4

<br>

### Project done in collaboration by students Rayen Mahjoub, Gaspard Daumas and Wesley Rodrigues Machado. Supervised by professor Arnaud Couairon.

