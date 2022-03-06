#!/usr/bin/env bash
main_tex_file="main.tex"
main_aux_file="main.aux"
latex_compiler="pdflatex"
# Exit when any command fails. A Tip from 
# https://intoli.com/blog/exit-on-errors-in-bash-scripts/ 
set -e
if [ $# -eq 0 ]; then
    cd tex
    noweave  -n -index -latex mainlitprog.nw > mainlitprog.tex
    $latex_compiler -interaction=nonstopmode -halt-on-error  -shell-escape $main_tex_file >> ../wtta.log 2>&1
    bibtex $main_aux_file                                                                 >> ../wtta.log 2>&1
    $latex_compiler -interaction=nonstopmode -halt-on-error  -shell-escape $main_tex_file >> ../wtta.log 2>&1
    $latex_compiler -interaction=nonstopmode -halt-on-error  -shell-escape $main_tex_file >> ../wtta.log 2>&1
    rm --force -f *-blx.bib *.aux *.bbl *.bcf *.blg *.ind *.idx *.ilg \
                  *.log *.out *.pbsdat *.prc *.pre *.run.xml *.tdo *.toc  *~
    printf " ( done )\n"
    printf "Finishing up ..."
    cd .. # Move back to home directory. 
    mv tex/main.pdf .                      # weaved document moved to home folder
    printf " ( done )\n"
fi
