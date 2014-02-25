#!/bin/bash

set -e
set -u

input="$1"

function plot() {
    input="$1"
    terminal="$2"
    output="$3"
    cat <(echo "file = \"$input\";set terminal $terminal; set output '$output'") gnuplot.plt | gnuplot
    # open "$output"
}

plot "$input" "canvas" "comparison.html"
# plot "$input" "png" "output.png"
