#!/bin/bash

set -e
set -u

input="$1"
meta="$2"

function plot() {
    input="$1"
    meta="$2"
    terminal="$3"
    output="$4"
    cat <(echo "file = \"$input\";set terminal $terminal; set output '$output'") $meta | gnuplot
}

plot "$input" "$meta" "canvas" "comparison.html"
# plot "$input" "$meta" "png" "output.png"
