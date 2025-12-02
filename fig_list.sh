#!/usr/bin/env bash

set -e

src="./images"
dst="./images_keep"
list="fig_list.txt"

mkdir -p "$dst"

while IFS= read -r f; do
    [ -z "$f" ] && continue
    cp "$src/$f" "$dst/"
done < "$list"

