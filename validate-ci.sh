#!/bin/bash
set -e
set -o pipefail

OUT_DIR=data

md5sum -c ${OUT_DIR}/md5s.txt