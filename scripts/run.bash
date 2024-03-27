#!/bin/bash

# This script is a master wrapper script for this project. It can be modified
# to port this project to other environments.

set -e

cmd=$1
shift 1

do_mafft() {
    exec mafft "$@"
}

do_rscript() {
    exec Rscript --vanilla "$@"
}

case $cmd in
  mafft)
    do_mafft "$@"
    ;;
  Rscript)
    do_rscript "$@"
    ;;
  *)
    echo This script does not know how to "$cmd".
    exit 1
    ;;
esac
