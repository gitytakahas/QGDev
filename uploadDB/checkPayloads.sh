#!/bin/bash

version=$1

echo "Payloads in QGL_$version.db:"
payloads="$(cmscond_list_iov -c sqlite_file:QGL_$version.db -a)"
for payload in $payloads; do
  echo
  cmscond_list_iov -c sqlite_file:QGL_$version.db -t $payload
done
