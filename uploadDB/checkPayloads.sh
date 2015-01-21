#!/bin/bash

version=$1

echo "Payloads in $version:"
payloads="$(cmscond_list_iov -c sqlite_file:$version -a)"
for payload in $payloads; do
  echo
  cmscond_list_iov -c sqlite_file:$version -t $payload
done
