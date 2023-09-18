#!/bin/bash

set -e

VERSION=`cat VERSION.txt`


docker push plantgenomics/entap:$VERSION
docker push plantgenomics/entap:latest