#!/bin/bash

set -e

VERSION=`cat VERSION.txt`


docker build -t plantgenomics/entap:$VERSION .
docker build -t plantgenomics/entap:latest .