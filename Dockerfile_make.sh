#!/usr/bin/bash

echo "[Dockerfile_make] creating source archive"
git archive --format=tar.gz --output=maki-HEAD.tar.gz --prefix maki-latest/ HEAD

echo "[Dockerfile_make] building Docker image"
docker build -t $(whoami)/maki-2.0.3 .

echo "[Dockerfile_make] removing source archive"
rm maki-HEAD.tar.gz
