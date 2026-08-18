#!/bin/bash

targets=(
  ubuntu24
  rocky9
)

for target in "${targets[@]}"
do
    docker build \
        -f recipes/docker/${target}.Dockerfile \
        -t maki:${target} .

    id=$(docker create maki:${target})

    mkdir -p artifacts/${target}

    docker cp \
        ${id}:/opt/maki/maki-${target} \
        artifacts/${target}/

    docker rm ${id}
done