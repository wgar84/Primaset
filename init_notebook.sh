#!/usr/bin/env bash

docker run \
    -it \
    --rm \
    --network primaset_network \
    -p 62000:8888 \
    -v "$(pwd)":/home/jovyan/app \
    --name 'primaset_kernel' \
    primaset_kernel \
    jupyter notebook \
    --no-browser \
    --ip=0.0.0.0
