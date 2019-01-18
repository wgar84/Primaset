#!/usr/bin/env bash

docker run \
    -it \
    --rm \
    --network prima_network \
    -v primaset
    --name 'primaset_postgres' \
    postgres:10.6
