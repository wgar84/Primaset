#!/usr/bin/env bash

docker network create primaset_network
docker volume create --opt type=none --opt device=/home/guilherme/Primaset/pg --opt o=bind primaset_pg

docker build kernel -t primaset_kernel
