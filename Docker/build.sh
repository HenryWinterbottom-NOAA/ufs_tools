#!/bin/sh

docker build -f ubuntu20.04-miniconda_ufs_pyutils.ufs_tools.docker --tag ubuntu20.04-miniconda_ufs_pyutils.ufs_tools:latest .
docker tag ubuntu20.04-miniconda_ufs_pyutils.ufs_tools:latest noaaufsrnr/ubuntu20.04-miniconda_ufs_pyutils.ufs_tools:latest
docker push noaaufsrnr/ubuntu20.04-miniconda_ufs_pyutils.ufs_tools:latest

