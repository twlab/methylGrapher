#!/bin/bash


tag="V0.1.2"

docker build --platform linux/amd64 -t wenjin27/methylgrapher:latest -t wenjin27/methylgrapher:$tag ./
# docker push wenjin27/methylgrapher:latest
# docker push wenjin27/methylgrapher:$tag


