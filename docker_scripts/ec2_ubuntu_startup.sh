#!/bin/bash

# Note: this script must be run with "sudo ./ec2_ubuntu_startup.sh <github access token>"
# Fine-grained (read-only access) token for private rnaseq-assets repo
GITHUB_TOKEN=$1
AWS_ACCESS_KEY=$2
AWS_SECRET_KEY=$3

SECONDS=0

# Install docker
cd
apt update
apt install -y --no-install-recommends git ca-certificates docker.io
systemctl start docker

# Pull code
git clone "https://oauth2:${GITHUB_TOKEN}@github.com/JJ-Pineda/rnaseq-assets.git"

# Build docker image and start up container
docker build --tag rnaseq_image rnaseq-assets/docker_images/rnaseq/.
docker run -d -t -p 8080:8080 --name rnaseq rnaseq_image

# Transfer necessary files to docker container
docker cp rnaseq-assets/rnaseq_scripts/. rnaseq:/root/scripts/
mkdir .aws
echo "[default]\naws_access_key_id = ${AWS_ACCESS_KEY}\naws_secret_access_key = ${AWS_SECRET_KEY}" > .aws/credentials
docker cp .aws rnaseq:/root/

duration=$SECONDS
echo "Set-up is complete after ~$(($duration / 60)) minutes"
