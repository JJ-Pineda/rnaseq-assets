#!/bin/bash

# Note: this script must be run with:
# "sudo ./ec2_startup_rstudio.sh <github access token> <AWS access key> <AWS secret key>"
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
docker build --tag jupyter_image rnaseq-assets/docker_images/jupyter/.
docker run -d -t -p 8888:8888 --name jupyter jupyter_image

# Transfer necessary files to docker container
mkdir .aws
echo "[default]\naws_access_key_id = ${AWS_ACCESS_KEY}\naws_secret_access_key = ${AWS_SECRET_KEY}" > .aws/credentials
docker cp .aws jupyter:/home/jovyan/

# Get localhost url to Jupyter notebook
TOKEN=$(docker exec jupyter jupyter notebook list | grep -o 'token=[a-z0-9]\+')
LOCAL_URL="http://localhost:8888/lab?$TOKEN"

duration=$SECONDS
echo "Set-up is complete after ~$(($duration / 60)) minutes."
echo "Set up an SSH tunnel locally with \"ssh -i <pem file path> -L 8888:localhost:8888 ubuntu@<ec2-instance>\""
echo "Then use the following link to access Jupyter:"
echo "$LOCAL_URL"
