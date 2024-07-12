#!/bin/bash

# Note: this script must be run with:
# "sudo ./ec2_startup_rnaseq_pre.sh <github access token> <AWS access key> <AWS secret key>"
# Fine-grained (read-only access) token for private rnaseq-assets repo
GITHUB_TOKEN=$1
AWS_ACCESS_KEY=$2
AWS_SECRET_KEY=$3

SECONDS=0

# First check for attached volume and mount if present
BLOCKS=$(lsblk -o NAME -n)
if [[ $BLOCKS == *"xvdf"* ]]
then
  # Open access to everybody
  chmod 777 /dev/xvdf

  # Build file system if needed (i.e. if we're using a brand new volume)
  FS_CHECK=$(file -s /dev/xvdf)
  if [[ $FS_CHECK != *"XFS filesystem"* ]]
  then
    mkfs -t xfs /dev/xvdf
  fi

  # Mount the volume
  mkdir ~/data
  mount /dev/xvdf ~/data
  chmod 777 ~/data
fi

# Install docker
cd
apt update
apt install -y --no-install-recommends git ca-certificates docker.io
systemctl start docker

# Pull code
git clone "https://oauth2:${GITHUB_TOKEN}@github.com/JJ-Pineda/rnaseq-assets.git"

# Build docker image and start up container
docker build --tag rnaseq_pre_image rnaseq-assets/docker_images/rnaseq_pre/.
docker run -d -t -v ~/data:/root/data -p 8080:8080 --name rnaseq_pre rnaseq_pre_image

# Transfer necessary files to docker container
docker cp rnaseq-assets/rnaseq_scripts/. rnaseq:/root/scripts/
mkdir .aws
echo "[default]\naws_access_key_id = ${AWS_ACCESS_KEY}\naws_secret_access_key = ${AWS_SECRET_KEY}" > .aws/credentials
docker cp .aws rnaseq_pre:/root/

duration=$SECONDS
echo "Set-up is complete after ~$(($duration / 60)) minutes"
