#!/bin/bash

# Note: this script must be run with:
# "sudo ./ec2_startup_rstudio.sh <github access token> <AWS access key> <AWS secret key>"
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
  mkdir /home/ubuntu/javier
  mount /dev/xvdf /home/ubuntu/javier
  chmod 777 /home/ubuntu/javier
fi

# Install docker
cd
apt update
apt install -y --no-install-recommends git ca-certificates docker.io
systemctl start docker

# Pull code
git clone "https://oauth2:${GITHUB_TOKEN}@github.com/JJ-Pineda/rnaseq-assets.git"

# Build docker image and start up container
docker build --tag rstudio_image rnaseq-assets/docker_images/rstudio/.
docker run -d -t -v /home/ubuntu/javier:/home/rstudio/javier -p 8787:8787 --name rstudio rstudio_image

# Give the container time to boot up before moving onward
sleep 10

# Transfer necessary files to docker container
mkdir /home/ubuntu/.aws
echo "[default]\naws_access_key_id = ${AWS_ACCESS_KEY}\naws_secret_access_key = ${AWS_SECRET_KEY}" > /home/ubuntu/.aws/credentials
docker cp /home/ubuntu/.aws rstudio:/home/rstudio/

duration=$SECONDS
echo "Set-up is complete after ~$(($duration / 60)) minutes"
