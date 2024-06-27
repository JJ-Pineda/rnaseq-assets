#!/bin/bash

# Note: this script must be run with "sudo ./ec2_ubuntu_startup.sh"

# Install docker
apt update
apt install -y --no-install-recommends docker.io
systemctl start docker

# Pull code
