#!/bin/bash

# Assuming only ONE jupyter server running
PORT=$(jupyter notebook list | grep -o ':[0-9]\+/?' | sed 's/://' | sed 's/\/?//')
TOKEN=$(jupyter notebook list | grep -o 'token=[a-z0-9]\+')
LOCAL_URL="http://localhost:$PORT/lab?$TOKEN"
echo $LOCAL_URL