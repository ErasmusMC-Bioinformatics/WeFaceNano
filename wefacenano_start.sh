#!/bin/bash

# wait just to be safe
sleep 5
cd /root/
python3 /WeFaceNano/manage.py makemigrations
python3 /WeFaceNano/manage.py migrate --noinput
python3 /WeFaceNano/manage.py runserver 0.0.0.0:8008