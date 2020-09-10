#!/bin/bash

# wait just to be safe
sleep 5
cd $HOME/WeFaceNano
python3 $HOME/WeFaceNano/manage.py makemigrations
python3 $HOME/WeFaceNano/manage.py migrate --noinput
python3 $HOME/WeFaceNano/manage.py runserver 0.0.0.0:8008
