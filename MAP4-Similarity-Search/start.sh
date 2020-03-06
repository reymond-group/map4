#!/bin/sh
cd Flask
exec gunicorn -b :5000 --access-logfile - --timeout 300 --error-logfile - wsgi:app