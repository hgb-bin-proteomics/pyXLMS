# Dockerfile for pyXLMS GUI
# author: Micha Birklbauer
# version: 1.0.0

FROM python:3.12

LABEL maintainer="micha.birklbauer@gmail.com"

RUN mkdir pyXLMS
COPY ./ pyXLMS/
WORKDIR pyXLMS

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install --no-cache-dir .[gui]

WORKDIR gui

CMD  ["streamlit", "run", "streamlit_app.py"]
