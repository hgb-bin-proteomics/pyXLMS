# Dockerfile for pyXLMS GUI
# author: Micha Birklbauer
# version: 1.0.1

FROM python:3.14

LABEL maintainer="micha.birklbauer@gmail.com"

RUN mkdir pyXLMS
COPY ./ pyXLMS/
WORKDIR pyXLMS

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install --no-cache-dir uv

WORKDIR gui

RUN uv sync --no-cache

CMD  ["uv", "run", "streamlit", "run", "streamlit_app.py"]
