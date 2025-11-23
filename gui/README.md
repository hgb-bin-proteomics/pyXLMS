# pyXLMS Web Application

Additionally to the python package, pyXLMS features an easy-to-use web app that offers most of the functionality of the python package. The web app
can either be accessed using our public server or run locally or self-hosted.

## Public Web App

Our web app is publicly accessible for free via [hgb-bin-proteomics.github.io/pyXLMS-app](https://hgb-bin-proteomics.github.io/pyXLMS-app).
No data is saved and any analysis only persists as long as your session is active in your browser!

You can read more about our publicly hosted instance here: [github.com/hgb-bin-proteomics/pyXLMS-app](https://github.com/hgb-bin-proteomics/pyXLMS-app).

## Running Locally or Self-Hosting

Additionally, the web application can also be run locally or self-hosted. **This is the preferred way if you are analyzing large amounts of data or**
**sensitive data.** You can run the web app locally or self-hosted via:

- Cloning the [pyXLMS GitHub repository](https://github.com/hgb-bin-proteomics/pyXLMS):
  - ```bash
    git clone https://github.com/hgb-bin-proteomics/pyXLMS.git
    ```
  - ```bash
    cd pyXLMS
    ```
- **\[Recommended\]** Running the [streamlit](https://streamlit.io/) app via [uv](https://docs.astral.sh/uv/) locally:
  - ```bash
    cd gui
    ```
  - ```bash
    pip install uv
    ```
  - ```bash
    uv run -- streamlit run streamlit_app.py
    ```
- _or_ running the [streamlit](https://streamlit.io/) app natively locally:
  - ```bash
    cd gui
    ```
  - ```bash
    pip install -r requirements.txt
    ```
  - ```bash
    streamlit run streamlit_app.py
    ```
- _or_ running via [Docker](https://docs.docker.com/get-started/get-docker/) by building the image yourself:
  - ```bash
    docker build . -f Dockerfile -t pyxlms
    ```
  - ```bash
    docker run -p 8501:8501 pyxlms
    ```
- _or_ by pulling the provided image from Docker Hub:
  - ```bash
    docker run -p 8501:8501 michabirklbauer/pyxlms:latest
    ```
