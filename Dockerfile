# Container image that runs run_postprocessing.py.
#
# Build:  docker build -t mintpy-postprocessing .
# Run:    see "Run with Docker" in README.md.
#
# The base image ships micromamba, a small conda installer, so the image uses the
# same environment.yaml as a local setup instead of a separate requirements list.
FROM mambaorg/micromamba:2.0.5

# Install the dependencies first. This layer is cached and only rebuilt when
# environment.yaml changes, not on every code change.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/environment.yaml
RUN micromamba install --yes --name base --file /tmp/environment.yaml \
    && micromamba clean --all --yes

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER run_postprocessing.py ./
COPY --chown=$MAMBA_USER:$MAMBA_USER postprocessing ./postprocessing
COPY --chown=$MAMBA_USER:$MAMBA_USER tests ./tests

# The git commit the image is built from (set by the GitHub workflow). It is
# stored in each output file's metadata, because the image has no .git folder.
ARG GIT_COMMIT=""
ENV POSTPROCESSING_GIT_COMMIT=$GIT_COMMIT

# The base image's entrypoint activates the conda environment, then runs this.
# Arguments given to `docker run` are appended, e.g. `--run-dir /data/run`.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/app/run_postprocessing.py"]
