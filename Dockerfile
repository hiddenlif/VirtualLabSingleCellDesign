FROM python:3.13-slim
RUN apt-get update && apt-get install -y --no-install-recommends procps util-linux && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.3 scipy==1.15.2
RUN useradd -m -u 1000 solver && mkdir -p /app
WORKDIR /app
COPY app/ /app/
RUN chown -R solver:solver /app
USER solver
CMD ["sleep", "infinity"]
