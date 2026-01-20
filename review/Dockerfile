# Sandboxed environment for Ralph Wiggum review loop
# Isolates from host filesystem and restricts network to Anthropic API only

FROM python:3.11-slim

# Install Node.js for claude CLI
RUN apt-get update && apt-get install -y \
    curl \
    git \
    && curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y nodejs \
    && rm -rf /var/lib/apt/lists/*

# Install claude CLI
RUN npm install -g @anthropic-ai/claude-code

# Install Python dependencies (add what's needed for your repo)
RUN pip install --no-cache-dir \
    pyyaml \
    biopython \
    numpy \
    pandas

# Create non-root user
RUN useradd -m -s /bin/bash reviewer
USER reviewer

WORKDIR /repo

# Default command
CMD ["bash"]
