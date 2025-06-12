# Use a base Python image
# We're using a slightly older but stable Python 3.9 on Debian Buster
# as it's often more compatible for building native extensions like Infomap.
FROM python:3.9-slim-buster

# Set the working directory inside the container
WORKDIR /app

# Install system dependencies required for Infomap and building Python packages
# `build-essential` provides compilers (gcc, g++)
# `cmake` is often used by build systems
# `git` might be needed if infomap is cloned from source
# `libxml2-dev`, `zlib1g-dev` are common dependencies for various Python packages or native libraries.
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libxml2-dev \
    zlib1g-dev \
    # Clean up apt caches to reduce image size
    && rm -rf /var/lib/apt/lists/*

# Install infomap Python library
# This assumes the `infomap` pip package successfully builds its C++ dependency.
# This is the most common scenario on Linux-based environments.
RUN pip install infomap

# Copy your requirements.txt and install Python dependencies
# Using --no-cache-dir to reduce image size
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy your application code into the container
# This includes app.py, backend.py, and the example_data folder
COPY . .

# Expose the port Streamlit runs on (default is 8501)
EXPOSE 8501

# Command to run your Streamlit application when the container starts
# --server.port=8501: Tells Streamlit to listen on this port
# --server.enableCORS=false & --server.enableXsrfProtection=false: Often necessary for deployment platforms
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.enableCORS=false", "--server.enableXsrfProtection=false"]