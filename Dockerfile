# Use Python 3.9 image
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first (for better caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application
COPY . .

# Build C++ engine
RUN cd cpp-engine && make || echo "C++ build skipped"

# Expose port
EXPOSE 8000

# Run the API
CMD ["python", "api.py"]
