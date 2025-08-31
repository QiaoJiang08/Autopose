#!/bin/bash

# AutoDock Animation Generator - Installation Script
# This script installs all dependencies and sets up the application

set -e  # Exit on any error

echo "AutoDock Animation Generator - Installation"
echo "=========================================="

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.7 or higher."
    exit 1
fi

echo "✓ Python 3 found: $(python3 --version)"

# Check if pip is installed
if ! command -v pip3 &> /dev/null; then
    echo "❌ pip3 is not installed. Please install pip."
    exit 1
fi

echo "✓ pip3 found"

# Install Python dependencies
echo ""
echo "Installing Python dependencies..."
pip3 install -r requirements.txt

# Check if FFmpeg is installed
if ! command -v ffmpeg &> /dev/null; then
    echo ""
    echo "⚠️  FFmpeg is not installed. Attempting to install..."
    
    # Detect OS and install FFmpeg
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        if command -v brew &> /dev/null; then
            echo "Installing FFmpeg via Homebrew..."
            brew install ffmpeg
        else
            echo "❌ Homebrew not found. Please install FFmpeg manually:"
            echo "   brew install ffmpeg"
            exit 1
        fi
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Linux
        if command -v apt-get &> /dev/null; then
            echo "Installing FFmpeg via apt..."
            sudo apt-get update
            sudo apt-get install -y ffmpeg
        elif command -v yum &> /dev/null; then
            echo "Installing FFmpeg via yum..."
            sudo yum install -y ffmpeg
        else
            echo "❌ Package manager not found. Please install FFmpeg manually."
            exit 1
        fi
    else
        echo "❌ Unsupported OS. Please install FFmpeg manually."
        exit 1
    fi
fi

echo "✓ FFmpeg found: $(ffmpeg -version | head -n1)"

# Create necessary directories
echo ""
echo "Creating directories..."
mkdir -p uploads outputs

# Make scripts executable
chmod +x run.py test_setup.py

echo ""
echo "✅ Installation completed successfully!"
echo ""
echo "To test the installation:"
echo "  python3 test_setup.py"
echo ""
echo "To start the application:"
echo "  python3 run.py"
echo ""
echo "Or alternatively:"
echo "  python3 app.py"
echo ""
echo "Then open http://localhost:5000 in your web browser." 