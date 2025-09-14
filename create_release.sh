#!/bin/bash

# Autopose Release Script
# This script helps create a GitHub release for Autopose

set -e

VERSION="2.0.0"
TAG="v${VERSION}"
REPO="yourusername/autopose"

echo "🚀 Creating Autopose Release ${VERSION}"
echo "=================================="

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo "❌ Error: Not in a git repository"
    exit 1
fi

# Check if tag already exists
if git tag -l | grep -q "^${TAG}$"; then
    echo "⚠️  Tag ${TAG} already exists"
    read -p "Do you want to delete and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        git tag -d ${TAG}
        git push origin :refs/tags/${TAG}
    else
        echo "❌ Aborting release"
        exit 1
    fi
fi

# Check if all files are committed
if ! git diff-index --quiet HEAD --; then
    echo "❌ Error: Uncommitted changes detected"
    echo "Please commit all changes before creating a release"
    git status
    exit 1
fi

# Create the tag
echo "📝 Creating tag ${TAG}..."
git tag -a ${TAG} -m "Release Autopose ${VERSION} - ClusPro Extension

🎉 Major Release: ClusPro Protein-Protein Docking Support

NEW FEATURES:
- ✅ ClusPro protein-protein docking animations
- ✅ Interface contact analysis and visualization
- ✅ Multiple model animation with clustering
- ✅ Enhanced overlay information display
- ✅ Professional homodimer visualization

IMPROVEMENTS:
- ✅ Improved AutoDock4 overlay parsing
- ✅ Enhanced Vina multi-pose handling
- ✅ Better error handling and validation
- ✅ Optimized memory management
- ✅ Comprehensive documentation

BUG FIXES:
- ✅ Fixed AutoDock4 overlay display issues
- ✅ Resolved Vina pose positioning problems
- ✅ Corrected ClusPro frame generation
- ✅ Fixed file upload validation

"

# Push the tag
echo "📤 Pushing tag to GitHub..."
git push origin ${TAG}

echo "✅ Tag ${TAG} created and pushed successfully!"
echo ""
echo "🎯 Next Steps:"
echo "1. Go to https://github.com/${REPO}/releases"
echo "2. Click 'Create a new release'"
echo "3. Select tag '${TAG}'"
echo "4. Use the content from RELEASE_NOTES.md for the description"
echo "5. Publish the release"
echo ""
echo "📋 Release Checklist:"
echo "- [ ] Tag created and pushed"
echo "- [ ] GitHub release created"
echo "- [ ] Release notes published"
echo "- [ ] Community notified"
echo ""
echo "🎉 Autopose ${VERSION} is ready for release!"
