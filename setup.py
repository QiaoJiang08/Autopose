from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="autopose",
    version="2.1.0",
    author="QiaoJiang08",
    author_email="qiaojiang08@example.com",
    description="A comprehensive tool for generating MP4 video animations of protein docking results",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/QiaoJiang08/Autopose",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
    },
    entry_points={
        "console_scripts": [
            "autopose=app:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["templates/*", "static/*"],
    },
    keywords="protein docking animation autodock vina cluspro molecular visualization",
    project_urls={
        "Bug Reports": "https://github.com/QiaoJiang08/Autopose/issues",
        "Source": "https://github.com/QiaoJiang08/Autopose",
        "Documentation": "https://github.com/QiaoJiang08/Autopose#readme",
    },
)
