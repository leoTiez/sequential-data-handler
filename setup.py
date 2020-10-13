#!/usr/bin/python3
import setuptools

setuptools.setup(
    name="seqDataHandler",
    version="0.0.1",
    author="Leo Zeitler",
    author_email="leo.zeitler@i2bc.paris-saclay.fr",
    description="Data handler for managing bigwig and bed files and their respective data",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: CECILL License"
    ],
    python_requires='>=3.6',
)

