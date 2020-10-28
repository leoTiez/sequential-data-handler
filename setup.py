#!/usr/bin/python3
import setuptools

setuptools.setup(
    name="seqDataHandler",
    version="0.0.2",
    author="Leo Zeitler",
    author_email="leo.zeitler@i2bc.paris-saclay.fr",
    description="Data handler for managing bigwig and bed files and their respective data",
    url="git@github.com:leoTiez/seqDataHandler.git",
    packages=setuptools.find_packages(),
    install_requires=[
        'roman>=3.3',
        'numpy>=1.18',
        'scipy>=1.5',
        'pyBigWig>=0.3',
        'pybedtools>=0.8'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: CECILL License"
    ],
    python_requires='>=3.6',
)

