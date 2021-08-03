import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scKGAM",
    version="0.9",
    author="Elvis Cui",
    author_email="elviscuihan@g.ucla.edu",
    description="Single-cell Gene Expression Kinetics Generalized Additive Model.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ElvisCuiHan/scKGAM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)