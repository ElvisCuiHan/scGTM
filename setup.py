import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scBCM",
    version="0.3",
    author="Elvis Cui",
    author_email="elviscuihan@g.ucla.edu",
    description="Bell-shape Curve Model.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ElvisCuiHan/scBCM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

