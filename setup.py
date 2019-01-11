import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cosmos",
    version="0.0.1",
    author="Peter Lustig",
    author_email="peter.lustig@physik.lmu.de",
    description="source matching with cosmos field sources and further analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pslustig/cosmos.git",
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'astropy', 'matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
