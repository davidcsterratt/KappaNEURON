import setuptools

VERSION = '0.3.2'
POSTVERSION = '' # rc2

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='KappaNEURON',
    version=VERSION + POSTVERSION,
    author="David C Sterratt",
    author_email="david.c.sterratt@ed.ac.uk",
    description="KappaNEURON integrates the SpatialKappa simulator with NEURON to allow rule-based simulations of molecular systems embedded in neurons.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/davidcsterratt/KappaNEURON",
    packages=setuptools.find_packages(),
    install_requires=['SpatialKappa >= 2.1.5', 'SciPy', 'matplotlib < 3.0'],
    package_data={'KappaNEURON': ['tests/*.ka', 'tests/*.mod', 'tests/*.inc', 'README.md']},
    classifiers=(
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
    python_requires='<3')
