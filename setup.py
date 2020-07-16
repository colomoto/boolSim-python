
from setuptools import setup, find_packages

NAME = 'boolsim'

setup(name=NAME,
    version='9999',
    author = "Loïc Paulevé",
    author_email = "loic.pauleve@labri.fr",
    url = "https://github.com/pauleve/boolSim-python",
    description = "Python interface to boolSim",
    install_requires = [
        "colomoto_jupyter",
        "pandas",
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="computational systems biology",

    py_modules = ["boolsim_setup", "boolsim"]
)

