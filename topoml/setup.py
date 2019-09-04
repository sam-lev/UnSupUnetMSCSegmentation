"""
      Setup script for topoml
"""
from setuptools import setup
import re


extra_args = {}


def get_property(prop, project):
    """
        Helper function for retrieving properties from a project's
        __init__.py file
        @In, prop, string representing the property to be retrieved
        @In, project, string representing the project from which we will
        retrieve the property
        @Out, string, the value of the found property
    """
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)


VERSION = get_property("__version__", "topoml")

# Consult here: https://packaging.python.org/tutorials/distributing-packages/
setup(
    name="topoml",
    packages=["topoml"],
    version=VERSION,
    description="A library for performing topology-based machine learning",
    long_description="TODO",
    author="Sam Leventhal, Steve Petruzza, Attila Gyulassy, Torin McDonald, and Dan Maljovec",
    author_email="maljovec002@gmail.com",
    license="BSD",
    url="https://bitbucket.org/cedmav/topoml/src/master/",
    keywords=[""],
    # Consult here: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: C++",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    # TODO: add versions where applicable
    install_requires=[
        "scipy",
        "numpy",
        "matplotlib",
        "sklearn",
        "scikit-image",
        "networkx",
    ],
    python_requires=">=2.7, <4",
)
