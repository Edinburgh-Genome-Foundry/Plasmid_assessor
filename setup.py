from setuptools import setup, find_packages

version = {}
with open("plasmid_assessor/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="plasmid_assessor",
    version=version["__version__"],
    author="Peter Vegh",
    author_email="egf-software@ed.ac.uk",
    description="Plasmid assessment for Golden Gate cloning",
    long_description=open("pypi-readme.rst").read(),
    long_description_content_type="text/x-rst",
    license="MIT",
    keywords="biology dna",
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    install_requires=[],
)