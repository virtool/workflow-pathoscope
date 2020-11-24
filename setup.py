from pathlib import Path
from setuptools import setup, find_packages

AUTHORS = ["Ian Boyes", "Blake Smith"]

CLASSIFIERS = [
    "Topic :: Software Development :: Libraries",
    "Programming Language:: Python:: 3.9",
]

PACKAGES = find_packages(exclude="tests")


INSTALL_REQUIRES = [
    "virtool-workflow==0.1.1"
]

setup(
    name="vt-workflow-pathoscope",
    version="0.1.2",
    description="A workflow for detecting viruses in high-throughput sequencing (HTS) libraries.",
    long_description=Path("README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://github.com/virtool/workflow-pathoscope",
    author=", ".join(AUTHORS),
    license="MIT",
    platforms="linux",
    packages=PACKAGES,
    install_requires=INSTALL_REQUIRES,
    python_requires=">=3.9",
)