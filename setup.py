from setuptools import setup, find_packages

setup(
    name="fullplasmid_demix",
    version="0.1.0",
    author="Rachel Weissman",
    author_email="your.email@example.com",
    description="A tool for demixing and assembling pooled plasmid sequences",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/rachelfweissman/fullplasmid_demix",
    packages=find_packages(),
    py_modules=["fullPlasmidSeq_demix_RFW"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "pandas",
        "biopython",
        "plotly",
        "openpyxl",
    ],
    entry_points={
        "console_scripts": [
            "fullplasmid-demix=fullPlasmidSeq_demix_RFW:main",
        ],
    },
) 