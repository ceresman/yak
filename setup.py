from setuptools import setup, find_packages

setup(
    name="yak-genetics",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'biopython>=1.81',
        'pandas>=2.0.0',
        'numpy>=1.24.0',
        'pyvcf3>=1.0.0',
        'scipy>=1.11.0',
        'scikit-bio>=0.5.8'
    ],
    python_requires='>=3.8',
    author="Ceresman",
    description="Yak genetic breeding research pipeline",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
