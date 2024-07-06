from setuptools import setup, find_packages

setup(
    name='lmm2d',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'libigl',
        'meshio'
    ],
    author='Gabriel Alkuino',
    author_email='gsalkuin@syr.edu',
    description='create LAMMPS data file from 2D triangular mesh',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/gsalkuin/lammps-magneto-mechanical-2d',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.0',
)
